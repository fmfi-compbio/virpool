import csv
import math
import argh
import numpy as np
from scipy.special import softmax
import pathos.multiprocessing as mp

import helpers
from helpers import load_alignments, apply_to_cigartuples, SIMPLE_ALPHABET, \
    get_genome_size_from_alignments, load_masked_positions


class PosteriorCounter:
    EPS = 1e-300

    def __init__(self, variants, subst_rate, alphabet,
                 clipping_start=0, clipping_end=math.inf,
                 pvalues=None, pvalue_threshold=None,
                 masked_positions=None):
        # values in `variants` are assumed to be normalised
        self.variants = variants
        self.subst_rate = subst_rate
        self.ces = [0 for _ in range(self.variants.shape[0])]
        self.L2N = {l: n for n, l in enumerate(alphabet)}
        self.clipping_start = clipping_start
        self.clipping_end = clipping_end
        self.pvalues = pvalues
        self.pvalue_theshold = np.log(pvalue_threshold) if pvalue_threshold is not None else None
        self.masked_positions = set(masked_positions) if masked_positions is not None else set()

    def _count_alongside_cigar(self, op, l, r, q, al):
        # this model ignores symbols outside the alphabet
        query = al.query_sequence
        if op == 0 or op == 7 or op == 8:
            for k in range(l):
                base = query[q + k]
                if base not in self.L2N:
                    continue
                qln = self.L2N[base]
                ref_pos = r + k

                if ref_pos < self.clipping_start or ref_pos >= self.clipping_end:
                    # ignore bases outside of [clipping_start, clipping_end)
                    continue

                if self.pvalues is not None and self.pvalues[ref_pos, qln] > self.pvalue_theshold:
                    continue

                if ref_pos in self.masked_positions:
                    continue

                for cnum in range(self.variants.shape[0]):
                    r_a = self.variants[cnum][ref_pos][qln]
                    r_a_prime = r_a
                    # q_a = (1 - self.subst_rate) * r_a_prime + \
                    #       self.subst_rate / 3 * (1 - r_a_prime)
                    total = sum(self.variants[cnum][ref_pos])
                    q_a = (1 - self.subst_rate) * r_a_prime + \
                          self.subst_rate / 3 * (total - r_a_prime)

                    ce = -math.log(self.EPS + q_a)
                    self.ces[cnum] += ce

    def dump_cross_entropies(self):
        return self.ces

    def dump_posteriors(self):
        return list(softmax([-x for x in self.ces]))


def eval_cross_entropies(al, variants, subst_rate, clipping_start, clipping_end, pvalues=None, pvalue_threshold=None, masked_positions=None):
    counter = PosteriorCounter(variants, subst_rate, SIMPLE_ALPHABET, clipping_start, clipping_end,
                               pvalues=pvalues, pvalue_threshold=pvalue_threshold,
                               masked_positions=masked_positions)
    apply_to_cigartuples(counter._count_alongside_cigar, al)
    result = [al.query_name] + list(counter.dump_cross_entropies())
    return result


def merge_paired_reads_ces(ces):
    data = dict()
    for line in ces:
        id, nums = line[0], line[1:]
        if id not in data:
            data[id] = nums
        else:
            # summing log-posteriors for paired reads
            data[id] = [o + n for o, n in zip(data[id], nums)]
    result = [[id] + nums for id, nums in data.items()]
    return result


class SimplifiedAlignment:  # is needed to use pathos.multiprocessing library
    def __init__(self, alignment):
        self.query_name = alignment.query_name
        self.reference_start = alignment.reference_start
        self.cigartuples = alignment.cigartuples
        self.query_sequence = alignment.query_sequence


def dump_posteriors(variant_names, posteriors, f):
    writer = csv.writer(f, delimiter="\t")
    header = ["read_id"] + variant_names
    writer.writerow(header)
    writer.writerows(posteriors)


def dump_posteriors_as_np(rows):
    identifiers = []
    posteriors = []
    for row in rows:
        posterior = list(map(float, row[1:]))
        identifiers.append(row[0])
        posteriors.append(posterior)
    posteriors_np = np.array(posteriors, dtype=np.float64)
    return identifiers, posteriors_np


def evaluate_ces_multiprocessing(
        alignments,
        threads_count,
        variants,
        subst_rate,
        clipping_start,
        clipping_end,
        pvalues=None,
        pvalue_threshold=None,
        masked_positions=None
    ):
    als = map(lambda al: SimplifiedAlignment(al), alignments)
    with mp.Pool(threads_count) as p:
        posteriors = p.map(lambda al: eval_cross_entropies(
            al,
            variants,
            subst_rate,
            clipping_start,
            clipping_end,
            pvalues=pvalues,
            pvalue_threshold=pvalue_threshold,
            masked_positions=masked_positions
        ),
        als)
    return posteriors


def posteriors_simple(alignments, variant_table, subst_rate: float, threads: int,
                      clipping_start: int, clipping_end: int, pvalues=None, pvalue_threshold=None,
                      masked_positions=None):
    ces = evaluate_ces_multiprocessing(
        alignments,
        threads,
        variant_table,
        subst_rate,
        clipping_start,
        clipping_end,
        pvalues=pvalues,
        pvalue_threshold=pvalue_threshold,
        masked_positions=masked_positions
    )
    merged_ces = merge_paired_reads_ces(ces)
    merged_posteriors = map(lambda x: [x[0]] + list(softmax([-a for a in x[1:]])), merged_ces)
    return merged_posteriors


def load_pvalues(f):
    lines = []
    reader = csv.reader(f, delimiter="\t")
    maxpos = 0
    for row in reader:
        pos, base, _, logpval = int(row[0]), helpers.L2N[row[1]], int(row[2]), float(row[3])
        lines.append((pos, base, logpval))
        maxpos = max(maxpos, pos)
    result = np.zeros((maxpos+1, len(helpers.SIMPLE_ALPHABET)))
    for pos, base, logpval in lines:
        result[pos, base] = logpval
    return result


def apply_bonferroni_correction(pvalues):
    pvalues += np.log(pvalues.shape[0]) + np.log(pvalues.shape[1])


@argh.arg("subst-rate", type=float)
@argh.arg('-t', '--threads', type=int)
@argh.arg("-cs", "--clipping-start", type=int, help="Number of ignored bases from the start of the ref. sequence")
@argh.arg("-ce", "--clipping-end", type=int, help="Number of ignored bases from the end of the ref. sequence")
@argh.arg("-m", "--masking-vcf", help="VCF file with masked positions")
@argh.arg("--percentile", help="Punish gaps with coverage < `percentile` quantile for the given variant")
def main(variants: str,
         alignments: str,
         output: str,
         subst_rate: float,
         threads: int = 1,
         clipping_start: int = 0,
         clipping_end: int = 0,
         pvalues: str = None,
         pvalue_threshold: float = 0.001,
         masking_vcf: str = None,
         percentile: float = 0.0
         ):
    variant_names, table = helpers.load_variants_with_simple_alphabet_punish_gaps(variants, percentile)

    # convert `clipping_end` from number of clipped bps at the end to the position of the first removed basepair
    genome_size = get_genome_size_from_alignments(alignments)
    clipping_end = genome_size - clipping_end

    alignments = load_alignments(alignments)

    if pvalues is not None:
        with open(pvalues) as f:
            pvalues = load_pvalues(f)
        apply_bonferroni_correction(pvalues)

    if masking_vcf is not None:
        with open(masking_vcf) as f:
            masked_positions = load_masked_positions(f)
    else:
        masked_positions = None

    merged_posteriors = posteriors_simple(alignments=alignments,
                                          variant_table=table,
                                          subst_rate=subst_rate,
                                          threads=threads,
                                          clipping_start=clipping_start,
                                          clipping_end=clipping_end,
                                          pvalues=pvalues,
                                          pvalue_threshold=pvalue_threshold,
                                          masked_positions=masked_positions)

    with open(output, "w") as f:
        dump_posteriors(variant_names, merged_posteriors, f)


if __name__ == "__main__":
    argh.dispatch_command(main)
