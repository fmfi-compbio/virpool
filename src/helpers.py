import csv
import itertools
import logging
import re
from collections import defaultdict
from typing import TextIO
import numpy as np
import pysam as pysam

FASTA_LABEL_SYMBOL = ">"
SIMPLE_ALPHABET = "ACGT"
L2N = {l: n for n, l in enumerate(SIMPLE_ALPHABET)}
FULL_ALPHABET = "ACGTN-"


def load_fasta(f: TextIO):
    global FASTA_LABEL_SYMBOL
    label, buffer = None, []
    for line in f:
        if len(line) > 0 and line.startswith(FASTA_LABEL_SYMBOL):
            if len(buffer) > 0:
                yield label, "".join(buffer)
            buffer = []
            label = line.strip()[1:]
        else:
            buffer.append(line.strip())
    if len(buffer) > 0:
        yield label, "".join(buffer)


def is_genomic_alphabet(s: str):
    return re.fullmatch(rf"[{SIMPLE_ALPHABET}]*", s) is not None


def load_virus_genome(filename: str):
    with open(filename) as f:
        full_genome = list(load_fasta(f))
    # virus genome is a single chromosome
    assert len(full_genome) == 1, \
        f"Virus genome length is {len(full_genome)}, not 1!"
    genome = full_genome[0][1].upper()

    assert is_genomic_alphabet(genome), \
        f"The loaded genome {genome[:12]}... contains letters outside " \
        f"of {SIMPLE_ALPHABET} alphabet!"
    return genome


def _load_variants_raw(f: TextIO, alphabet: str):
    reader = csv.reader(f, delimiter='\t')
    raw_data = defaultdict(lambda: defaultdict(lambda: [0.0 for _ in alphabet]))
    for row_num, row in enumerate(reader):
        variant, position, letter_num, count = \
            row[0], int(row[1]), alphabet.index(row[2].upper()), float(row[3])
        assert letter_num != -1, \
            f"Letter {row[2]} at line {row_num + 1} is not in the alphabet {alphabet}!"
        assert position >= 0, \
            f"Position {position} at line {row_num + 1} is below zero!"
        raw_data[variant][position][letter_num] = count
    if len(raw_data) == 0:
        raise Exception("The variants' file is empty!")
    variant_names = list(raw_data.keys())
    length = 1 + max(max(variant_data.keys()) for variant_data in raw_data.values())
    table = np.array([[d[pos] for pos in np.arange(length)] for d in raw_data.values()])
    return variant_names, table


def _normalized_variants(table, percentile=0.0):
    variant_num, length, alphabet_size = table.shape
    result = np.zeros(table.shape)

    for v in np.arange(variant_num):
        quantile_coverage = np.quantile(np.sum(table[v], 1), percentile)

        for p in np.arange(length):
            total = np.sum(table[v][p])
            if total > 0:
                result[v][p] = table[v][p] / max(total, quantile_coverage)
            else:
                result[v][p] = np.array([1 / alphabet_size for _ in range(alphabet_size)])
    return result


def _remove_letters_from_variants(table, positions):
    result = np.delete(table, positions, 2)
    return result


def _positions_of_missing_letters(a, b):
    result = []
    for pos, letter in enumerate(a):
        if letter not in b:
            result.append(pos)
    return result


def load_variants_with_simple_alphabet(filename: str):
    """load normalized variants over an alphabet {A,C,G,T}.

    :returns: variant names: list[str], variant frequencies: numpy array (variant, position, letter)
    """
    # this is the function we would use for some time
    with open(filename, newline='') as f:
        variant_names, table = _load_variants_raw(f, FULL_ALPHABET)
    positions_to_remove = _positions_of_missing_letters(FULL_ALPHABET, SIMPLE_ALPHABET)
    filtered_table = _remove_letters_from_variants(table, positions_to_remove)
    normalized_table = _normalized_variants(filtered_table)
    return variant_names, normalized_table


def load_variants_with_simple_alphabet_punish_gaps(filename: str, gap_percentile: float = 0.25):
    """load normalized variants over an alphabet {A,C,G,T}.
    positions with coverage less than `gap_percentile` quantile would be normalised to quantile coverage
    (so that their sum would be coverage / quantile)

    :returns: variant names: list[str], variant frequencies: numpy array (variant, position, letter)
    """
    # this is the function we would use for some time
    with open(filename, newline='') as f:
        variant_names, table = _load_variants_raw(f, FULL_ALPHABET)
    positions_to_remove = _positions_of_missing_letters(FULL_ALPHABET, SIMPLE_ALPHABET)
    filtered_table = _remove_letters_from_variants(table, positions_to_remove)
    normalized_table = _normalized_variants(filtered_table, gap_percentile)
    return variant_names, normalized_table


def _load_alignments_raw(filename: str):
    # `until_eof=True` turns off the usage of an index file `bam.bai`
    return pysam.AlignmentFile(filename).fetch(until_eof=True)


def _filter_good_alignments(alignments):
    def is_good(al: pysam.AlignedSegment):
        return not any((al.is_secondary,
                        al.is_supplementary,
                        al.is_unmapped,
                        al.mate_is_unmapped,))

    filtered_count, total_count = 0, 0
    for al in alignments:
        total_count += 1
        if is_good(al):
            yield al
        else:
            filtered_count += 1
    if filtered_count > 0:
        logging.info(f"{filtered_count} alignments filtered (total {total_count} alignments)")


def load_alignments(filename: str):
    """load alignments while filtering secondary, supplementary and (partially) unaligned"""
    all_alignments = _load_alignments_raw(filename)
    good_alignments = _filter_good_alignments(all_alignments)
    return good_alignments


def get_genome_size_from_alignments(filename: str):
    lengths = pysam.AlignmentFile(filename).lengths
    assert len(lengths) == 1, "Alignment file contains more than one reference sequence!"
    return lengths[0]


def apply_to_cigartuples(fun, alignment, *args, **kwargs):
    """
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9 (????!)
    """
    query_pos = 0
    reference_pos = alignment.reference_start
    for op, length in alignment.cigartuples:
        fun(op, length, reference_pos, query_pos, alignment, *args, **kwargs)
        if op == 0 or op == 7 or op == 8:
            reference_pos += length
            query_pos += length
        elif op == 1 or op == 4:
            query_pos += length
        elif op == 2 or op == 3:
            reference_pos += length
        elif op == 5 or op == 6:
            pass
        else:
            raise Exception(f"Operation code of cigar tuple is outside of range [0-8]: "
                            f"op={op}, length={length}")


def load_posteriors(f: TextIO):
    """Load variant names and posteriors without identifiers"""
    reader = csv.reader(f, delimiter="\t")
    header = reader.__next__()
    variant_names = header[1:]
    posteriors = []
    for row in reader:
        posterior = list(map(float, row[1:]))
        posteriors.append(posterior)
    posteriors_np = np.array(posteriors, dtype=np.float64)
    return variant_names, posteriors_np


def load_posteriors_with_identifiers(f: TextIO):
    """Load variant names and posteriors with identifiers"""
    reader = csv.reader(f, delimiter="\t")
    header = reader.__next__()
    variant_names = header[1:]
    posteriors = []
    identifiers = []
    for row in reader:
        posterior = list(map(float, row[1:]))
        identifiers.append(row[0])
        posteriors.append(posterior)
    posteriors_np = np.array(posteriors,
                             dtype=np.float64)
    return variant_names, identifiers, posteriors_np


def count_simple_coverage(genome_length: int, alignments):
    diffs = [0 for _ in range(genome_length + 1)]
    for al in alignments:
        assert 0 <= al.reference_start < genome_length, \
            f"alignment start {al.reference_start} is outside of range [0, {genome_length})"
        diffs[al.reference_start] += 1
        assert 0 <= al.reference_end <= genome_length, \
            f"alignment end {al.reference_end} is outside of range [0, {genome_length})"
        assert al.reference_start <= al.reference_end, \
            f"alignment start {al.reference_start} is after the alignment end {al.reference_end}"
        diffs[al.reference_end] -= 1

    result = list(itertools.accumulate(diffs))[:-1]
    return result


def load_profile(filename):
    raw_lines = []
    max_pos = 0
    with open(filename, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            pos, base_num, value = int(row[0]), L2N[row[1]], int(row[2])
            raw_lines.append((pos, base_num, value))
            max_pos = max(max_pos, pos)
    result = np.zeros((max_pos+1, 4), dtype=int)
    for pos, base_num, value in raw_lines:
        result[pos, base_num] = value
    return result


def apply_subst_noise(p, subst_rate):
    return (1 - subst_rate) * p + (1 - p) * subst_rate/3


def load_masked_positions(f):
    POS_COLNUM = 1
    FILTER_COLNUM = 6
    MASKING = "mask"
    masked_positions = []
    for line in f:
        if line.startswith("##"):
            # commented line
            continue
        row = line.strip().split("\t")
        if row[FILTER_COLNUM] == MASKING:
            pos = int(row[POS_COLNUM])-1  # convert from 1-based to 0-based
            masked_positions.append(pos)
    return masked_positions