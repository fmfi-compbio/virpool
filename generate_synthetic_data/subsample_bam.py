import logging
import random

import argh
import numpy as np
from pysam import AlignmentFile

from helpers import load_alignments, \
    get_genome_size_from_alignments, count_simple_coverage


@argh.arg("coverage", type=float)
def main(bam_filename, output_filename, coverage: float,
         sampled_reads_count_filename=None):
    sampling_ratio = eval_sampling_ratio(bam_filename, coverage)
    alignments = load_alignments(bam_filename)

    is_read_chosen = dict()
    with AlignmentFile(output_filename, "wb", template=AlignmentFile(bam_filename)) as f:
        first_alignment = None
        written_alignment_count = 0
        for al in alignments:
            if first_alignment is None:
                first_alignment = al
            if al.query_name not in is_read_chosen:
                is_read_chosen[al.query_name] = (random.random() < sampling_ratio)
            if is_read_chosen[al.query_name]:
                f.write(al)
                written_alignment_count += 1
        if written_alignment_count == 0:
            f.write(first_alignment)

    if sampled_reads_count_filename is not None:
        with open(sampled_reads_count_filename, "w") as f:
            print(written_alignment_count, file=f)


def eval_sampling_ratio(bam_filename, coverage):
    average_coverage = get_average_coverage(bam_filename)
    sampling_ratio = min(1, coverage / average_coverage)
    return sampling_ratio


def get_average_coverage(bam_filename):
    genome_length = get_genome_size_from_alignments(bam_filename)
    alignments = load_alignments(bam_filename)
    base_coverage = count_simple_coverage(genome_length, alignments)
    average_coverage = np.average(base_coverage)
    return average_coverage


if __name__ == "__main__":
    argh.dispatch_command(main)
