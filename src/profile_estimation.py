import csv
import itertools
import json
from collections import defaultdict

import argh
import numba
import numpy as np
import pysam
import logging, logging.handlers

import helpers


@numba.njit
def expand_variant_name(variant_name: str, aliases, separator: str = "."):
    if len(variant_name) == 0:
        return ""

    separator_position = variant_name.find(separator)
    if separator_position == -1:
        # the separator is not present in the variant_name
        separator_position = len(variant_name)

    root = variant_name[:separator_position]
    root_expansion = aliases.get(root, root)
    if root_expansion == "":
        # forbidden root
        return ""

    expanded_name = root_expansion + variant_name[separator_position:]
    return expanded_name


@numba.njit
def determine_variant_group(sample_variant, variant_groups_list, aliases):
    expanded_sample_variant = expand_variant_name(sample_variant, aliases)
    if expanded_sample_variant == "":
        return -1

    result = len(variant_groups_list)
    best_match_length = 0
    for num, name in enumerate(variant_groups_list):
        #logger.debug(f"Checking variant #{num} {name}...")
        if expanded_sample_variant == name or expanded_sample_variant.startswith(name + "."):
            #logger.debug(f"group {name} matches {expanded_sample_name}")
            if best_match_length < len(name):
                #logger.debug(f"Length of {name} is superior!")
                result = num
                best_match_length = len(name)
    #logger.debug(f"{result=}")
    return result


@argh.arg("-o", "--output")
@argh.arg("--date-min", help="Format YYYY-MM-DD")
@argh.arg("--date-max", help="Format YYYY-MM-DD")
def main(variant_groups,
         aliases,
         metadata,
         alignments,
         output="variants.tsv",
         reference=None,
         sufficient_matching_length: int = 0,
         process_supplementary_alignments: bool = True,
         log=None,
         date_min: str = None,
         date_max: str = None,
         without_other: bool = False
         ):
    set_root_logger(log)
    logger = logging.getLogger("profile_estimation")

    with open(variant_groups) as f:
        variant_groups_list = load_variant_groups_list(f)
    logger.debug(f"{variant_groups_list=}")

    with open(aliases) as f:
        aliases_data = json.load(f)
        del aliases_data["A"]
        del aliases_data["B"]
        for name in aliases_data:
            if isinstance(aliases_data[name], list):
                aliases_data[name] = ""
        aliases_numbed = numba.typed.Dict()
        for name, expansion in aliases_data.items():
            aliases_numbed[name] = expansion
    logger.debug(f"{aliases_data=}")

    variant_groups_list_expanded = list(map(
        lambda vname: expand_variant_name(vname, aliases_numbed),
        variant_groups_list))
    logger.debug(f"{variant_groups_list_expanded=}")
    variant_groups_list_numbed = numba.typed.List(variant_groups_list_expanded)

    with open(metadata, newline="") as f:
        sample_variants = load_sample_variants_with_date(f)

    genome_length = helpers.get_genome_size_from_alignments(alignments)
    logging.debug(f"Genome length: {genome_length}")

    profiles = np.zeros((len(variant_groups_list)+1, genome_length * 4), dtype=int)

    with pysam.AlignmentFile(alignments, reference_filename=reference) as bam:
        alignments_processed_count = 0
        prev_id = None
        sample_names_processed_count = 0
        accepted_alignments_count = 0

        skipped_time = 0
        skipped_no_alignment = 0
        skipped_supplementary = 0
        skipped_length = 0
        skipped_unknown_variant_group = 0
        skipped_unknown_query_seq = 0
        unknown_variants = defaultdict(int)
        variant_groups_counts = defaultdict(int)

        if without_other:
            variant_group_names = variant_groups_list
        else:
            variant_group_names = variant_groups_list + ['other']

        def debug_message():
            nonlocal alignments_processed_count, sample_names_processed_count, accepted_alignments_count
            nonlocal skipped_no_alignment, skipped_supplementary, skipped_length, skipped_time
            nonlocal skipped_unknown_query_seq, skipped_unknown_variant_group
            message = f"Processed: {alignments_processed_count} " \
                      f"unique ids: {sample_names_processed_count} " \
                      f"accepted: {accepted_alignments_count} " \
                      f"skipped no_aln: {skipped_no_alignment} " \
                      f"supplement: {skipped_supplementary} " \
                      f"length: {skipped_length} " \
                      f"time: {skipped_time} " \
                      f"unknown group: {skipped_unknown_variant_group} " \
                      f"unknown query seq: {skipped_unknown_query_seq}"
            logger.debug(message)
            logger.debug(f"Variant groups {dict(variant_groups_counts)}")
            # if len(unknown_variants) > 0:
            #     logger.debug(f"Unknown variants: {dict(unknown_variants)}")

        for alignment in bam.fetch(until_eof=True):
            if alignments_processed_count % 10000 == 0:
                debug_message()
            alignments_processed_count += 1

            if alignment.query_name != prev_id:
                sample_names_processed_count += 1
                prev_id = alignment.query_name

            if alignment.cigartuples is None:
                skipped_no_alignment += 1
                continue

            if not process_supplementary_alignments and alignments.is_supplementary:
                skipped_supplementary += 1
                continue

            if matching_length(alignment) < sufficient_matching_length:
                skipped_length += 1
                continue

            sample_variant, sample_date = sample_variants[alignment.query_name]
            if date_min is not None and sample_date < date_min:
                skipped_time += 1
                continue
            if date_max is not None and sample_date > date_max:
                skipped_time += 1
                continue

            #logger.debug(f"{sample_variant=}")

            variant_group = determine_variant_group(sample_variant, variant_groups_list_numbed, aliases_numbed)

            if variant_group not in range(len(variant_group_names)):
                skipped_unknown_variant_group += 1
                unknown_variants[sample_variant] += 1
                continue

            query_sequence = alignment.query_sequence

            if query_sequence is None:
                skipped_unknown_query_seq += 1
                continue

            accepted_alignments_count += 1
            variant_groups_counts[variant_group_names[variant_group]] += 1
            query_pos = 0

            reference_pos = alignment.reference_start
            for op, length in alignment.cigartuples:
                if op == 0 or op == 7 or op == 8:
                    add_counts(query_pos, length, reference_pos, query_sequence, profiles[variant_group])
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
        else:
            debug_message()

    with open(output, "w") as f:
        dump_profiles(f, genome_length, profiles, variant_group_names)


def dump_profiles(f, genome_length, profiles, variant_group_names):
    for vnum, vname in enumerate(variant_group_names):
        for pos in range(genome_length):
            for lnum, letter in enumerate("ACGT"):
                print(f"{vname}\t{pos}\t{letter}\t{profiles[vnum][pos * 4 + lnum]}", file=f)


def load_variant_groups_list(f):
    variant_groups_list = []
    for line in f:
        name = line.strip()
        if name:
            variant_groups_list.append(name)
    return variant_groups_list


def load_sample_variants(f):
    reader = csv.reader(f, delimiter="\t")
    header = reader.__next__()
    sample_name_colnum = header.index("Virus_name")
    variant_colnum = header.index("Pango lineage")

    result = dict()

    for row in reader:
        sample_name = row[sample_name_colnum]
        variant = row[variant_colnum]
        result[sample_name] = variant

    return result


def load_sample_variants_with_date(f):
    logger = logging.getLogger("profile_estimation")
    reader = csv.reader(f, delimiter="\t")
    header = reader.__next__()

    try:
        sample_name_colnum = header.index("Virus_name")
    except ValueError:
        sample_name_colnum = header.index("Virus name")

    variant_colnum = header.index("Pango lineage")
    date_colnum = header.index("Collection date")

    result = dict()

    for rnum, row in enumerate(reader):
        if rnum % 100000 == 0:
            logger.debug(f"Loaded {rnum} rows from metadata...")
        sample_name = row[sample_name_colnum]
        variant = row[variant_colnum]
        date = row[date_colnum]
        result[sample_name] = (variant, date)

    return result


def matching_length(alignment):
    match_length = 0
    for op, length in alignment.cigartuples:
        if op == 0 or op == 7 or op == 8:  # matching
            match_length += length
    return match_length


@numba.njit
def add_counts(query_pos, length, reference_pos, query_sequence, profile):
    for i in range(length):
        letter = query_sequence[query_pos + i]
        lnum = -1
        if letter == "A":
            lnum = 0
        elif letter == "C":
            lnum = 1
        elif letter == "G":
            lnum = 2
        elif letter == "T":
            lnum = 3
        else:
            continue

        profile_coord = 4 * (reference_pos+i) + lnum
        profile[profile_coord] += 1


def set_root_logger(log_filename):
    logger = logging.getLogger("profile_estimation")
    add_console_handler(logger)
    logger.setLevel(logging.DEBUG)

    # if log_filename is not None:
    #     add_file_handler(log_filename, logger, level=logging.DEBUG)


def add_file_handler(log_filename, logger, level=logging.INFO):
    handler = logging.handlers.WatchedFileHandler(log_filename)
    handler.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def add_console_handler(logger, level=logging.DEBUG):
    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


if __name__ == "__main__":
    argh.dispatch_command(main)
