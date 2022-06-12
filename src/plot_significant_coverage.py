from pprint import pprint

import argh
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from helpers import load_alignments, load_variants_with_simple_alphabet, \
    load_posteriors_with_identifiers, load_virus_genome, L2N, count_simple_coverage
from plot_posterior_coverage import filter_alignments_by_posterior


def eval_permutation(old_order: list, new_order: list):
    result = [-1 for _ in old_order]
    for i, item in enumerate(old_order):
        result[i] = new_order.index(item)
    return result


def rearrange_posterior_columns(table, old_order, new_order):
    permutation = eval_permutation(old_order, new_order)
    result = table[:, permutation]
    return result


def find_significant_positions(genome, variants_table, threshold=0.25):
    """Return positions where the referential frequency is below the threshold"""
    variant_count = variants_table.shape[0]
    results = [[] for _ in range(variant_count)]
    for vnum in range(variant_count):
        for i, letter in enumerate(genome):
            lnum = L2N[letter]
            if variants_table[vnum][i][lnum] < threshold:
                results[vnum].append(i)
    return results


def main(reads_filename,
         posteriors_filename,
         variants_filename,
         genome_filename,
         output_filename):
    alignments = list(load_alignments(reads_filename))
    with open(posteriors_filename) as f:
        variant_names_posteriors, identifiers, posteriors = load_posteriors_with_identifiers(f)
    variant_names, variants_table = load_variants_with_simple_alphabet(variants_filename)
    posteriors = rearrange_posterior_columns(posteriors, variant_names_posteriors, variant_names)
    genome = load_virus_genome(genome_filename)
    significant_positions = find_significant_positions(genome, variants_table)

    fig, ax = plot_significant_coverage(alignments, genome, identifiers, posteriors, significant_positions, variant_names)
    fig.savefig(output_filename, bbox_inches="tight")


def plot_significant_coverage(alignments, genome, identifiers, posteriors, significant_positions, variant_names):
    genome_length = len(genome)
    levels = [0, 0.5, 0.75, 0.95]
    labels = ["all", ">50%", ">75%", ">95%"]
    colors = ["grey", "blue", "orange", "red"]
    variant_count = len(variant_names)
    fig_height = (variant_count // 2 + variant_count % 2) * 1.5 + 1
    fig_width = 15
    matplotlib.rcParams.update({'font.size': 10})
    fig, ax = plt.subplots(figsize=(fig_width, fig_height),
                           ncols=2, nrows=variant_count // 2 + variant_count % 2,
                           sharey="all")

    # if there are only two variants, `nrows` variable is 1, and therefore
    # ax is not a table with two dimensions, but a list.
    # in order to preserve the original code, we add another layer of array.
    # yep, it's not **clean**, but so what. it's a bioinformatic tool after all.
    if variant_count <= 2:
        ax = [ax]

    for v in reversed(np.arange(variant_count)):
        ps = significant_positions[v]
        if len(ps) == 0:
            continue

        for lnum in range(len(levels)):
            level = levels[lnum]
            filtered_als = filter_alignments_by_posterior(alignments, identifiers, posteriors[:, v], level)
            coverage = count_simple_coverage(genome_length, filtered_als)
            cs = np.array([coverage[p] for p in ps])
            ax[v // 2][v % 2].bar(list(map(lambda x: str(x + 1), ps)), cs, color=colors[lnum])
            ax[v // 2][v % 2].tick_params(axis="x", rotation=-45)
        ax[v // 2][v % 2].set_ylim(bottom=0)
        ax[v // 2][v % 2].set_ylabel(variant_names[v], rotation=0, labelpad=20)
    if variant_count % 2 == 1:
        fig.delaxes(ax[-1][1])

    legend_elements = [Patch(facecolor=colors[i], label=labels[i], edgecolor=colors[i])
                       for i in range(len(colors))]
    fig.legend(handles=legend_elements,
               ncol=len(colors),
               loc="upper center",
               borderaxespad=0.)

    fig.tight_layout(pad=1)
    return fig, ax


if __name__ == "__main__":
    argh.dispatch_command(main)
