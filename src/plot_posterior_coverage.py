import argh
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from helpers import load_alignments, load_posteriors_with_identifiers, \
    get_genome_size_from_alignments, count_simple_coverage


def filter_alignments_by_posterior(alignments, identifiers, posteriors, threshold):
    good_ids = {identifiers[i] for i in range(len(identifiers)) if posteriors[i] >= threshold}
    for al in alignments:
        if al.query_name in good_ids:
            yield al


def main(bam_filename, posteriors_filename, output_filename, ylog=False):
    genome_length = get_genome_size_from_alignments(bam_filename)
    alignments = list(load_alignments(bam_filename))
    with open(posteriors_filename, newline="") as f:
        variant_names, identifiers, posteriors = load_posteriors_with_identifiers(f)

    fig, ax = plot_posterior_coverage(alignments, genome_length, identifiers, posteriors, variant_names, ylog=ylog)

    fig.savefig(output_filename, bbox_inches="tight")


def plot_posterior_coverage(alignments, genome_length, identifiers, posteriors, variant_names, ylog=False):
    variant_count = len(variant_names)
    levels = [0, 1/len(variant_names), 0.5, 0.75, 0.95]
    labels = ["all", ">prior", ">50%", ">75%", ">95%"]
    colors = ["grey", "purple", "blue", "orange", "red"]
    linestyles = ["solid", "dotted", "solid", "solid", "solid"]
    fig_height = variant_count * 0.7 + 1
    fig_width = 10
    ylimbottom = 0 if not ylog else 10

    fig, ax = plt.subplots(figsize=(fig_width, fig_height),
                           nrows=variant_count,
                           sharex="all",
                           sharey="all")
    for v in np.arange(variant_count):
        for lnum in range(len(levels)):
            level = levels[lnum]
            filtered_als = filter_alignments_by_posterior(alignments, identifiers, posteriors[:, v], level)
            coverage = count_simple_coverage(genome_length, filtered_als)
            ax[v].plot(np.arange(1, genome_length + 1), coverage, color=colors[lnum], linestyle=linestyles[lnum])
        ax[v].set_ylim(bottom=ylimbottom)
        ax[v].set_ylabel(variant_names[v], rotation=0, labelpad=20)
        if ylog:
            ax[v].set_yscale('log')
    ax[-1].set_xlabel("Reference (bp)")

    legend_elements = [Patch(facecolor=colors[i], label=labels[i], edgecolor=colors[i])
                       for i in range(len(colors))]
    ax[0].legend(handles=legend_elements,
                 ncol=len(colors),
                 loc="lower center",
                 bbox_to_anchor=(0., 1.02, 1., .102),
                 borderaxespad=0.
                 )

    return fig, ax


if __name__ == "__main__":
    argh.dispatch_command(main)
