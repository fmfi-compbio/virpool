import argh
import pysam

import helpers


def main(*filenames, output="output.bam"):
    with pysam.AlignmentFile(output, "wb", template=pysam.AlignmentFile(filenames[0])) as output_bam:
        for filename in filenames:
            for alignment in helpers.load_alignments(filename):
                output_bam.write(alignment)


if __name__ == "__main__":
    argh.dispatch_command(main)