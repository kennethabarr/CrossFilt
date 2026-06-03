"""
crossfilt-merge: Merge N consensus genomes (all in the same coordinate system) into one.

Takes two or more consensus genomes produced by consensus-pair and outputs a new
consensus where any position that is N in any input genome is set to N in the output.
"""
import argparse
import pysam
from timeit import default_timer as timer
import importlib.metadata

__version__ = importlib.metadata.version('crossfilt')


def main():
    parser = argparse.ArgumentParser(
        prog='consensus-merge',
        description='Takes consensus genomes (in the same coordinate system) and builds a new consensus.')

    parser.add_argument("-c", "--consensus",     required=True, help="A consensus genome. Use argument multiple times for an N-way consensus.", action='append')
    parser.add_argument("-o", "--output-prefix", required=True, help="The prefix for the output files")
    parser.add_argument('--version', action='version', version='CrossFilt v{version}'.format(version=__version__))

    args = parser.parse_args()

    genome_files = args.consensus
    output       = args.output_prefix

    if len(genome_files) < 2:
        raise Exception("Error: You must specify at least two genomes with option --consensus (-c)")

    consensus_genome = pysam.Fastafile(genome_files[0])
    other_files = genome_files[1:]

    genomes = [pysam.Fastafile(g) for g in other_files]

    with open(output + "_consensus.fa", "w") as ofile:
        for chrom in consensus_genome.references:
            start = timer()
            print("Processing chromosome " + chrom, file=sys.stderr)

            consensus_chr = list(consensus_genome.fetch(chrom).upper())
            write = True

            for genome in genomes:
                if not write:
                    continue

                if chrom not in genome.references:
                    print("Warning: chromosome " + chrom + " not found in all genomes. " +
                          "Were these all produced using consensus-pair using the" +
                          " same reference?", file=sys.stderr)
                    write = False
                    continue

                this_chr   = genome.fetch(chrom).upper()
                chr_length = len(consensus_chr)

                if len(this_chr) != chr_length:
                    raise Exception("Error: chromosome " + chrom + " is not the same length. " +
                                    "Were these all produced using consensus-pair using the" +
                                    " same reference?")

                for i in range(chr_length):
                    if consensus_chr[i] == 'N':
                        continue
                    if this_chr[i] == 'N':
                        consensus_chr[i] = 'N'

            consensus_chr = ''.join(consensus_chr)
            end = timer()
            print("Completed in", round(end - start, 2), "seconds\n", file=sys.stderr)
            if write:
                ofile.write(">" + chrom + "\n" + consensus_chr + "\n")


if __name__ == '__main__':
    main()
