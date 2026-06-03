"""
crossfilt-consensus: Build a 2-way consensus genome in the coordinate system of the reference.

Compares reference and alternate genomes base-by-base across chain-mapped intervals.
Positions where both genomes agree are filled in; positions with mismatches are left as N
in the consensus and filled separately in the reference and alternate outputs.
"""
import sys
import argparse
import pysam
from timeit import default_timer as timer
from .lib import liftover_functions as lift
import importlib.metadata

__version__ = importlib.metadata.version('crossfilt')


def main():
    parser = argparse.ArgumentParser(
        prog='consensus-pair',
        description='Builds a 2-way consensus genome in the coordinate system of the reference genome')

    parser.add_argument("-r", "--reference",     required=True, help="The fasta sequence of the reference genome")
    parser.add_argument("-a", "--alternate",     required=True, help="The fasta sequence of the alternate genome")
    parser.add_argument("-c", "--chain",         required=True, help="The reference to alternate chain file")
    parser.add_argument('-l', "--length",        required=False, type=int, default=6, help='The number of bases around indels to mask')
    parser.add_argument("-o", "--output-prefix", required=True, help="The prefix for the output files")
    parser.add_argument('--version', action='version', version='CrossFilt v{version}'.format(version=__version__))

    args = parser.parse_args()

    reference   = args.reference
    alternate   = args.alternate
    chainfile   = args.chain
    mask_length = args.length
    output      = args.output_prefix

    start = timer()
    print("Reading chain file", file=sys.stderr)

    REFFASTA = pysam.Fastafile(reference)
    ALTFASTA = pysam.Fastafile(alternate)

    maps, tSizes, qSizes = lift.read_chain_file(chainfile, REFFASTA.references, ALTFASTA.references)

    end = timer()
    print("Completed in", round(end - start, 2), "seconds\n", file=sys.stderr)

    def update_interval_consensus(interval, ref_chr, consensus_chr, reference_chr, alternate_chr):
        ref_start  = interval.begin
        ref_end    = interval.end

        alt_chr    = interval.data[0]
        alt_start  = interval.data[1]
        alt_end    = interval.data[2]
        alt_strand = interval.data[3]

        ref_seq = REFFASTA.fetch(ref_chr, ref_start, ref_end).upper()
        alt_seq = ALTFASTA.fetch(alt_chr, alt_start, alt_end).upper()

        if alt_strand == "-":
            alt_seq = lift.revcomp_DNA(alt_seq)

        l = len(ref_seq)

        ndiff = 0
        nsame = 0

        for i in range(mask_length, l - mask_length):
            this_base = ref_seq[i]
            alt_base  = alt_seq[i]
            if this_base == alt_base:
                nsame += 1
                consensus_chr[i + ref_start] = this_base
                reference_chr[i + ref_start] = this_base
                alternate_chr[i + ref_start] = this_base
            else:
                ndiff += 1
                reference_chr[i + ref_start] = this_base
                alternate_chr[i + ref_start] = alt_base

        tot = ndiff + nsame

        if tot > 100:
            if nsame / tot < 0.5:
                print("Warning: The interval " +
                      ref_chr + ":" + str(ref_start) + "-" + str(ref_end) +
                      " had " + str(round(100 * nsame / tot)) +
                      " identity! There may be a problem with your sequence or chain files", file=sys.stderr)
                print(ref_seq)
                print(alt_seq)
        return 1

    con_out = open(output + "_consensus.fa", "w")
    ref_out = open(output + "_reference.fa", "w")
    alt_out = open(output + "_alternate.fa", "w")

    for this_chr, this_len in tSizes.items():
        if "alt" in this_chr:
            continue

        start = timer()
        print("Processing chromosome " + this_chr, file=sys.stderr)

        consensus_chr = list("N" * this_len)
        reference_chr = list("N" * this_len)
        alternate_chr = list("N" * this_len)

        this_chains = lift.get_chr_chains(maps, this_chr)
        chains = sorted(this_chains, key=lambda chain: chain.data['score'])

        for this_chain in chains:
            n = this_chain.end - this_chain.begin
            consensus_chr[this_chain.begin:this_chain.end] = list("N" * n)
            reference_chr[this_chain.begin:this_chain.end] = list("N" * n)
            alternate_chr[this_chain.begin:this_chain.end] = list("N" * n)

            for interval in this_chain.data['mapTree'].items():
                update_interval_consensus(interval,
                                          this_chr,
                                          consensus_chr,
                                          reference_chr,
                                          alternate_chr)

        end = timer()
        print("Completed in", round(end - start, 2), "seconds\n", file=sys.stderr)

        start = timer()
        print("Validating this chromosome", file=sys.stderr)
        ref_seq = REFFASTA.fetch(this_chr, 0, this_len).upper()

        if len(ref_seq) != this_len:
            print("Error: reference sequence and consensus sequence were not the same length for " + this_chr, file=sys.stderr)
            con_out.close(); ref_out.close(); alt_out.close()
            sys.exit(1)

        ref_ns       = ref_seq.count('N')
        consensus_seq = ''.join(consensus_chr)
        consensus_ns  = consensus_seq.count('N')

        for i in range(this_len):
            if consensus_seq[i] != 'N' and consensus_seq[i] != ref_seq[i]:
                print("Error: reference sequence and consensus sequence not identical for " +
                      this_chr + " " + str(i), file=sys.stderr)

        print("Reference " + this_chr + " is " + str(round(100 * ref_ns / this_len)) + "% N", file=sys.stderr)
        print("Consensus " + this_chr + " is " + str(round(100 * consensus_ns / this_len)) + "% N", file=sys.stderr)

        end = timer()
        print("Completed in", round(end - start, 2), "seconds\n", file=sys.stderr)

        con_out.write(">" + this_chr + "\n" + consensus_seq + "\n")
        ref_out.write(">" + this_chr + "\n" + ''.join(reference_chr) + "\n")
        alt_out.write(">" + this_chr + "\n" + ''.join(alternate_chr) + "\n")

    con_out.close()
    ref_out.close()
    alt_out.close()


if __name__ == '__main__':
    main()
