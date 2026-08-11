"""
crossfilt-vcf: Generate a VCF from consensus-coordinate genomes.

Takes a consensus genome and N species genomes (all in consensus coordinates) and
writes a biallelic VCF with diploid GT calls for every SNP position in the consensus.
"""
import sys
import argparse
import pysam
import importlib.metadata

__version__ = importlib.metadata.version('crossfilt')


def main():
    parser = argparse.ArgumentParser(
        prog='consensus-vcf',
        description='Takes a consensus genome and consensus-converted alternate genomes and makes a VCF for identifying species.')

    parser.add_argument("-c", "--consensus",     required=True, help="The N-way consensus genome")
    parser.add_argument("-g", "--genomes",       required=True, help="Genomes in consensus coordinates, separated by commas")
    parser.add_argument("-l", "--labels",        required=True, help="Column headers to use in the VCF file for these genomes, separated by commas")
    parser.add_argument("-o", "--output-prefix", required=True, help="The prefix for the output files")
    parser.add_argument('--version', action='version', version='CrossFilt v{version}'.format(version=__version__))

    args = parser.parse_args()

    genome_files   = args.genomes.split(",")
    labels         = args.labels.split(",")
    output         = args.output_prefix
    consensus_file = args.consensus

    if len(labels) != len(genome_files):
        sys.exit("Error: --labels has " + str(len(labels)) + " entries but --genomes has " +
                 str(len(genome_files)) + ". They must correspond one-to-one.")

    consensus_genome = pysam.Fastafile(consensus_file)
    genomes = [pysam.Fastafile(g) for g in genome_files]

    vcfh = pysam.VariantHeader()
    for label in labels:
        vcfh.add_sample(label)
    for chrom in consensus_genome.references:
        vcfh.add_meta('contig', items=[('ID', chrom)])
    vcfh.add_meta('FORMAT', items=[('ID', "GT"), ('Number', 1), ('Type', 'String'),
                                   ('Description', 'Genotype')])

    with pysam.VariantFile(output + ".vcf", "w", header=vcfh) as vcf:
        for chrom in consensus_genome.references:
            consensus_chr = consensus_genome.fetch(chrom).upper()
            other_chrs = [genome.fetch(chrom).upper() for genome in genomes]
            l = len(consensus_chr)

            for i in range(1, l - 1):
                if consensus_chr[i] != 'N':
                    continue
                if consensus_chr[i - 1] == 'N':
                    continue
                if consensus_chr[i + 1] == 'N':
                    continue

                bases  = [this_chr[i] for this_chr in other_chrs]
                ubases = sorted(set(bases))
                if len(ubases) != 2:
                    continue

                r = vcf.new_record(contig=chrom, start=i, stop=i + 1,
                                   alleles=(ubases[0], ubases[1]))

                for j, label in enumerate(labels):
                    r.samples[label]['GT'] = (0, 0) if bases[j] == ubases[0] else (1, 1)

                vcf.write(r)


if __name__ == '__main__':
    main()
