"""
consensus-demux: Build a joint species+individual demultiplexing VCF.

Takes the species-level VCF from consensus-vcf, one individual VCF per species,
and a chain file for each non-reference species. Writes a single VCF with one
sample column per donor (species x individual) for use with demuxlet or vireo.

Species-level positions (N in the consensus) get fixed homozygous alleles per
species. Individual-level positions (non-N) get actual genotypes for the
relevant species and REF homozygous for all others.
"""
import sys
import argparse
import pysam
import importlib.metadata

from crossfilt.lib.liftover_functions import read_chain_file, get_chr_chains

__version__ = importlib.metadata.version('crossfilt')

_complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def _comp(base):
    return base.translate(_complement)


def liftover_position(chrom, start0, maps):
    """
    Lift a single 0-based SNP position to consensus coordinates.
    Returns (cons_chrom, cons_start0, is_reverse) or None if unmappable.
    """
    chr_chains = get_chr_chains(maps, chrom)
    if chr_chains is None:
        return None

    end0 = start0 + 1
    top_chains = sorted(
        [c for c in chr_chains.find(start0, end0)
         if c.value['tStart'] <= start0 and c.value['tEnd'] >= end0],
        key=lambda c: -c.value['score'])
    if not top_chains:
        return None

    sub_ivs = [iv for iv in top_chains[0].value['mapTree'].find(start0, end0)
               if iv.start <= start0 and iv.end >= end0]
    if not sub_ivs:
        return None

    iv = sub_ivs[0]
    offset     = start0 - iv.start
    is_reverse = iv.value[3] == '-'
    q_start0   = iv.value[2] - 1 - offset if is_reverse else iv.value[1] + offset
    return (iv.value[0], q_start0, is_reverse)


def _vcf_contigs(vcf_path):
    """Return contig names from VCF header; fall back to scanning records."""
    with pysam.VariantFile(vcf_path) as vcf:
        contigs = list(vcf.header.contigs)
        if contigs:
            return contigs
        seen = set()
        for rec in vcf:
            seen.add(rec.chrom)
    return list(seen)


def load_species_vcf(vcf_path):
    """
    Parse the consensus-vcf output.
    Returns (samples, variants) where variants maps (chrom, start0) to a dict
    containing 'alleles' and a per-species allele index (int).
    """
    variants = {}
    with pysam.VariantFile(vcf_path) as vcf:
        samples = list(vcf.header.samples)
        for rec in vcf:
            entry = {'alleles': rec.alleles}
            for sp in samples:
                entry[sp] = rec.samples[sp]['GT'][0]
            variants[(rec.chrom, rec.start)] = entry
    return samples, variants


def load_individual_vcf(vcf_path, maps, consensus_fasta):
    """
    Read an individual VCF, optionally liftover to consensus coordinates, and
    drop any position whose consensus base is N (species-discriminating site).

    Returns (samples, variants) where variants maps (chrom, start0) to a dict
    containing 'alleles' and per-sample GT tuples. Positions with allele
    conflicts across multiple records are silently dropped.
    """
    variants = {}
    with pysam.VariantFile(vcf_path) as vcf:
        samples = list(vcf.header.samples)
        for rec in vcf:
            if len(rec.alleles) != 2:
                continue
            ref, alt = rec.alleles
            if len(ref) != 1 or len(alt) != 1:
                continue

            chrom      = rec.chrom
            start0     = rec.start
            is_reverse = False

            if maps is not None:
                result = liftover_position(chrom, start0, maps)
                if result is None:
                    continue
                chrom, start0, is_reverse = result

            base = consensus_fasta.fetch(chrom, start0, start0 + 1).upper()
            if base == 'N':
                continue

            if is_reverse:
                ref = _comp(ref)
                alt = _comp(alt)

            key = (chrom, start0)
            if key not in variants:
                variants[key] = {'alleles': (ref, alt)}
            elif variants[key] is None:
                continue
            elif variants[key]['alleles'] != (ref, alt):
                variants[key] = None
                continue

            for sample in samples:
                gt = rec.samples[sample].get('GT')
                if gt is not None and None not in gt:
                    variants[key][sample] = gt

    return samples, {k: v for k, v in variants.items() if v is not None}


def main():
    parser = argparse.ArgumentParser(
        prog='consensus-demux',
        description=(
            'Build a joint species+individual demultiplexing VCF '
            'from consensus-aligned data.'))

    parser.add_argument('-c', '--consensus-vcf', required=True,
                        help='Species-level VCF from consensus-vcf (consensus coordinates)')
    parser.add_argument('-v', '--vcfs',           required=True,
                        help='Comma-separated individual VCFs, one per species')
    parser.add_argument('-s', '--species',        required=True,
                        help='Comma-separated species labels matching columns in the consensus VCF')
    parser.add_argument('-C', '--chains',         required=True,
                        help='Comma-separated chain files; use "none" for the reference species')
    parser.add_argument('-r', '--reference',      required=True,
                        help='Consensus genome FASTA (indexed with samtools faidx)')
    parser.add_argument('-o', '--output-prefix',  required=True,
                        help='Output prefix')
    parser.add_argument('--version', action='version',
                        version='CrossFilt v{version}'.format(version=__version__))

    args = parser.parse_args()

    vcf_files    = args.vcfs.split(',')
    species_list = args.species.split(',')
    chain_files  = args.chains.split(',')

    if len(vcf_files) != len(species_list) or len(chain_files) != len(species_list):
        sys.exit('Error: --vcfs, --species, and --chains must have the same number of entries.')

    consensus_fasta = pysam.Fastafile(args.reference)

    # Step 1: Load species-level VCF
    sp_samples, sp_variants = load_species_vcf(args.consensus_vcf)
    for sp in species_list:
        if sp not in sp_samples:
            sys.exit(f"Error: species '{sp}' not found in consensus VCF samples: {sp_samples}")

    # Step 2: Load individual VCFs, liftover non-reference species
    donors   = []   # list of (species_label, sample_id)
    ind_vars = {}   # (chrom, start0) -> {'alleles': ..., (species, sample): gt, ...}

    for species, vcf_path, chain_path in zip(species_list, vcf_files, chain_files):
        if chain_path.lower() == 'none':
            maps = None
        else:
            target_contigs = _vcf_contigs(vcf_path)
            maps, _, _ = read_chain_file(chain_path, target_contigs,
                                         list(consensus_fasta.references))

        samples, variants = load_individual_vcf(vcf_path, maps, consensus_fasta)

        for sample in samples:
            donors.append((species, sample))

        for key, data in variants.items():
            if key not in ind_vars:
                ind_vars[key] = {'alleles': data['alleles']}
            elif ind_vars[key] is None:
                continue
            elif ind_vars[key]['alleles'] != data['alleles']:
                ind_vars[key] = None
                continue
            for sample in samples:
                if sample in data:
                    ind_vars[key][(species, sample)] = data[sample]

    ind_vars = {k: v for k, v in ind_vars.items() if v is not None}

    # Step 3: Write combined VCF
    donor_ids   = [f"{sp}__{samp}" for sp, samp in donors]
    chrom_order = {c: i for i, c in enumerate(consensus_fasta.references)}
    all_pos     = sorted(set(sp_variants) | set(ind_vars),
                         key=lambda k: (chrom_order.get(k[0], 999999), k[1]))

    vcfh = pysam.VariantHeader()
    for did in donor_ids:
        vcfh.add_sample(did)
    for chrom in consensus_fasta.references:
        vcfh.add_meta('contig', items=[('ID', chrom)])
    vcfh.add_meta('FORMAT', items=[('ID', 'GT'), ('Number', 1), ('Type', 'String'),
                                   ('Description', 'Genotype')])

    n_sp = n_ind = 0
    with pysam.VariantFile(args.output_prefix + '.vcf', 'w', header=vcfh) as out:
        for key in all_pos:
            chrom, start0 = key
            if key in sp_variants:
                entry = sp_variants[key]
                r = out.new_record(contig=chrom, start=start0, stop=start0 + 1,
                                   alleles=entry['alleles'])
                for sp, samp in donors:
                    idx = entry[sp]
                    r.samples[f"{sp}__{samp}"]['GT'] = (idx, idx)
                out.write(r)
                n_sp += 1
            else:
                entry = ind_vars[key]
                r = out.new_record(contig=chrom, start=start0, stop=start0 + 1,
                                   alleles=entry['alleles'])
                for sp, samp in donors:
                    r.samples[f"{sp}__{samp}"]['GT'] = entry.get((sp, samp), (0, 0))
                out.write(r)
                n_ind += 1

    print(f"Wrote {n_sp + n_ind} variants ({n_sp} species-level, {n_ind} individual-level) "
          f"for {len(donors)} donors.")


if __name__ == '__main__':
    main()
