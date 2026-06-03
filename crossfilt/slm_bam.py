#!/usr/bin/python3
"""
crossfilt-slm: Split, lift, and merge a BAM file between genomes.

Splits the input BAM into N chromosome-balanced chunks using the BAM index
(no per-read iteration during splitting), lifts each chunk in a separate
process (so chain parsing and liftover run in parallel), then merges the
results into a single coordinate-sorted output.

Each worker reads only the chain entries relevant to its chromosomes, keeping
both memory and chain-parse time proportional to its share of the data.
"""

import sys
import os
import argparse
import logging
import multiprocessing
import pysam
import tempfile
from timeit import default_timer as timer
from .lib import liftover_functions as lift
import importlib.metadata

__version__ = importlib.metadata.version('crossfilt')

logger = logging.getLogger(__name__)


def _chunk_paths(out_prefix, bin_idx):
    """Return the (chunk_path, lifted_path) pair for a given bin index."""
    return (f"{out_prefix}.chunk.{bin_idx}.bam",
            f"{out_prefix}.chunk.{bin_idx}.lifted.bam")


def _remove_if_exists(path):
    """Remove a file silently if it exists."""
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


# ------------------------------------------------------------------ #
#  Top-level worker (must be module-level to be picklable)            #
# ------------------------------------------------------------------ #

def _lift_chunk(args):
    """
    Extract, lift, and sort one chromosome bin.

    Runs in a worker process spawned by multiprocessing.Pool.  All arguments
    are passed as a plain dict so that adding fields never silently misaligns
    positional unpacking.

    Parameters
    ----------
    args : dict
        Keys: bin_idx, chroms, chrom_ord, infile, out_prefix,
        chainfile, target_fasta, query_fasta, old_header, new_header,
        name_to_id, is_paired, best, convert_seq, bam_threads, sort_threads,
        keep_chunks

    Thread allocation
    -----------------
    bam_threads : int
        BAM I/O threads for reading/writing during the lift phase.  Should be
        1 when N workers run simultaneously — liftover is CPU-bound Python and
        extra I/O threads only waste cores.
    sort_threads : int
        Threads passed to ``pysam.sort`` inside each worker.  Set to
        ``total_threads // n_workers`` so the parallel sorts together consume
        roughly the full thread budget.

    Returns
    -------
    dict
        bin_idx, lifted_path, n_total, n_lifted, n_no_chain, n_no_match,
        n_boundary_indel, n_unpaired
    """
    bin_idx      = args['bin_idx']
    chroms       = args['chroms']
    chrom_ord    = args['chrom_ord']
    infile       = args['infile']
    out_prefix   = args['out_prefix']
    chainfile    = args['chainfile']
    target_fasta = args['target_fasta']
    query_fasta  = args['query_fasta']
    old_header   = args['old_header']
    new_header   = args['new_header']
    name_to_id   = args['name_to_id']
    is_paired    = args['is_paired']
    best         = args['best']
    convert_seq  = args['convert_seq']
    bam_threads  = args['bam_threads']
    sort_threads = args['sort_threads']
    keep_chunks  = args['keep_chunks']

    chunk_path, lifted_path = _chunk_paths(out_prefix, bin_idx)

    # -- Extract chromosomes for this bin (already in header order) ------
    sam = pysam.AlignmentFile(infile, "rb", threads=bam_threads)
    ordered = sorted(chroms, key=lambda c: chrom_ord.get(c, float('inf')))
    chunk_out = pysam.AlignmentFile(chunk_path, "wb", header=old_header,
                                    threads=bam_threads,
                                    format_options=[b'level=1'])
    for chrom in ordered:
        for read in sam.fetch(chrom):
            chunk_out.write(read)
    chunk_out.close()
    sam.close()
    pysam.index(chunk_path)

    # -- Read chain file filtered to this bin's chromosomes only ---------
    TARGETFILE = pysam.Fastafile(target_fasta)
    QUERYFILE  = pysam.Fastafile(query_fasta)
    maps, tSizes, _ = lift.read_chain_file(chainfile, chroms,
                                            QUERYFILE.references)

    # -- Lift (temp file guarded so it is always cleaned up) -------------
    # bam_threads=1: liftover is CPU-bound Python; extra I/O threads waste cores
    chunk_sam = pysam.AlignmentFile(chunk_path, "rb", threads=bam_threads)
    _fd, tmp_path = tempfile.mkstemp(suffix='.bam')
    os.close(_fd)
    try:
        tmp_bam = pysam.AlignmentFile(tmp_path, "wb", header=new_header,
                                      threads=bam_threads)
        kwds = dict(
            SAMFILE      = chunk_sam,
            outfile      = tmp_bam,
            old_header   = old_header,
            new_header   = new_header,
            target_fasta = TARGETFILE,
            query_fasta  = QUERYFILE,
            maps         = maps,
            tSizes       = tSizes,
            qSizes       = {},        # not used inside process_se / process_pe
            name_to_id   = name_to_id,
            best         = best,
            convert_seq  = convert_seq,
        )
        result = lift.process_pe(**kwds) if is_paired else lift.process_se(**kwds)
        tmp_bam.close()
    finally:
        chunk_sam.close()
        TARGETFILE.close()
        QUERYFILE.close()
        if not keep_chunks:
            _remove_if_exists(chunk_path)
            _remove_if_exists(chunk_path + ".bai")

    # -- Sort lifted chunk (sort_threads = total_threads // n_workers) ----
    pysam.sort("-@", str(sort_threads), "-o", lifted_path, tmp_path)
    _remove_if_exists(tmp_path)

    # Unpack result tuple by name so callers never use magic indices.
    # process_se / process_pe return (total, lifted_ok, no_chain, no_match,
    #                                  boundary_indel, unpaired)
    n_total, n_lifted, n_no_chain, n_no_match, n_boundary_indel, n_unpaired = result

    return dict(
        bin_idx          = bin_idx,
        lifted_path      = lifted_path,
        n_total          = n_total,
        n_lifted         = n_lifted,
        n_no_chain       = n_no_chain,
        n_no_match       = n_no_match,
        n_boundary_indel = n_boundary_indel,
        n_unpaired       = n_unpaired,
    )


# ------------------------------------------------------------------ #
#  Main entry point                                                   #
# ------------------------------------------------------------------ #

def main():
    parser = argparse.ArgumentParser(
        prog='crossfilt-slm',
        description=(
            'Split a BAM by chromosome bins, lift each bin in parallel, '
            'and merge the results.  Uses the BAM index for fast chromosome '
            'extraction; each worker reads only the chain entries for its '
            'chromosomes.'
        )
    )

    parser.add_argument("-i", "--input",         required=True,
                        help="Input BAM file (coordinate-sorted; indexed or indexable)")
    parser.add_argument("-o", "--output",        required=True,
                        help="Output BAM prefix (result written to <prefix>.bam)")
    parser.add_argument("-c", "--chain",         required=True,
                        help="UCSC chain file for the liftover")
    parser.add_argument("-t", "--target-fasta",  required=True,
                        help="FASTA for the target genome (reads are currently "
                             "aligned to this)")
    parser.add_argument("-q", "--query-fasta",   required=True,
                        help="FASTA for the query genome (reads will be lifted to "
                             "this)")
    parser.add_argument("-n", "--nfiles",        required=True, type=int,
                        help="Number of parallel chromosome chunks")
    parser.add_argument("-p", "--paired",        action="store_true",
                        help="Input contains paired-end reads")
    parser.add_argument("-b", "--best",          action="store_true",
                        help="Only use the highest-scoring chain for each read")
    parser.add_argument("-@", "--threads",       type=int, default=1,
                        help="Total threads available.  During the parallel "
                             "lift phase each worker uses 1 BAM I/O thread "
                             "(liftover is CPU-bound); sort threads per worker "
                             "= threads // n_chunks.  The final merge uses all "
                             "threads.")
    parser.add_argument("--no-seq",             action="store_true",
                        help="Skip sequence conversion; only update coordinates "
                             "and CIGAR")
    parser.add_argument("--keep-chunks",        action="store_true",
                        help="Keep intermediate per-chunk BAM files after merging")
    parser.add_argument('--version', action='version',
                        version='CrossFilt v{version}'.format(version=__version__))

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(message)s',
                        stream=sys.stderr)

    infile       = args.input
    out_prefix   = args.output
    chainfile    = args.chain
    target_fasta = args.target_fasta
    query_fasta  = args.query_fasta
    n_files      = args.nfiles
    is_paired    = args.paired
    best         = args.best
    convert_seq  = not args.no_seq
    threads      = args.threads
    keep_chunks  = args.keep_chunks

    # ------------------------------------------------------------------ #
    #  Index check and chromosome packing                                  #
    # ------------------------------------------------------------------ #
    if not os.path.exists(infile + ".bai"):
        logger.info("Indexing input BAM...")
        pysam.index(infile)

    sam        = pysam.AlignmentFile(infile, "rb", threads=threads)
    idx_stats  = sam.get_index_statistics()
    old_header = sam.header.to_dict()
    chrom_ord  = {sq['SN']: i for i, sq in enumerate(old_header['SQ'])}
    sam.close()

    # Validate contigs against target FASTA
    TARGETFILE = pysam.Fastafile(target_fasta)
    for stat in idx_stats:
        if stat.total > 0 and stat.contig not in TARGETFILE.references:
            TARGETFILE.close()
            sys.exit(f"Contig '{stat.contig}' not found in target FASTA. "
                     "Did you use the right contig names?")
    TARGETFILE.close()

    t0 = timer()
    logger.info("Packing chromosomes into %d chunks...", n_files)
    chrom_bins = lift.pack_chromosomes(idx_stats, n_files)
    if not chrom_bins:
        sys.exit("No mapped reads found in BAM index.")

    # Warn if fewer bins were produced than requested (fewer chroms than n_files)
    if len(chrom_bins) < n_files:
        logger.info("  Note: only %d non-empty chromosome bin(s) available "
                    "(requested %d)", len(chrom_bins), n_files)

    nreads_by_chrom = {s.contig: s.total for s in idx_stats}
    for i, chroms in enumerate(chrom_bins):
        bin_total = sum(nreads_by_chrom.get(c, 0) for c in chroms)
        logger.info("  chunk %d: %d chrom(s), %s reads  [%s]",
                    i, len(chroms), f"{bin_total:,}", ', '.join(chroms))
    logger.info("  Done in %ss\n", round(timer() - t0, 2))

    # ------------------------------------------------------------------ #
    #  Build shared output BAM header from query FASTA index              #
    #  (done here so every worker gets the full set of query chromosomes) #
    # ------------------------------------------------------------------ #
    QUERYFILE  = pysam.Fastafile(query_fasta)
    qSizes     = {r: QUERYFILE.get_reference_length(r)
                  for r in QUERYFILE.references}
    QUERYFILE.close()

    comments = ['ORIGINAL_BAM_FILE=' + infile]
    new_header, name_to_id = lift.bam_header_generator(
        orig_header=old_header,
        chrom_size=qSizes,
        prog_name="CrossFilt",
        prog_ver=__version__,
        format_ver=__version__,
        sort_type='coordinate',
        co=comments,
    )

    # ------------------------------------------------------------------ #
    #  Dispatch workers                                                    #
    # ------------------------------------------------------------------ #
    n_workers    = len(chrom_bins)
    # BAM I/O during liftover: 1 thread per worker — liftover is CPU-bound
    # Python; extra I/O threads would just compete with each other.
    bam_threads  = 1
    # Sort after liftover: divide remaining threads evenly across workers.
    sort_threads = max(1, threads // n_workers)

    logger.info("Thread allocation: %d worker(s) × 1 BAM I/O thread + "
                "%d sort thread(s); %d thread(s) for final merge",
                n_workers, sort_threads, threads)

    worker_args = [
        dict(
            bin_idx      = bin_idx,
            chroms       = chroms,
            chrom_ord    = chrom_ord,
            infile       = infile,
            out_prefix   = out_prefix,
            chainfile    = chainfile,
            target_fasta = target_fasta,
            query_fasta  = query_fasta,
            old_header   = old_header,
            new_header   = new_header,
            name_to_id   = name_to_id,
            is_paired    = is_paired,
            best         = best,
            convert_seq  = convert_seq,
            bam_threads  = bam_threads,
            sort_threads = sort_threads,
            keep_chunks  = keep_chunks,
        )
        for bin_idx, chroms in enumerate(chrom_bins)
    ]

    logger.info("Launching %d parallel worker(s)...\n", n_workers)
    t_all = timer()

    # Wrap pool.map so that partial chunk files are cleaned up if a worker
    # raises before completing.
    try:
        with multiprocessing.Pool(processes=n_workers) as pool:
            chunk_results = pool.map(_lift_chunk, worker_args)
    except Exception:
        logger.info("Worker failure — cleaning up partial chunk files...")
        for idx in range(len(chrom_bins)):
            chunk_path, lifted_path = _chunk_paths(out_prefix, idx)
            for p in [chunk_path, chunk_path + ".bai", lifted_path]:
                _remove_if_exists(p)
        raise

    # sort by bin_idx so merge order is deterministic
    chunk_results.sort(key=lambda r: r['bin_idx'])

    # ------------------------------------------------------------------ #
    #  Summarise per-chunk results                                         #
    # ------------------------------------------------------------------ #
    lifted_paths    = []
    total_processed = total_lifted = 0

    for r in chunk_results:
        pct = round(100 * r['n_lifted'] / max(r['n_total'], 1), 2)
        logger.info("  chunk %d: %s/%s (%s%%) lifted",
                    r['bin_idx'],
                    f"{r['n_lifted']:,}", f"{r['n_total']:,}", pct)
        total_processed += r['n_total']
        total_lifted    += r['n_lifted']
        lifted_paths.append(r['lifted_path'])

    elapsed = round(timer() - t_all, 2)
    logger.info("\nAll chunks completed in %ss", elapsed)

    if total_lifted == 0:
        logger.info("Zero reads successfully lifted across all chunks.")
        for p in lifted_paths:
            _remove_if_exists(p)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    #  Merge                                                               #
    # ------------------------------------------------------------------ #
    t0 = timer()
    logger.info("Merging chunks...")
    final_bam = out_prefix + ".bam"

    if len(lifted_paths) == 1:
        os.rename(lifted_paths[0], final_bam)
    else:
        pysam.merge("-f", "-@", str(threads), final_bam, *lifted_paths)
        if not keep_chunks:
            for p in lifted_paths:
                _remove_if_exists(p)

    pysam.index("-@", str(threads), final_bam)
    logger.info("  Done in %ss", round(timer() - t0, 2))

    pct_total = round(100 * total_lifted / max(total_processed, 1), 2)
    logger.info("\nProcessed %s reads total; %s (%s%%) successfully lifted.",
                f"{total_processed:,}", f"{total_lifted:,}", pct_total)
    logger.info("Output: %s", final_bam)


if __name__ == '__main__':
    main()
