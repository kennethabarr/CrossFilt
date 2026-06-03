"""
Core functions for lifting BAM alignments between genomes using UCSC chain files.

Terminology used throughout: "target" is the source genome (reads are currently
aligned to it); "query" is the destination genome (reads will be lifted to it).
This follows UCSC chain file convention, which is the inverse of how many tools
use these terms.
"""
import bx.intervals.intersection as bx
import sys
import pysam
import array
import bz2
import gzip
import heapq
import math
import numpy as np
from collections import defaultdict

MAX_LIFTED_SEQ_LEN = 500   # reads whose lifted sequence exceeds this are discarded
MAX_INSERT_SIZE    = 10000  # paired-end pairs whose insert size exceeds this are discarded

_complement_table = str.maketrans('ACGTNXacgtnx', 'TGCANXtgcanx')

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
    
  
def nopen(f, mode="rb"):
    """Open a file path, stdin/stdout, gzip, or bzip2 transparently."""
    if not isinstance(f, str):
        return f
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
        else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
        else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
        else open(f, mode)

def reader(fname):
    """Yield decoded and stripped text lines from any source supported by nopen."""
    for l in nopen(fname):
        yield l.decode('utf8').strip().replace("\r", "")
        
def bam_header_generator(orig_header, chrom_size, prog_name, prog_ver, co, format_ver=1.0, sort_type = 'coordinate'):
    """
    Build a new BAM header with updated chromosome sizes and program metadata.

    Parameters
    ----------
    orig_header : dict
        Original BAM header from pysam header.to_dict().
    chrom_size : dict
        Chromosome name → size for the query (destination) genome.
    prog_name : str
        Program name for the PG tag.
    prog_ver : str
        Program version for the PG tag.
    co : list of str
        Comments to append to the CO tag.
    format_ver : str or float
        BAM format version for the HD VN field.
    sort_type : str
        Sort order for the HD SO field.

    Returns
    -------
    tuple of (dict, dict)
        Updated header dict and a chromosome name → integer ID mapping.
    """
    bamHeaderLine=orig_header.copy()
    name2id={}
    chrom_id = 0
    # replace 'HD'
    bamHeaderLine['HD'] = {'VN':format_ver,'SO':sort_type}

    # replace SQ
    tmp=[]
    for ref_name in sorted(chrom_size):
        tmp.append({'LN':chrom_size[ref_name],'SN':ref_name})
        name2id[ref_name] = chrom_id
        chrom_id += 1
    bamHeaderLine['SQ'] =  tmp
    if 'PG' in bamHeaderLine:
        bamHeaderLine['PG'] .append( {'ID':prog_name,'VN':prog_ver})
    else:
        bamHeaderLine['PG'] = [{'ID':prog_name,'VN':prog_ver}]

    for comment in co:
        if 'CO' in bamHeaderLine:
            bamHeaderLine['CO'].append(comment)
        else:
            bamHeaderLine['CO'] = [comment]
    return (bamHeaderLine, name2id)

def revcomp_DNA(dna):
    """Return the reverse complement of a DNA string (supports A/C/G/T/N/X)."""
    return dna[::-1].translate(_complement_table)

def read_chain_file(chain_file,target_contig_list, query_contig_list):
    """
    Parse a UCSC chain file into per-chromosome interval trees for fast lookup.

    Each node in the tree stores the corresponding query coordinates and strand.
    Builds one bx-python Intersecter per target chromosome.

    Parameters
    ----------
    chain_file : str
        Path to a chain file (plain, gzip, or bzip2).
    target_contig_list : list of str
        Contigs present in the target (source) genome.
    query_contig_list : list of str
        Contigs present in the query (destination) genome.

    Returns
    -------
    tuple of (dict, dict, dict)
        maps : target chrom → Intersecter of chain intervals
        tSizeDict : target chrom → chromosome size
        qSizeDict : query chrom → chromosome size
    """
    chainnames = ["score","tName","tSize","tStrand","tStart","tEnd","qName","qSize","qStrand","qStart","qEnd","id"]
    maps = {}
    this_chain = {}
    last_nfields = 1
    tSizeDict = {}
    qSizeDict = {}
    skip = False
    
    # Note target is the reference, query is the genome to map to. This terminology is confusing to me
    # and apparently also the writer of Crossmap, but I will try to be consistent
    
    for line in reader(chain_file):
        
        if not line.strip():
            continue
        sline=line.strip()
        if sline.startswith(('#',' ')):continue
        
        fields = line.rstrip().split()
        
        nfields = len(fields)
        
        if (last_nfields == 1 and fields[0] != 'chain'):
                raise Exception("Chain file has incorrect number of fields 1")

        if (last_nfields != 1 and fields[0] == 'chain'):
                raise Exception("Chain file has incorrect number of fields 2")
                                
        if fields[0] == 'chain' and nfields in [12, 13]:
            last_nfields = nfields
            # convert fields to the appropriate class and remove the 'chain' field
            fields = [t[0](t[1]) for t in zip([int, str, int, str, int, int, str, int, str, int, int, str], fields[1:])]
                
            # convert this to a dictionary
            this_chain = dict(zip(chainnames, fields))
            skip = False
            
            if this_chain['tName'] not in target_contig_list:
                skip = True

            if this_chain['qName'] not in query_contig_list:
                skip = True
                
            if skip: continue
            
            this_chain['mapTree'] = bx.Intersecter()
        
            if this_chain['tName'] not in maps:
                maps[this_chain['tName']] = bx.Intersecter()
                
            tfrom, qfrom = this_chain['tStart'], this_chain['qStart']
            
            tSizeDict[this_chain['tName']] = this_chain['tSize']
            qSizeDict[this_chain['qName']] = this_chain['qSize']
        
        elif (nfields == 3): # this is a data field 
            last_nfields = nfields
            if skip: continue
            size, tgap, qgap = int(fields[0]), int(fields[1]), int(fields[2])
            
            if this_chain['qStrand'] == '+':
                this_chain['mapTree'].add_interval( bx.Interval(tfrom, tfrom+size,(this_chain['qName'],qfrom, qfrom+size,this_chain['qStrand'])))
            elif this_chain['qStrand'] == '-':
                this_chain['mapTree'].add_interval( bx.Interval(tfrom, tfrom+size,(this_chain['qName'],this_chain['qSize'] - (qfrom+size), this_chain['qSize'] - qfrom, this_chain['qStrand'])))
                
            tfrom += size + tgap
            qfrom += size + qgap
        elif (nfields == 1): # this is a data field and the last in a chain
            last_nfields = nfields
            if skip: continue
            size = int(fields[0])

            if this_chain['qStrand'] == '+':
                this_chain['mapTree'].add_interval( bx.Interval(tfrom, tfrom+size,(this_chain['qName'],qfrom, qfrom+size,this_chain['qStrand'])))
            elif this_chain['qStrand'] == '-':
                this_chain['mapTree'].add_interval( bx.Interval(tfrom, tfrom+size,(this_chain['qName'],this_chain['qSize'] - (qfrom+size), this_chain['qSize'] - qfrom, this_chain['qStrand'])))
            maps[this_chain['tName']].add_interval(bx.Interval(this_chain['tStart'],this_chain['tEnd'], this_chain))
            
        else:
            raise Exception("Invalid chain format. (%s)" % line)
    return (maps, tSizeDict, qSizeDict)
        
def inside(s1,e1,s2,e2):
    """Return True if [s1, e1) is fully contained within [s2, e2)."""
    if (s1 >= s2 and e1 <= e2):
        return True
    else:
        return False
    
def get_chr_chains(maps, chrom):
    """Return the Intersecter for chrom, or None if chrom has no chains."""
    if (chrom not in maps):
        return None
    else:
        return maps[chrom]

def get_chains(chr_chains, start, end):
    """
    Return chains that fully contain [start, end), sorted by descending score.

    Parameters
    ----------
    chr_chains : bx.intervals.intersection.Intersecter
    start, end : int
        Read coordinates in target genome (0-based, half-open).
    """
    out = [c for c in chr_chains.find(start, end)
           if inside(start, end, c.value['tStart'], c.value['tEnd'])]
    if len(out) <= 1:
        return out
    return sorted(out, key=lambda c: -c.value['score'])
  
def get_chains_pe(chr_chains, start1, end1, start2, end2):
    """
    Return chains that fully contain both read1 [start1, end1) and read2 [start2, end2).

    Both reads must be fully enclosed in the same chain for a pair to be lifted.
    """
    out = [c for c in chr_chains.find(start1, end1)
           if inside(start1, end1, c.value['tStart'], c.value['tEnd'])
           and inside(start2, end2, c.value['tStart'], c.value['tEnd'])]
    return sorted(out, key=lambda c: -c.value['score'])

def add_solid_interval(out, read_chr, intervals, this_absolute_start, this_absolute_end,
                       this_relative_start, this_relative_end, this_add, target_fasta,
                       query_fasta, read_seq, tup, read_quality):
    """
    Append a single uninterrupted alignment block to the liftover output dict.

    Fetches target and query sequences for the block and substitutes read-specific
    variants: positions where the read differs from the target are kept as-is;
    positions that match the target are replaced by the query base.
    """
    query_chr  = intervals[0].value[0]
    # grab this section of read from the genome so we can find mismatches
    target_tmp    = target_fasta.fetch(read_chr, this_absolute_start, this_absolute_end).upper()
    
    # the query start is the interval plus the offset
    offset = abs(intervals[0].start - this_absolute_start)
    if (out['is_reverse']):
        query_start = intervals[0].value[2] - offset - this_add
    else:
        query_start = intervals[0].value[1] + offset

    out['segments'].append((query_start,query_start+this_add))
    
    if out['query_pos'] is None:
        out['query_pos'] = query_start
    else:
        out['query_pos'] = min(out['query_pos'], query_start)
                                
    query_tmp = query_fasta.fetch(query_chr, query_start, query_start+this_add).upper()
    read_add  = read_seq[this_relative_start:this_relative_end]

    t = np.frombuffer(target_tmp.encode(), dtype=np.uint8)
    q = np.frombuffer(query_tmp.encode(),  dtype=np.uint8).copy()
    r = np.frombuffer(read_add.encode(),   dtype=np.uint8)
    q[t != r] = r[t != r]
    query_add = q.tobytes().decode()

    if out['is_reverse']:
        query_add = revcomp_DNA(query_add)
    
    out['query_sequence'] += query_add
    out['cigartuples'].append(tup)
    out['qualityscores'].extend(read_quality[this_relative_start:this_relative_end])
    
    this_absolute_start += this_add
    this_relative_start += this_add 
    
    return out, this_absolute_start, this_relative_start

def add_gapped_interval(out, read_chr, intervals, this_absolute_start, this_absolute_end,
                        this_relative_start, this_relative_end, this_add, target_fasta,
                        query_fasta, read_seq, tup, read_quality):
    """
    Append an alignment block that spans multiple chain sub-intervals (indels between species).

    Called when a read's aligned segment crosses one or more gaps in the chain —
    i.e., insertions or deletions in the target relative to the query. Handles
    coordinate remapping, sequence conversion, and quality tracking across all
    sub-intervals and the gaps between them.
    """
    query_chr     = intervals[0].value[0]
    nintervals    = len(intervals)
    is_reverse    = out['is_reverse']
    new_qualities = array.array('B')
    # we dont want to modify the intervals in the intervaltree, 
    # so we need to copy the values to a mutable object
    target_ranges = []
    query_ranges  = []
    for i in range(nintervals):
        target_ranges.append([])
        query_ranges.append([])
        target_ranges[i].append(intervals[i].start)
        target_ranges[i].append(intervals[i].end)
        query_ranges[i].append(intervals[i].value[1])
        query_ranges[i].append(intervals[i].value[2]) 
          
    # Now we want to trim the intervals to match the coverage of the read and calculate interval lengths
    
    # trim the first range to the span here
    offset = abs(this_absolute_start - intervals[0].start)
    target_ranges[0][0] += offset
    if is_reverse:
        query_ranges[0][1] -= offset
    else:
        query_ranges[0][0] += offset
    # trim the last range to the span here
    offset = abs(this_absolute_end - intervals[-1].end)
    target_ranges[-1][1] -= offset
    if is_reverse:
        query_ranges[-1][0] += offset
    else:
        query_ranges[-1][1] -= offset
        
    # calculate interval lengths
    for i in range(nintervals):
        target_ranges[i].append(target_ranges[i][1]-target_ranges[i][0])
        query_ranges[i].append(query_ranges[i][1]-query_ranges[i][0])

    # find the first position in the query
    if is_reverse:
        query_start = query_ranges[-1][0]
        query_end   = query_ranges[0][1]
    else:
        query_start = query_ranges[0][0]
        query_end   = query_ranges[-1][1]
    
    # update the start position of this read
    if out['query_pos'] is None:
        out['query_pos'] = query_start
    else:
        out['query_pos'] = min(out['query_pos'], query_start)
    
    out['segments'].append((query_start,query_end))
    
    query_len  = query_end - query_start
        
    # how many insertions or deletions are there between ranges?
    # an insertion in the query shows up as a gap in start/end of query sequences
    # a deltion shows up as a gap in start/end of target sequences

    insertions = []
    deletions = []
    for i in range(nintervals-1):
        deletions.append(target_ranges[i+1][0] - target_ranges[i][1])
        if is_reverse:
            insertions.append(query_ranges[i][0]-query_ranges[i+1][1])
        else: 
            insertions.append(query_ranges[i+1][0]-query_ranges[i][1])

    query_parts = []
    tmp_relative_start = this_relative_start

    # now we can add to the sequence and qualities for all but the last interval
    for i in range(nintervals-1):
        tmp_relative_end    = tmp_relative_start+target_ranges[i][2]

        target_tmp = target_fasta.fetch(read_chr, target_ranges[i][0], target_ranges[i][1]).upper()
        query_tmp  = query_fasta.fetch(query_chr, query_ranges[i][0], query_ranges[i][1]).upper()
        read_add   = read_seq[tmp_relative_start:tmp_relative_end]

        new_qualities += read_quality[tmp_relative_start:tmp_relative_end]
        last_qual = read_quality[tmp_relative_end]

        t = np.frombuffer(target_tmp.encode(), dtype=np.uint8)
        q = np.frombuffer(query_tmp.encode(),  dtype=np.uint8).copy()
        r = np.frombuffer(read_add.encode(),   dtype=np.uint8)
        q[t != r] = r[t != r]
        add_tmp = q.tobytes().decode()

        query_parts.append(revcomp_DNA(add_tmp) if is_reverse else add_tmp)

        if (insertions[i] > 0):
            out['has_insertion'] = True
            if is_reverse:
                query_parts.append(revcomp_DNA(query_fasta.fetch(query_chr, query_ranges[i][0]-insertions[i], query_ranges[i][0]).upper()))
            else:
                query_parts.append(query_fasta.fetch(query_chr, query_ranges[i][1], query_ranges[i][1]+insertions[i]).upper())
            new_qualities += array.array('B', [last_qual]*insertions[i])

        tmp_relative_start += target_ranges[i][2]

        if (deletions[i] > 0):
            tmp_relative_start += deletions[i]
            out['has_deletion'] = True

    # finally add the last interval
    tmp_relative_end    = tmp_relative_start+target_ranges[-1][2]
    target_tmp = target_fasta.fetch(read_chr, target_ranges[-1][0], target_ranges[-1][1]).upper()
    query_tmp  = query_fasta.fetch(query_chr, query_ranges[-1][0], query_ranges[-1][1]).upper()
    read_add   = read_seq[tmp_relative_start:tmp_relative_end]

    t = np.frombuffer(target_tmp.encode(), dtype=np.uint8)
    q = np.frombuffer(query_tmp.encode(),  dtype=np.uint8).copy()
    r = np.frombuffer(read_add.encode(),   dtype=np.uint8)
    q[t != r] = r[t != r]
    add_tmp = q.tobytes().decode()

    new_qualities += read_quality[tmp_relative_start:tmp_relative_end]
    query_parts.append(revcomp_DNA(add_tmp) if is_reverse else add_tmp)

    out['query_sequence'] += ''.join(query_parts)
    out['cigartuples'].append((0, query_len))
    out['qualityscores'].extend(new_qualities)
    
    this_absolute_start += this_add
    this_relative_start += this_add 
    
    return out, this_absolute_start, this_relative_start



# error codes
# 0 = no error
# 1 = No chain contains target sequence (this gets called before this function)
# 2 = No match found
# 3 = The start or end position of the read is an insertion in target relative to query
# 4 = There are internal insertions or deletions in the target relative to the query

def _coords_only_gapped_interval(out, intervals, this_absolute_start, this_absolute_end, this_add):
    """
    Compute lifted coordinates for a gapped interval without any FASTA access.

    Used by liftover_segment when convert_seq=False. Updates out['segments'],
    out['query_pos'], and out['cigartuples'] only.
    """
    is_reverse = out['is_reverse']
    t_starts = [iv.start      for iv in intervals]
    t_ends   = [iv.end        for iv in intervals]
    q_starts = [iv.value[1]   for iv in intervals]
    q_ends   = [iv.value[2]   for iv in intervals]

    # trim first and last intervals to the read boundaries
    offset = abs(this_absolute_start - intervals[0].start)
    t_starts[0] += offset
    if is_reverse: q_ends[0]    -= offset
    else:          q_starts[0]  += offset

    offset = abs(this_absolute_end - intervals[-1].end)
    t_ends[-1] -= offset
    if is_reverse: q_starts[-1] += offset
    else:          q_ends[-1]   -= offset

    if is_reverse:
        query_start = q_starts[-1]
        query_end   = q_ends[0]
    else:
        query_start = q_starts[0]
        query_end   = q_ends[-1]

    out['segments'].append((query_start, query_end))
    out['query_pos'] = query_start if out['query_pos'] is None else min(out['query_pos'], query_start)
    # Use query_end - query_start (correct query genome M length). If this differs from
    # this_add (the original M length), _build_lifted_alignment will pad/trim the
    # sequence with N's to keep the BAM valid.
    out['cigartuples'].append((0, query_end - query_start))
    out['has_indel'] = True

    return out, this_absolute_start + this_add


def liftover_segment(chain, old_alignment, target_fasta, query_fasta, read_chr,
                     convert_seq=True):
    """
    Convert a single BAM alignment from target to query genome coordinates.

    Walks the CIGAR string and maps each match block through the chain interval
    tree. Returns a result dict; check result['pass'] before using it.

    Parameters
    ----------
    chain : bx.intervals.intersection.Interval
        Chain interval that fully contains this read.
    old_alignment : pysam.AlignedSegment
        Original alignment in the target genome.
    target_fasta : pysam.FastaFile
        Indexed FASTA for the target (source) genome.
    query_fasta : pysam.FastaFile
        Indexed FASTA for the query (destination) genome.
    read_chr : str
        Target chromosome name for the read.

    Returns
    -------
    dict
        Keys: 'pass' (bool), 'error type' (0=ok, 2=no chain match,
        3=read boundary falls in an indel), 'query_sequence', 'query_chrom',
        'query_pos', 'cigartuples', 'qualityscores', 'is_reverse'.
    """
    read_start   = old_alignment.reference_start
    read_seq     = old_alignment.query_sequence     # only used when convert_seq=True
    read_quality = old_alignment.query_qualities    # only used when convert_seq=True
    cigar_tuples = old_alignment.cigartuples
   
    out       = {}
    intervals = []
                      
    this_absolute_start = read_start
    this_relative_start = 0
    
    out['query_sequence']     = ''
    out['segments']           = []
    out['query_chrom']        = chain.value['qName']
    out['query_pos']          = None 
    out['cigartuples']        = []
    out['qualityscores']      = array.array('B')
    out['pass']               = True
    out['is_reverse']         = True if chain.value['qStrand'] == "-" else False
    
    out['error type']         = 0
    out['has_indel']          = False
    out['has_insertion']      = False
    out['has_deletion']       = False
    
    for tup in cigar_tuples:
        op = tup[0]
        if op == 0:  # match/mismatch
            this_add = tup[1]

            if not out['pass']: continue

            this_absolute_end = this_absolute_start + this_add
            this_relative_end = this_relative_start + this_add

            intervals = chain.value['mapTree'].find(this_absolute_start, this_absolute_end)

            if (len(intervals) == 0):
                out['pass'] = False
                out['error type'] = 2
                continue
            elif (len(intervals) == 1):
                # I will require the start and end position of the read to match in the query
                if (intervals[0].start > this_absolute_start or intervals[0].end < this_absolute_end):
                    out['pass'] = False
                    out['error type'] = 3
                    continue
                elif convert_seq:
                    out, this_absolute_start, this_relative_start = add_solid_interval(out,
                            read_chr, intervals, this_absolute_start, this_absolute_end,
                            this_relative_start, this_relative_end, this_add, target_fasta,
                            query_fasta, read_seq, tup, read_quality)
                else:
                    # coordinate-only: compute query position without FASTA access
                    offset = abs(intervals[0].start - this_absolute_start)
                    query_start = (intervals[0].value[2] - offset - this_add if out['is_reverse']
                                   else intervals[0].value[1] + offset)
                    out['segments'].append((query_start, query_start + this_add))
                    out['query_pos'] = (query_start if out['query_pos'] is None
                                        else min(out['query_pos'], query_start))
                    out['cigartuples'].append(tup)
                    this_absolute_start += this_add
                    this_relative_start += this_add

            else:  # we have gaps in the alignment
                out['has_indel'] = True
                if (intervals[0].start > this_absolute_start or intervals[-1].end < this_absolute_end):
                    out['pass'] = False
                    out['error type'] = 3
                    continue
                elif convert_seq:
                    out, this_absolute_start, this_relative_start = add_gapped_interval(out,
                        read_chr, intervals, this_absolute_start, this_absolute_end,
                        this_relative_start, this_relative_end, this_add, target_fasta,
                        query_fasta, read_seq, tup, read_quality)
                else:
                    out, this_absolute_start = _coords_only_gapped_interval(
                        out, intervals, this_absolute_start, this_absolute_end, this_add)
                    this_relative_start += this_add

        elif op in (1, 4):  # insertion or soft clip
            if not out['pass']: continue
            this_add = tup[1]
            this_relative_end = this_relative_start + this_add
            if convert_seq:
                out['query_sequence'] += read_seq[this_relative_start:this_relative_end]
                out['qualityscores'].extend(read_quality[this_relative_start:this_relative_end])
            out['cigartuples'].append(tup)
            out['segments'].append(None)
            this_relative_start = this_relative_end

        elif op in (2, 3):  # deletion or skip
            if not out['pass']: continue
            this_absolute_start += tup[1]
            out['cigartuples'].append(tup)
            out['segments'].append(None)
            
    
    # update the length of Ns
    if out['pass']:
        if convert_seq and (len(out['query_sequence']) > MAX_LIFTED_SEQ_LEN):
          out['pass'] = False
          out['error type'] = 2
        for i in range(len(out['cigartuples'])):
            if out['cigartuples'][i][0] == 3:
                if out['is_reverse']:
                    out['cigartuples'][i] = (3,out['segments'][i-1][0] - out['segments'][i+1][1])
                    if out['cigartuples'][i][1] < 0: raise Exception("Splice distance cannot be negative")
                else:
                    out['cigartuples'][i] = (3,out['segments'][i+1][0] - out['segments'][i-1][1])
                    if out['cigartuples'][i][1] < 0: raise Exception("Splice distance cannot be negative")
      
    return out

def _build_lifted_alignment(new_header, old_alignment, new_read, name_to_id,
                             convert_seq=True):
    """
    Construct a pysam AlignedRead in the query genome from a liftover result.

    Copies query name and tags from old_alignment, sets coordinates from
    new_read, handles strand reversal, and preserves the RG tag.

    When convert_seq=False the original sequence and quality scores are kept
    (only coordinates, orientation flag, and CIGAR are updated).
    """
    aln = pysam.AlignedRead(new_header)

    aln.query_name = old_alignment.query_name

    aln.next_reference_id    = -1
    aln.next_reference_start = 0
    aln.template_length      = 0

    aln.reference_id    = name_to_id[new_read['query_chrom']]
    aln.reference_start = new_read['query_pos']
    aln.mapping_quality = old_alignment.mapping_quality

    aln.flag = 0x0
    if old_alignment.is_reverse != new_read['is_reverse']:
        aln.flag = aln.flag | 0x10

    if convert_seq:
        if new_read['is_reverse']:
            aln.query_sequence  = revcomp_DNA(new_read['query_sequence'])
            aln.query_qualities = new_read['qualityscores'][::-1]
            aln.cigartuples     = new_read['cigartuples'][::-1]
        else:
            aln.query_sequence  = new_read['query_sequence']
            aln.query_qualities = new_read['qualityscores']
            aln.cigartuples     = new_read['cigartuples']
    else:
        cigtuples = (new_read['cigartuples'][::-1] if new_read['is_reverse']
                     else new_read['cigartuples'])
        aln.cigartuples = cigtuples

        # Required sequence length from the new CIGAR (sum of M/I/S ops)
        cigar_seq_len = sum(l for op, l in cigtuples if op in (0, 1, 4))
        orig_seq  = old_alignment.query_sequence
        orig_qual = old_alignment.query_qualities

        if len(orig_seq) != cigar_seq_len:
            # Species-specific indel within a match block changed the M length.
            # Replace with N's so the BAM remains valid; the coordinate is correct.
            orig_seq  = 'N' * cigar_seq_len
            orig_qual = array.array('B', [0] * cigar_seq_len)

        if new_read['is_reverse'] != old_alignment.is_reverse:
            aln.query_sequence  = revcomp_DNA(orig_seq)
            aln.query_qualities = orig_qual[::-1]
        else:
            aln.query_sequence  = orig_seq
            aln.query_qualities = orig_qual

    # Copy all tags in one pass, fixing the RG value type inline to avoid
    # a second get_tag/set_tag round-trip.
    tags = old_alignment.get_tags(with_value_type=True)
    aln.set_tags([(t, str(v) if t == 'RG' else v, tp) for t, v, tp in tags])

    return aln


def process_se(SAMFILE        = None,
               outfile        = None,
               old_header     = None,
               new_header     = None,
               target_fasta   = None,
               query_fasta    = None,
               maps           = None,
               tSizes         = None,
               qSizes         = None,
               chainfile      = None,
               name_to_id     = None,
               best           = None,
               convert_seq    = True):
    """
    Lift all single-end reads in SAMFILE to the query genome and write to outfile.

    Parameters
    ----------
    SAMFILE : pysam.AlignmentFile
        Input BAM coordinate-sorted against the target genome.
    outfile : pysam.AlignmentFile
        Output BAM handle for successfully lifted reads.
    old_header : dict
        Original BAM header (used for reference name lookup).
    new_header : dict
        New BAM header with query genome chromosome sizes.
    target_fasta : pysam.FastaFile
        Target genome sequence.
    query_fasta : pysam.FastaFile
        Query genome sequence.
    maps : dict
        Per-chromosome interval trees from read_chain_file.
    tSizes, qSizes : dict
        Target / query chromosome size dicts.
    name_to_id : dict
        Query chromosome name → integer BAM reference ID.
    best : bool
        If True, only attempt liftover with the highest-scoring chain.

    Returns
    -------
    tuple of int
        (total_reads, lifted_ok, no_chain, no_match, boundary_indel, unpaired)
    """
    OUT_FILE_QUERY = outfile

    new_header = pysam.AlignmentHeader.from_dict(new_header)

    nreads = n0 = n1 = n2 = n3 = n4 = 0
    for old_alignment in SAMFILE.fetch():

        nreads += 1
        
        read_chr       = SAMFILE.get_reference_name(old_alignment.reference_id)
        read_start     = old_alignment.reference_start
        read_end       = old_alignment.reference_end
            
        chr_chains = get_chr_chains(maps, read_chr)
        if chr_chains is None: continue
        
        chains = get_chains(chr_chains, read_start, read_end)
        nchains = len(chains)
        if (nchains == 0): 
            n1 += 1
            continue
        
        if best: 
          chains  = [chains[0]]
          nchains = 1
        
        error_type = 0
        for i in range(nchains):
          new_read = liftover_segment(chains[i], old_alignment, target_fasta, query_fasta, read_chr,
                                       convert_seq=convert_seq)
          
          if i == 0:
            error_type = new_read['error type']
            
          if new_read['pass']:
            error_type = new_read['error type']
            break
          
        if (error_type == 0):
            n0 += 1
        elif (error_type == 2):
            n2 += 1
            continue
        elif (error_type == 3):
            n3 += 1
            continue


        OUT_FILE_QUERY.write(_build_lifted_alignment(new_header, old_alignment, new_read, name_to_id,
                                                      convert_seq=convert_seq))
    
    return (nreads, n0, n1, n2, n3, n4)
  
def process_pe(SAMFILE        = None,
               outfile        = None,
               old_header     = None,
               new_header     = None,
               target_fasta   = None,
               query_fasta    = None,
               maps           = None,
               tSizes         = None,
               qSizes         = None,
               chainfile      = None,
               name_to_id     = None,
               best           = None,
               convert_seq    = True):
    """
    Lift all paired-end reads in SAMFILE to the query genome and write to outfile.

    Both reads in a pair must be lifted successfully by the same chain;
    pairs that end up on the same strand or exceed 10 kb insert size are discarded.

    Parameters
    ----------
    SAMFILE : pysam.AlignmentFile
        Input BAM coordinate-sorted against the target genome.
    outfile : pysam.AlignmentFile
        Output BAM handle for successfully lifted read pairs.
    old_header : dict
        Original BAM header (used for reference name lookup).
    new_header : dict
        New BAM header with query genome chromosome sizes.
    target_fasta : pysam.FastaFile
        Target genome sequence.
    query_fasta : pysam.FastaFile
        Query genome sequence.
    maps : dict
        Per-chromosome interval trees from read_chain_file.
    tSizes, qSizes : dict
        Target / query chromosome size dicts.
    name_to_id : dict
        Query chromosome name → integer BAM reference ID.
    best : bool
        If True, only attempt liftover with the highest-scoring chain.

    Returns
    -------
    tuple of int
        (total_reads, lifted_ok, no_chain, no_match, boundary_indel, failed_pair)
    """
    new_header = pysam.AlignmentHeader.from_dict(new_header)
    OUT_FILE_QUERY = outfile

    nreads = n0 = n1 = n2 = n3 = n4 = 0
    for old1, old2 in read_pair_generator(SAMFILE):
        nreads += 2

        read1_chr       = SAMFILE.get_reference_name(old1.reference_id)
        read1_start     = old1.reference_start
        read1_end       = old1.reference_end
        
        read2_chr       = SAMFILE.get_reference_name(old2.reference_id)
        read2_start     = old2.reference_start
        read2_end       = old2.reference_end
            
        if read1_chr != read2_chr: 
          n4 += 2
          continue
        
        chr_chains = get_chr_chains(maps, read1_chr)
        if chr_chains is None: 
          n4 += 2
          continue
        
        chains = get_chains_pe(chr_chains, read1_start, read1_end, read2_start, read2_end)
        nchains = len(chains)
        
        if (nchains == 0): 
            n1 += 2
            continue
        
        if best: 
          chains  = [chains[0]]
          nchains = 1
          
        error_type1 = 0
        error_type2 = 0
        for i in range(nchains):
          new1 = liftover_segment(chains[i], old1, target_fasta, query_fasta, read1_chr,
                                   convert_seq=convert_seq)
          new2 = liftover_segment(chains[i], old2, target_fasta, query_fasta, read2_chr,
                                   convert_seq=convert_seq)
        
          if i == 0:
            error_type1 = new1['error type']
            error_type2 = new2['error type']
          
          if new1['pass'] and new2['pass']:
            error_type1 = new1['error type']
            error_type2 = new2['error type']
            break

        if (error_type1 == 2) or (error_type2 == 2):
            n2 += 2
            continue
        elif (error_type1 == 3) or (error_type2 == 3):
            n3 += 2
            continue
        
        new_alignment1 = _build_lifted_alignment(new_header, old1, new1, name_to_id,
                                                  convert_seq=convert_seq)
        new_alignment2 = _build_lifted_alignment(new_header, old2, new2, name_to_id,
                                                  convert_seq=convert_seq)

        # make sure forward and reverse orientation are preserved
        if new_alignment1.is_reverse == new_alignment2.is_reverse:
          n4 += 2
          continue
        
        new_tlen = 0
        if new_alignment1.is_reverse:
          new_tlen = new_alignment2.reference_start - new_alignment1.reference_end
          new_alignment2.mate_is_reverse = True
        else:
          new_tlen = new_alignment2.reference_end - new_alignment1.reference_start
          new_alignment1.mate_is_reverse = True
        
        if (abs(new_tlen) > MAX_INSERT_SIZE):
          n4 += 2
          continue
        
        n0 += 2
        new_alignment1.template_length = new_tlen
        new_alignment2.template_length = new_tlen
        
        new_alignment1.is_paired = True
        new_alignment2.is_paired = True
        
        new_alignment1.is_proper_pair = True
        new_alignment2.is_proper_pair = True
        
        new_alignment1.is_read1 = True
        new_alignment2.is_read2 = True
        
        new_alignment1.mate_is_unmapped = False
        new_alignment2.mate_is_unmapped = False

        OUT_FILE_QUERY.write(new_alignment1)
        OUT_FILE_QUERY.write(new_alignment2)

    return (nreads, n0, n1, n2, n3, n4)


def pack_chromosomes(index_stats, n_files):
    """
    Distribute chromosomes into n_files bins as evenly as possible by read count.

    Uses the LPT (Longest Processing Time First) greedy algorithm: sort
    chromosomes by read count descending, then assign each one to the bin
    with the current smallest total.  This produces near-optimal balance
    without solving the NP-hard bin-packing problem exactly.

    Parameters
    ----------
    index_stats : list
        Sequence of pysam index-statistics objects (from
        AlignmentFile.get_index_statistics()).  Each element must expose
        ``.contig`` (str) and ``.total`` (int) attributes.
    n_files : int
        Desired number of output bins.

    Returns
    -------
    list of list of str
        Each inner list contains the chromosome names assigned to that bin,
        in the order they were assigned.  Bins are ordered by their internal
        index (stable assignment order).  Fewer than *n_files* bins are
        returned when there are fewer non-empty chromosomes than *n_files*.
    """
    chroms = [(s.total, s.contig) for s in index_stats if s.total > 0]
    if not chroms:
        return []

    n = min(n_files, len(chroms))
    chroms.sort(reverse=True)           # LPT: assign largest chromosomes first

    # heap entries: (bin_total, bin_index, chrom_list)
    heap = [(0, i, []) for i in range(n)]
    heapq.heapify(heap)

    for count, chrom in chroms:
        total, idx, chrom_list = heapq.heappop(heap)
        chrom_list.append(chrom)
        heapq.heappush(heap, (total + count, idx, chrom_list))

    heap.sort(key=lambda x: x[1])      # restore stable bin-index order
    return [cl for _, _, cl in heap if cl]
