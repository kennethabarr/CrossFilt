import bx.intervals.intersection as bx
import sys
import pysam
import logging
import array
import bz2
import gzip
import urllib
from subprocess import Popen
import math
from collections import defaultdict

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
    if not isinstance(f, str):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r": return p.stdout
        return p
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
        else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
        else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
        else urllib.urlopen(f) if f.startswith(("http://", "https://","ftp://")) \
        else open(f, mode)

def reader(fname):
    for l in nopen(fname):
        yield l.decode('utf8').strip().replace("\r", "")
        
def bam_header_generator(orig_header, chrom_size, prog_name, prog_ver, co, format_ver=1.0, sort_type = 'coordinate'):
    '''
    generates header section for BAM file
    '''
    bamHeaderLine=orig_header.copy()
    name2id={}
    id = 0
    # replace 'HD'
    bamHeaderLine['HD'] = {'VN':format_ver,'SO':sort_type}

    # replace SQ
    tmp=[]
    for ref_name in sorted(chrom_size):
        tmp.append({'LN':chrom_size[ref_name],'SN':ref_name})
        name2id[ref_name] = id
        id += 1
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

complement = {'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
def revcomp_DNA(dna):
    return ''.join([complement[base] for base in reversed(dna)])

def read_chain_sizes(chain_file,target_contig_list, query_contig_list):
    chainnames = ["score","tName","tSize","tStrand","tStart","tEnd","qName","qSize","qStrand","qStart","qEnd","id"]
    last_nfields = 1
    tSizeDict = {}
    qSizeDict = {}
    skip = False
    
    # Note target is the reference, query is the genome to map to. This terminology is confusing to me
    # and apparently also the writer of Crossmap, but I will try to be consistent
    
    for line in reader(chain_file):
        
        if not line.strip(): continue
        sline=line.strip()
        if sline.startswith(('#',' ')): continue
        
        fields = line.rstrip().split()
        
        nfields = len(fields)
        
        if fields[0] == 'chain' and nfields in [12, 13]:
            # convert fields to the appropriate class and remove the 'chain' field
            fields = [t[0](t[1]) for t in zip([int, str, int, str, int, int, str, int, str, int, int, str], fields[1:])]
                
            # convert this to a dictionary
            this_chain = dict(zip(chainnames, fields))
            skip = False
            
            if this_chain['tName'] not in target_contig_list:
                skip = True
                if (this_chain['tName'].find("alt") == -1 and
                   this_chain['tName'].find("fix") == -1 and
                   this_chain['tName'].find("chrUn") == -1 and
                   this_chain['tName'].find("random") == -1):
                  print("Message: Contig " + this_chain['tName'] + " not in target file. Skipping chain", file=sys.stderr)
            
            if this_chain['qName'] not in query_contig_list:
                skip = True
                if (this_chain['qName'].find("alt") == -1 and
                   this_chain['qName'].find("fix") == -1 and
                   this_chain['qName'].find("chrUn") == -1 and
                   this_chain['qName'].find("random") == -1):
                  print("Message: Contig " + this_chain['qName'] + " not in query file. Skipping chain", file=sys.stderr)
                
            if skip: continue
            
            tSizeDict[this_chain['tName']] = this_chain['tSize']
            qSizeDict[this_chain['qName']] = this_chain['qSize']

    return (tSizeDict, qSizeDict)
  
def read_chain_file(chain_file,target_contig_list, query_contig_list):
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
                raise Exception("Chain file has incorrect number of fields 1", file=sys.stderr)
                
        if (last_nfields != 1 and fields[0] == 'chain'):
                raise Exception("Chain file has incorrect number of fields 2", file=sys.stderr)
                                
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
        
def inside(s1,e1,s2,e2): # is range 1 inside range 2?
    if (s1 >= s2 and e1 <= e2):
        return True
    else:
        return False
    
def get_chr_chains(maps, chrom): # we dont want to keep having to keep checking chr since the reads will be sorted by chr
    if (chrom not in maps):
        return None
    else:
        return maps[chrom]

def get_chains(chr_chains, start, end):
    chains = sorted(chr_chains.find(start, end), key=lambda chain: -chain.value['score'])
    out = []
    for chain in chains:
        if inside(start, end, chain.value['tStart'],chain.value['tEnd']):
            out.append(chain)
    return(out)
  
def get_chains_pe(chr_chains, start1, end1, start2, end2):
    chains = sorted(chr_chains.find(start1, end1), key=lambda chain: -chain.value['score'])
    out = []
    for chain in chains:
        if inside(start1, end1, chain.value['tStart'],chain.value['tEnd']):
          if inside(start2, end2, chain.value['tStart'],chain.value['tEnd']):
            out.append(chain)
    return(out)

def string_ident(str1, str2): 
    s = sum(1 for a, b in zip(str1, str2) if a != b)
    return s/len(str1)

def add_solid_interval(out, read_chr, intervals, this_absolute_start, this_absolute_end, 
                       this_relative_start, this_relative_end, this_add, target_fasta, 
                       query_fasta, read_seq, tup, read_quality):

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
    read_add = read_seq[this_relative_start:this_relative_end]
                
    query_add  = list(query_tmp)
    for a,b,j in zip(target_tmp, 
                     read_add,
                     range(len(target_tmp))):
      if a != b:
        query_add[j]  = read_add[j]
    query_add = ''.join(query_add)

    #query_add  = ''
    #for j in range(len(read_add)):
    #    if (target_tmp[j] == read_add[j]):
    #        query_add  += query_tmp[j]
    #    else:
    #        query_add  += read_add[j]
    
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

    query_add = ''
    tmp_relative_start = this_relative_start
        
    # now we can add to the sequence and qualities for all but the last interval
    for i in range(nintervals-1):
        #add_tmp             = ''
        tmp_relative_end    = tmp_relative_start+target_ranges[i][2]
    
        target_tmp = target_fasta.fetch(read_chr, target_ranges[i][0], target_ranges[i][1]).upper()
        query_tmp  = query_fasta.fetch(query_chr, query_ranges[i][0], query_ranges[i][1]).upper()
        read_add   = read_seq[tmp_relative_start:tmp_relative_end]

        new_qualities += read_quality[tmp_relative_start:tmp_relative_end] 
        last_qual = read_quality[tmp_relative_end]
        
        add_tmp  = list(query_tmp)
        for a,b,j in zip(target_tmp, 
                         read_add,
                         range(target_ranges[i][2])):
          if a != b:
            add_tmp[j]  = read_add[j]
        add_tmp = ''.join(add_tmp)

        #for j in range(target_ranges[i][2]):
        #    if (target_tmp[j] == read_add[j]):
        #        #target_add += target_tmp[j]
        #        add_tmp  += query_tmp[j]
        #    else:
        #        #target_add += read_add[j]
        #        add_tmp  += read_add[j]
            
        if is_reverse:
            query_add += revcomp_DNA(add_tmp)
        else:
            query_add += add_tmp

        if (insertions[i] > 0):
            out['has_insertion'] = True
            if is_reverse:
                query_add += revcomp_DNA(query_fasta.fetch(query_chr, query_ranges[i][0]-insertions[i], query_ranges[i][0]).upper())
            else:
                query_add += query_fasta.fetch(query_chr, query_ranges[i][1], query_ranges[i][1]+insertions[i]).upper()
                
            new_qualities += array.array('B', [last_qual]*insertions[i])
        
        tmp_relative_start += target_ranges[i][2]
        
        if (deletions[i] > 0):
            tmp_relative_start += deletions[i]
            out['has_deletion'] = True
    
    # finally add the last interval
    tmp_relative_end    = tmp_relative_start+target_ranges[-1][2]
    #add_tmp    = ''
    target_tmp = target_fasta.fetch(read_chr, target_ranges[-1][0], target_ranges[-1][1]).upper()
    query_tmp  = query_fasta.fetch(query_chr, query_ranges[-1][0], query_ranges[-1][1]).upper()
    read_add   = read_seq[tmp_relative_start:tmp_relative_end] 

    add_tmp  = list(query_tmp)
    for a,b,j in zip(target_tmp, 
                     read_add,
                     range(target_ranges[-1][2])):
      if a != b:
        add_tmp[j]  = read_add[j]
    add_tmp = ''.join(add_tmp)
        
    #for j in range(target_ranges[-1][2]):
    #    if (target_tmp[j] == read_add[j]):
    #        #target_add += target_tmp[j]
    #        add_tmp  += query_tmp[j]
    #    else:
    #        #target_add += read_add[j]
    #        add_tmp  += read_add[j]
        
    new_qualities += read_quality[tmp_relative_start:tmp_relative_end] 

    if is_reverse:
        query_add += revcomp_DNA(add_tmp)
    else:
        query_add += add_tmp
    
    out['query_sequence'] += query_add
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

def liftover_segment(chain, old_alignment, target_fasta, query_fasta, read_chr):
    # we will build a list of possible matches with every chain that covers this read
    # we will then return the one in the end that has the best score
    
    read_start   = old_alignment.reference_start
    read_seq     = old_alignment.query_sequence
    read_quality = old_alignment.query_qualities
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
        match(tup[0]):
            case (0):
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
                    else: # we have a perfect 1:1 matching
                        out, this_absolute_start, this_relative_start = add_solid_interval(out,
                                read_chr, intervals, this_absolute_start, this_absolute_end, 
                                this_relative_start, this_relative_end, this_add, target_fasta, 
                                query_fasta, read_seq, tup, read_quality)
                                                                                                  
                else: # we have gaps in the alignment
                    out['has_indel'] = True
                    if (intervals[0].start > this_absolute_start or intervals[-1].end < this_absolute_end):
                        out['pass'] = False
                        out['error type'] = 3
                        continue
                    else:
                        out, this_absolute_start, this_relative_start = add_gapped_interval(out,
                            read_chr, intervals, this_absolute_start, this_absolute_end, 
                            this_relative_start, this_relative_end, this_add, target_fasta, 
                            query_fasta, read_seq, tup, read_quality)
                        
 
            case (1 | 4): # insertion or soft clip. We add these sequences and add to the relative read position but not genomic
                if not out['pass']: continue
                this_add = tup[1]
                this_relative_end = this_relative_start + this_add
                seq_add  = read_seq[this_relative_start:this_relative_end]
                qual_add = read_quality[this_relative_start:this_relative_end]
                out['query_sequence']  += seq_add
                out['cigartuples'].append(tup)
                out['qualityscores'].extend(qual_add)
                out['segments'].append(None)
                this_relative_start = this_relative_end
            case (2 | 3): # deletion or skip. Add to genomic position but add no sequence
                if not out['pass']: continue
                this_add = tup[1]
                this_absolute_start += this_add   
                out['cigartuples'].append(tup)
                out['segments'].append(None)
            
    
    # update the length of Ns
    if out['pass']:
        if (len(out['query_sequence']) > 500):
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

# There doesnt seem to be a good way of dividing the genome into reasonable chunks. I want a function that 
# will divide the genome in to N equal bp chunks

# divide the genome up into equal sized chunks
def get_genome_chunks(tSizes,n):
    out = []
    out.append([])
    
    # calculate the genome size, and then chunk size
    gsize = 0
    for key, value in tSizes.items():
        gsize += value
        
    chunk_size = math.ceil(gsize/n)
    current_chunk_size = chunk_size
    i = 0
    
    for chrom, size in tSizes.items():
        cur_pos = 0
        while (1):
            next_pos = cur_pos + current_chunk_size
            if next_pos > size:
                next_pos = size
                out[i].append((chrom, cur_pos, next_pos))
                current_chunk_size -= (next_pos - cur_pos) 
                break
            else:
                out[i].append((chrom, cur_pos, next_pos))
                current_chunk_size = chunk_size
                cur_pos = next_pos
                i += 1
                out.append([])       
    return out
    
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
               best           = None):
    
    OUT_FILE_QUERY = outfile
    
    new_header = pysam.AlignmentHeader.from_dict(new_header)
    
    index_stats = SAMFILE.get_index_statistics()
    target_contig_list = []
    for i in index_stats:
      if i[3] != 0: target_contig_list.append(i[0])
      
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
          new_read = liftover_segment(chains[i], old_alignment, target_fasta, query_fasta, read_chr)
          
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


        # build the new alignment
        new_alignment = pysam.AlignedRead(new_header) # create AlignedRead object
        
        new_alignment.query_name = old_alignment.query_name                
        new_alignment.set_tags(old_alignment.get_tags() )                

        new_alignment.next_reference_id = -1
        new_alignment.next_reference_start = 0
        new_alignment.template_length = 0
        
        new_alignment.reference_id    = name_to_id[new_read['query_chrom']]
        new_alignment.reference_start = new_read['query_pos']
        new_alignment.mapping_quality = old_alignment.mapping_quality
                             
        # the is_reverse flag tells us if the chain is on the reverse strand, NOT if the read sequence changed
        # if 
        # we need to add the bam flag 
        new_alignment.flag = 0x0
        
        #if (new_read['is_reverse'] and old_alignment.is_forward): print("here")
        
        if old_alignment.is_reverse != new_read['is_reverse']:
            #print(old_alignment.is_reverse, new_read['is_reverse'], "reversing")
            new_alignment.flag = new_alignment.flag | 0x10
            
        if new_read['is_reverse']:
            new_alignment.query_sequence  = revcomp_DNA(new_read['query_sequence'])
            new_alignment.query_qualities = new_read['qualityscores'][::-1]
            new_alignment.cigartuples     = new_read['cigartuples'][::-1]
        else:
            new_alignment.query_sequence  = new_read['query_sequence']
            new_alignment.query_qualities = new_read['qualityscores']
            new_alignment.cigartuples     = new_read['cigartuples']

        # fix tags that get renamed
        try:
            rg, rgt = old_alignment.get_tag("RG", with_value_type=True)
        except KeyError:
            pass
        else:
            new_alignment.set_tag("RG", str(rg), rgt)
            
        OUT_FILE_QUERY.write(new_alignment)
    
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
               best           = None):

    new_header = pysam.AlignmentHeader.from_dict(new_header)
    OUT_FILE_QUERY = outfile
    
    index_stats = SAMFILE.get_index_statistics()
    target_contig_list = []
    for i in index_stats:
      if i[3] != 0: target_contig_list.append(i[0])
    
    nreads = n0 = n1 = n2 = n3 = n4 = 0
    for old1, old2 in read_pair_generator(SAMFILE):
        #print(old1)
        #print(old2)

        nreads += 2
        
        #if nreads > 100000: break
        #print(nreads, file=sys.stderr)
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
          new1 = liftover_segment(chains[i], old1, target_fasta, query_fasta, read1_chr)
          new2 = liftover_segment(chains[i], old2, target_fasta, query_fasta, read2_chr)
        
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
        
        # build the new alignment
        new_alignment1 = pysam.AlignedRead(new_header) # create AlignedRead object
        new_alignment2 = pysam.AlignedRead(new_header)
        
        new_alignment1.query_name = old1.query_name                
        new_alignment1.set_tags(old1.get_tags() )                

        new_alignment1.next_reference_id = -1
        new_alignment1.next_reference_start = 0
        new_alignment1.template_length = 0
        
        new_alignment2.query_name = old2.query_name                
        new_alignment2.set_tags(old2.get_tags() )                

        new_alignment2.next_reference_id = -1
        new_alignment2.next_reference_start = 0
        new_alignment2.template_length = 0
        
        new_alignment1.reference_id    = name_to_id[new1['query_chrom']]
        new_alignment1.reference_start = new1['query_pos']
        new_alignment1.mapping_quality = old1.mapping_quality
        
        new_alignment2.reference_id    = name_to_id[new2['query_chrom']]
        new_alignment2.reference_start = new2['query_pos']
        new_alignment2.mapping_quality = old2.mapping_quality
                             
        # the is_reverse flag tells us if the chain is on the reverse strand, NOT if the read sequence changed
        # we need to add the bam flag 
        new_alignment1.flag = 0x0
        new_alignment2.flag = 0x0
        
        if old1.is_reverse != new1['is_reverse']:
            new_alignment1.flag = new_alignment1.flag | 0x10
        
        if old2.is_reverse != new2['is_reverse']:
            new_alignment2.flag = new_alignment2.flag | 0x10
            
        if new1['is_reverse']:
            new_alignment1.query_sequence  = revcomp_DNA(new1['query_sequence'])
            new_alignment1.query_qualities = new1['qualityscores'][::-1]
            new_alignment1.cigartuples     = new1['cigartuples'][::-1]
        else:
            new_alignment1.query_sequence  = new1['query_sequence']
            new_alignment1.query_qualities = new1['qualityscores']
            new_alignment1.cigartuples     = new1['cigartuples']
            
        if new2['is_reverse']:
            new_alignment2.query_sequence  = revcomp_DNA(new2['query_sequence'])
            new_alignment2.query_qualities = new2['qualityscores'][::-1]
            new_alignment2.cigartuples     = new2['cigartuples'][::-1]
        else:
            new_alignment2.query_sequence  = new2['query_sequence']
            new_alignment2.query_qualities = new2['qualityscores']
            new_alignment2.cigartuples     = new2['cigartuples']

        # fix tags that get renamed
        try:
            rg, rgt = old1.get_tag("RG", with_value_type=True)
        except KeyError:
            pass
        else:
            new_alignment1.set_tag("RG", str(rg), rgt)
        
        try:
            rg, rgt = old2.get_tag("RG", with_value_type=True)
        except KeyError:
            pass
        else:
            new_alignment2.set_tag("RG", str(rg), rgt)
        
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
        
        # get the new insert length, rejecting reads longer than 10kb
        if (abs(new_tlen) >  10000): 
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

        #print(old_alignment)
        OUT_FILE_QUERY.write(new_alignment1)
        OUT_FILE_QUERY.write(new_alignment2)
        
    #OUT_FILE_TARGET.close()
    #OUT_FILE_QUERY.close()
    
    return (nreads, n0, n1, n2, n3, n4)
