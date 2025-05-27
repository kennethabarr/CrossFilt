#!/usr/bin/python3

import sys
import argparse
import pysam
import logging
import array
from timeit import default_timer as timer
import math
import os
from collections import defaultdict
import importlib.metadata

__version__ = importlib.metadata.version('crossfilt')

def main():
  parser = argparse.ArgumentParser(
                      prog='crossfilt-filter',
                      description='Outputs reads from bam1 that that have identical contig, position, CIGAR string, and XF tag (optional) in bam2')
  
  parser.add_argument("-x", "--xf", required=False, help="Require identical XF tag", action="store_true")
  parser.add_argument("bam1", help="Input bam files.")
  parser.add_argument("bam2", help="Input bam files.")
  parser.add_argument('--version', action='version',
                    version='CrossFilt v{version}'.format(version=__version__))
  
  args = parser.parse_args()
  use_xf = args.xf
  
  name_sorted = False
  
  if not os.path.exists(args.bam1 + ".bai"):
    print("Warning: " + args.bam1 + ".bai not found. Assuming files are filtered and sorted by read name.", file=sys.stderr)
    name_sorted = True
    
  if not os.path.exists(args.bam2 + ".bai"):
    print("Warning: " + args.bam2 + ".bai not found. Assuming files are filtered and sorted by read name.", file=sys.stderr)
    name_sorted = True
    
  
  # This function will pull read pairs from two position sorted files. It will 
  # cache reads from the second argument 
  
  def read_pair_generator(bam1, bam2, chrom):
    
    read2_iter = bam2.fetch(chrom)
    read2_dict = defaultdict()
    last_read2_pos = 0
    
    for read1 in bam1.fetch(chrom):
      
      qname = read1.query_name
      pos   = read1.reference_start
      r1    = read1.is_read1
      
      read1_id = qname + str(r1)
      
      # add reads from file2 till the position in that file is greater than this read
      while last_read2_pos <= pos:
        try:
          read2 = next(read2_iter)
          last_read2_pos = read2.reference_start
          read2_dict[read2.query_name + str(read2.is_read1)] = read2
        except StopIteration:
          break
      
      if read1_id in read2_dict:
        yield read1, read2_dict[read1_id]
        del read2_dict[read1_id]
        
      # clear out the cache of reads that we have passed
      for read2_id, read2 in read2_dict.copy().items():
        if read2.reference_start < pos:
          del read2_dict[read2_id]
    
  
  def get_read_count(file):
    contig_list = []
    total_reads = 0
    index_stats = file.get_index_statistics()
    
    for i in index_stats:
      if i[3] != 0:
        contig_list.append(i[0])
        total_reads += i[3]
  
    return total_reads, contig_list
  
  
  SAMFILE1 = pysam.AlignmentFile(args.bam1, "rb")
  SAMFILE2 = pysam.AlignmentFile(args.bam2, "rb")
  OUTFILE  = pysam.AlignmentFile('-', "wb", template=SAMFILE1)
  
  if not name_sorted:
    file1_total_reads, file1_contigs = get_read_count(SAMFILE1)
    file2_total_reads, file2_contigs = get_read_count(SAMFILE2)
     
    # this will be more efficient if bam2 is the smaller file
    i = matched = 0
    for contig in file1_contigs:
      for read1, read2 in read_pair_generator(SAMFILE1, SAMFILE2, contig):
        i += 1
        if not read1.reference_start == read2.reference_start: continue
        if not read1.reference_name  == read2.reference_name: continue
        if not read1.cigarstring     == read2.cigarstring: continue
        
        if use_xf:
          if not read1.has_tag("XF"): continue
          if not read2.has_tag("XF"): continue
          if not read1.get_tag("XF") ==  read2.get_tag("XF"): continue
    
        matched += 1
        OUTFILE.write(read1)
    
    print(str(matched) + ' (' + str(round(100*matched/i,2)) + '%) successfully matched', file=sys.stderr)
  
  else:
    iter1 = SAMFILE1.fetch(until_eof = True)
    iter2 = SAMFILE2.fetch(until_eof = True)
      
    i = matched = 0
    for read1, read2 in zip(iter1, iter2):
      i += 1
      # check read names to make sure they match
      if not read1.query_name == read2.query_name:
        sys.exit("Error: Read number " + str(i) + " query names are not identical (" + read1.query_name + " and " + read2.query_name + ")\nUse position sorted files or filter and sort your bam files by name.")
      
      if use_xf:
        if not read1.has_tag("XF"): continue
        if not read2.has_tag("XF"): continue
        if not read1.get_tag("XF") ==  read2.get_tag("XF"): continue
    
      matched += 1
      OUTFILE.write(read1)
    
    print(matched,"(",round(100*matched/i,2),"%) successfully matched",file=sys.stderr)
  
if __name__ == '__main__':
    main()
