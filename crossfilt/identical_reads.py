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
from line_profiler import profile

__version__ = importlib.metadata.version('crossfilt')


# iterate through a bam file base-by-base and return a 
# dictionary of reads at that position.
def dict_generator(bam, chrom): 

  read_dict  = defaultdict()
  last_pos = 0
  
  for read in bam.fetch(chrom):
    this_pos = read.reference_start
    if this_pos > last_pos:
      yield last_pos, read_dict
      read_dict = defaultdict()
      
    read_dict[read.query_name + str(read.is_read1)] = read
    last_pos = this_pos
  yield last_pos, read_dict

  
# Read through two position-sorted bam files and return reads that have the
# same chromosome, position, and pair mate (read1/2)
def read_pair_generator(bam1, bam2, chrom):
    
  bam2_generator = dict_generator(bam2, chrom)
  
  # dict_pos holds the position in the dictionary
  dict_pos, read2_dict = next(bam2_generator)
  
  for read1 in bam1.fetch(chrom):
    
    pos = read1.reference_start

    # update the dictionary till the position is greater than or equal to the current position
    while dict_pos < pos:
      try:
        dict_pos, read2_dict = next(bam2_generator)
      except StopIteration:
        break

    # if the dictionary is greater we can keep moving with bam1
    if dict_pos != pos: continue 
    
    qname    = read1.query_name
    r1       = read1.is_read1
    read1_id = qname + str(r1)
    
    if read1_id in read2_dict:
     yield read1, read2_dict[read1_id]
        
def get_read_count(file):
  contig_list = []
  total_reads = 0
  index_stats = file.get_index_statistics()
  
  for i in index_stats:
    if i[3] != 0:
      contig_list.append(i[0])
      total_reads += i[3]

  return total_reads, contig_list
  
@profile
def tags_equal(read1, read2, tags):
  for tag in tags:
    if not read1.has_tag(tag): return False
    if not read2.has_tag(tag): return False
    if not read1.get_tag(tag) ==  read2.get_tag(tag): return False
  return True
  
@profile
def main():
  parser = argparse.ArgumentParser(
                      prog='crossfilt-filter',
                      description='Outputs reads from bam1 that that have identical contig, position, CIGAR string, and tag values (optional) in bam2')
  
  parser.add_argument("-t", "--tag",     required=False, help="Tag values to compare. Can be specified multiple times to compare multiple tags.", default=[], action="append")
  parser.add_argument("-x", "--xf",      required=False, help="Compare the XF tag. Equivalent to --tag XF", action="store_true")
  parser.add_argument("-@", "--threads", required=False, help="Number of compression/decrpression threads when reading/writing bam files.", default=1, type=int)
  parser.add_argument("bam1", help="Input bam file 1.")
  parser.add_argument("bam2", help="Input bam file 2.")
  parser.add_argument('--version', action='version',
                    version='CrossFilt v{version}'.format(version=__version__))
  
  args = parser.parse_args()
  tags = args.tag
  if args.xf: tags.append("XF")
  
  name_sorted = False
  
  if not os.path.exists(args.bam1 + ".bai"):
    print("Warning: " + args.bam1 + ".bai not found. Assuming files are filtered and sorted by read name.", file=sys.stderr)
    name_sorted = True
    
  if not os.path.exists(args.bam2 + ".bai"):
    print("Warning: " + args.bam2 + ".bai not found. Assuming files are filtered and sorted by read name.", file=sys.stderr)
    name_sorted = True
    
  
  # This function will pull read pairs from two position sorted files. It will 
  # cache reads from the second argument 
  
  SAMFILE1 = pysam.AlignmentFile(args.bam1, "rb", threads = args.threads)
  SAMFILE2 = pysam.AlignmentFile(args.bam2, "rb", threads = args.threads)
  OUTFILE  = pysam.AlignmentFile('-', "wb", template=SAMFILE1, threads = args.threads)
  
  if not name_sorted:
    start = timer()
    file1_total_reads, file1_contigs = get_read_count(SAMFILE1)
    file2_total_reads, file2_contigs = get_read_count(SAMFILE2)
     
    matched = 0
    for contig in file1_contigs:
      if contig in file2_contigs:
        for read1, read2 in read_pair_generator(SAMFILE1, SAMFILE2, contig):
          
          # we dont need to check if contig or position are the same because our 
          # generator only works if they are the same
          if not read1.cigarstring == read2.cigarstring: continue
          
          if not tags_equal(read1, read2, tags): continue
      
          matched += 1
          OUTFILE.write(read1)
    
    print(str(matched) + ' (' + str(round(100*matched/file1_total_reads,2)) + '%) successfully matched', file=sys.stderr)
    end = timer()
    print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)
  else:
    start = timer()
    iter1 = SAMFILE1.fetch(until_eof = True)
    iter2 = SAMFILE2.fetch(until_eof = True)
      
    i = matched = 0
    for read1, read2 in zip(iter1, iter2):
      i += 1
      # check read names to make sure they match
      if not read1.query_name == read2.query_name:
        sys.exit("Error: Read number " + str(i) + " query names are not identical (" + read1.query_name + " and " + read2.query_name + ")\nUse position sorted files or filter and sort your bam files by name.")
      
      if not tags_equal(read1, read2, tags): continue
    
      matched += 1
      OUTFILE.write(read1)
    
    print("processed ", i, " reads.",file=sys.stderr)
    print(matched,"(",round(100*matched/i,2),"%) successfully matched",file=sys.stderr)
    end = timer()
    print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)
  
if __name__ == '__main__':
    main()
