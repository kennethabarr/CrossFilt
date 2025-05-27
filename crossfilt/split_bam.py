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

def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        if not read.is_paired or read.is_secondary or read.is_supplementary:
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
       
       
def main():
  parser = argparse.ArgumentParser(
                      prog='crossfilt-split',
                      description='Splits a bam file into equal sized chunks, keeping paired reads together. This may return fewer files than expected if many reads are missing a pair.')
  
  parser.add_argument("-i", "--input",         required=True,  help="The input BAM file to split")
  parser.add_argument("-o", "--output",        required=True,  help="Prefix for the output files")
  parser.add_argument("-n", "--ncpu",          required=False, type=int, default=1, help="The number of CPU cores to use")
  parser.add_argument("-p", "--paired",        required=False, help="Add this flag if the reads are paired", action="store_true")
  
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument("-f", "--nfiles",    type=int, help="The number of files to split this into")
  group.add_argument("-s", "--file-size", type=int, help="The number of reads per file")
  
  parser.add_argument('--version', action='version',
                    version='CrossFilt v{version}'.format(version=__version__))
                    
  args = parser.parse_args()
  
  infile         = args.input
  outfile_prefix = args.output
  ncpu           = args.ncpu
  is_paired      = args.paired
  
  SAMFILE     = pysam.AlignmentFile(infile, "rb", threads=ncpu)
  old_header  = SAMFILE.header.to_dict()
  bamiter     = SAMFILE.fetch(until_eof=True)
  
  if args.file_size is None:
    total_reads = 0
    
    if not os.path.exists(infile + ".bai"): pysam.index(infile)
    
    index_stats = SAMFILE.get_index_statistics()
    for i in index_stats:
      total_reads += i[3]
          
    chunk_size = math.ceil(total_reads/args.nfiles) + 1
    
  else:
    chunk_size = args.file_size
  
  
  
  # set a list of contigs
  print("Splitting input file into files with " + str(chunk_size) + " reads...",file=sys.stderr)
  start = timer()
  
  chunk_iter = 0
  file_iter = 0
  print("Splitting chunk " + str(file_iter),file=sys.stderr)
  
  this_file = pysam.Samfile(outfile_prefix + "." + str(file_iter) + ".bam", "wb", header = old_header, threads=ncpu)
  
  if is_paired:
    for read1, read2 in read_pair_generator(SAMFILE):
      chunk_iter += 2
      
      this_file.write(read1)
      this_file.write(read2)
      if chunk_iter >= chunk_size:
        this_file.close()
        pysam.sort("-@", str(ncpu), "-o", outfile_prefix + ".sorted." + str(file_iter) + ".bam", outfile_prefix + "." + str(file_iter) + ".bam")
        pysam.index(outfile_prefix + ".sorted." + str(file_iter) + ".bam")
        os.remove(outfile_prefix + "." + str(file_iter) + ".bam")
        chunk_iter = 0
        file_iter += 1
        print("Splitting chunk " + str(file_iter),file=sys.stderr)
        this_file = pysam.Samfile(outfile_prefix + ".input." + str(file_iter) + ".bam", "wb", header = old_header, threads = max(round(ncpu/2), 1))
  else:
    for read in SAMFILE.fetch(until_eof=True):
      chunk_iter += 1
      this_file.write(read)
      if chunk_iter >= chunk_size:
        this_file.close()
        pysam.sort("-@", str(ncpu), "-o", outfile_prefix + ".sorted." + str(file_iter) + ".bam", outfile_prefix + "." + str(file_iter) + ".bam")
        pysam.index(outfile_prefix + ".sorted." + str(file_iter) + ".bam")
        os.remove(outfile_prefix + "." + str(file_iter) + ".bam")
        chunk_iter = 0
        file_iter += 1
        print("Splitting chunk " + str(file_iter),file=sys.stderr)
        this_file = pysam.Samfile(outfile_prefix + "." + str(file_iter) + ".bam", "wb", header = old_header, threads = ncpu)
              
          
  this_file.close()
  pysam.sort("-@", str(ncpu), "-o", outfile_prefix + ".sorted." + str(file_iter) + ".bam", outfile_prefix + "." + str(file_iter) + ".bam")
  pysam.index(outfile_prefix + ".sorted." + str(file_iter) + ".bam")
  os.remove(outfile_prefix + "." + str(file_iter) + ".bam")
  SAMFILE.close()
  
  end = timer()
  print("Complete in",round(end-start,2),"seconds\n",file=sys.stderr)
  
if __name__ == '__main__':
    main()
