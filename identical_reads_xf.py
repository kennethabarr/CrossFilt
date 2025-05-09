import sys
import argparse
import pysam
import logging
import array
import numpy as np
from timeit import default_timer as timer
from tqdm import tqdm
import math
import os

# total arguments
n = len(sys.argv)

if not n==3:
  raise Exception("Usage: identical_reads_xf.py bam1 bam2 > reads.txt")

SAMFILE1 = pysam.AlignmentFile(sys.argv[1], "rb")
SAMFILE2 = pysam.AlignmentFile(sys.argv[2], "rb")

iter1 = SAMFILE1.fetch(until_eof = True)
iter2 = SAMFILE2.fetch(until_eof = True)
    
i = matched = 0
for read1, read2 in zip(iter1, iter2):
  i += 1
  # check read names to make sure they match
  if not read1.query_name == read2.query_name:
    raise Exception("Read number" + str(i) + ": query names are not identical (" + read1.query_name + " and " + read2.query_name + ")\nFilter and sort your bam files by name.")
  
  if not read1.has_tag("XF"): continue
  if not read2.has_tag("XF"): continue
  if not read1.get_tag("XF") ==  read2.get_tag("XF"): continue
  
  matched += 1
  sys.stdout.write(read1.query_name + "\n")
  
print(matched,"(",round(100*matched/i,2),"%) successfully matched",file=sys.stderr)
