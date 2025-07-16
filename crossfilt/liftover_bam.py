import sys
import argparse
import pysam
import logging
import array
from timeit import default_timer as timer
from .lib import liftover_functions as lift
import math
import os
from random import randint
import importlib.metadata

__version__ = importlib.metadata.version('crossfilt')

def main():
  parser = argparse.ArgumentParser(
                        prog='crossfilt-lift',
                        description='Converts genome coordinates and nucleotide sequence for othologous segments in a BAM file')
    
  parser.add_argument("-i", "--input",         required=True,  help="The input BAM file to convert")
  parser.add_argument("-o", "--output",        required=True,  help="Name prefix for the output file")
  parser.add_argument("-c", "--chain",         required=True,  help="The UCSC chain file")
  parser.add_argument("-t", "--target-fasta",  required=True,  help="The genomic sequence of the target (the species we are converting from)")
  parser.add_argument("-q", "--query-fasta",   required=True,  help="The genomic sequence of the query (the species we are converting to)")
  parser.add_argument("-p", "--paired",        required=False, help="Add this flag if the reads are paired", action="store_true")
  parser.add_argument("-b", "--best",          required=False, help="Only attempt to lift using the best chain", action="store_true")

  parser.add_argument('--version', action='version',
                    version='CrossFilt v{version}'.format(version=__version__))
    
  #######################################################################
  ####################        Parse Arguments        ####################
  #######################################################################
  
  args = parser.parse_args()
  
  infile         = args.input
  outfile_prefix = args.output
  chainfile      = args.chain
  target_fasta   = args.target_fasta
  query_fasta    = args.query_fasta
  is_paired      = args.paired   
  best           = args.best
  
  TARGETFILE  = pysam.Fastafile(target_fasta)
  QUERYFILE   = pysam.Fastafile(query_fasta)
  
  if not os.path.exists(infile + ".bai"):
    pysam.index(infile)
  
  SAMFILE     = pysam.AlignmentFile(infile, "rb")
  index_stats = SAMFILE.get_index_statistics()
  old_header  = SAMFILE.header.to_dict()
  
  nreads_dict = {}
  target_contig_list = []
  total_reads = 0
  for i in index_stats:
    if i[3] != 0:
      target_contig_list.append(i[0])
      nreads_dict[i[0]] = i[3]
      total_reads += i[3]
      if i[0] not in TARGETFILE.references:
        sys.exit("Contig " + i[0] + " not found in fasta file. Did you use the right contig names?")
  
  
  #######################################################################
  ####################        Read Chain File        ####################
  #######################################################################
  
  start = timer()
  print("Reading chain file",file=sys.stderr)
  
  maps, tSizes, qSizes = lift.read_chain_file(chainfile, target_contig_list, QUERYFILE.references)
  
  # recalculate total reads, but only for reads that are in the chain file
  total_reads = 0
  for i in index_stats:
    if i[3] != 0:
      if i[0] in tSizes.keys():
        total_reads += i[3]
      
  if total_reads == 0:
    print("There were zero reads remaining after filtering out contigs without chains\n", file=sys.stderr)
    os.remove(tempname)
    exit()
    
  end = timer()
  print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)
  
  
  #######################################################################
  ####################        Lift Bam               ####################
  #######################################################################
  
  start = timer()
  print("Lifting",file=sys.stderr)
  
  comments=['ORIGINAL_BAM_FILE=' + infile]
  (new_header, name_to_id) = lift.bam_header_generator(orig_header = old_header, 
                                                       chrom_size  = qSizes, 
                                                       prog_name   = "CrossFilt",
                                                       prog_ver    = 1.0, 
                                                       format_ver  = 1.0,
                                                       sort_type   = 'coordinate',
                                                       co          = comments)
  
  
  
  tempname = outfile_prefix + '.temp' + str(randint(0,99999)) + '.bam'
  tempfile = pysam.Samfile(tempname, "wb", header = new_header)
  
  kwds = {'SAMFILE'        : SAMFILE, 
          'outfile'        : tempfile, 
          'old_header'     : old_header, 
          'new_header'     : new_header, 
          'target_fasta'   : TARGETFILE, 
          'query_fasta'    : QUERYFILE, 
          'maps'           : maps, 
          'tSizes'         : tSizes,
          'qSizes'         : qSizes,
          'name_to_id'     : name_to_id,
          'best'           : best}
    
  if is_paired:
    result = lift.process_pe(**kwds)
  else:
    result = lift.process_se(**kwds)
          
  tempfile.close()
  if result[1] == 0:
    print("Zero reads successfully lifted.\n", file=sys.stderr)
    
    exit()
  
  end = timer()
  print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)
  
  print("Processed",result[0],"reads",file=sys.stderr)
  print(result[1],"(",round(100*result[1]/result[0],2),"%) successfully mapped",file=sys.stderr)
  print(result[2],"had no enveloping chain",file=sys.stderr)
  print(result[3],"had no match in the chain",file=sys.stderr)
  print(result[4],"had a start or end that is an insertion/deletion",file=sys.stderr)
  print(result[5],"did not properly pair\n",file=sys.stderr)
  
  
  SAMFILE.close()
  TARGETFILE.close()
  QUERYFILE.close()
  
  
  #######################################################################
  ####################        Sort Output            ####################
  #######################################################################
  
  start = timer()
  print("Sorting",file=sys.stderr)
  
  pysam.sort("-o", outfile_prefix + ".bam", tempname)
  pysam.index(outfile_prefix + ".bam")
  os.remove(tempname)
  
  end = timer()
  print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)



if __name__ == '__main__':
    main()
