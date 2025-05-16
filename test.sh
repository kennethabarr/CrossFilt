#!/bin/bash

python liftover_bam.py \
  -i test/human.test.bam \
  -o test/human.lifted \
  -c test/chains/human2chimp.chain.gz \
  -t test/genomes/human.fa.gz \
  -q test/genomes/chimp.fa.gz 
  
python liftover_bam.py \
  -i test/human.lifted.bam \
  -o test/human.liftedback \
  -c test/chains/chimp2human.chain.gz \
  -t test/genomes/chimp.fa.gz \
  -q test/genomes/human.fa.gz 

python identical_reads.py test/human.liftedback.bam test/human.test.bam > /dev/null

