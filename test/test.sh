#!/bin/bash

python liftover_bam.py \
  -i human.test.bam \
  -o human.lifted \
  -c chains/human2chimp.chain.gz \
  -t genomes/human.fa.gz \
  -q genomes/chimp.fa.gz 
  
python liftover_bam.py \
  -i human.lifted.bam \
  -o human.liftedback \
  -c chains/chimp2human.chain.gz \
  -t genomes/chimp.fa.gz \
  -q genomes/human.fa.gz 

python identical_reads.py human.liftedback.bam human.test.bam > /dev/null

