#!/bin/bash

python liftover_bam.py \
  -i test/human.test.bam \
  -o test/human.lifted \
  -c test/chains/human2chimp.chain.gz \
  -t test/genomes/human.fa.gz \
  -q test/genomes/chimp.fa.gz 
  
python liftover_bam.py \
  -i test/human.lifted.query.sorted.bam \
  -o test/human.liftedback \
  -c test/chains/chimp2human.chain.gz \
  -t test/genomes/chimp.fa.gz \
  -q test/genomes/human.fa.gz 

samtools sort -n test/human.liftedback.query.sorted.bam > test/human.liftedback.byname.bam
samtools sort -n test/human.test.bam > test/human.test.byname.bam

python identical_reads.py test/human.liftedback.byname.bam test/human.test.byname.bam > /dev/null

