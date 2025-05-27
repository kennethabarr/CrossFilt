#!/bin/bash 

# set up the conda environment with:
# conda create -n crossfilt bioconda::crossfilt bioconda::star bioconda::samtools bioconda::htseq 'python>=3.10'
# conda activate crossfilt

# This will run a short test of the CrossFilt pipeline. This example problem
# has a small set of 50bp single-end reads. For efficiency the provided
# reads, chains, and genomes only include chr22. In this example I am filtering
# the reads for those that do not show alignment or annotation bias in a 
# comparisson between humans and chimpanzees.


############ INPUT FILES #################

TARGET_FA=genomes/human.fa.gz
QUERY_FA=genomes/chimp.fa.gz

TARGET_GTF=gtfs/human.gtf.gz
QUERY_GTF=gtfs/chimp.gtf.gz

TARGET2QUERY=chains/human2chimp.chain.gz
QUERY2TARGET=chains/chimp2human.chain.gz

############ BUILD STAR REFERENCES #################

# In this example I am aligning with STAR. First we need to build the
# references.

mkdir STARrefs
mkdir STARrefs/human
mkdir STARrefs/chimp

#STAR cant use compressed files, so uncompress first
gunzip < $TARGET_FA > tmp.fa
gunzip < $TARGET_GTF > tmp.gtf

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir STARrefs/human \
  --genomeFastaFiles tmp.fa \
  --sjdbGTFfile tmp.gtf \
  --sjdbOverhang 49 \
  --sjdbGTFtagExonParentGene gene_name \
  --genomeSAindexNbases 11 \
  --readFilesCommand zcat

gunzip < $QUERY_FA > tmp.fa
gunzip < $QUERY_GTF > tmp.gtf

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir STARrefs/chimp \
  --genomeFastaFiles tmp.fa \
  --sjdbGTFfile tmp.gtf \
  --sjdbOverhang 49 \
  --sjdbGTFtagExonParentGene gene_name \
  --genomeSAindexNbases 11 \
  --readFilesCommand zcat

rm tmp.fa
rm tmp.gtf


############ RUN INITIAL ALIGNMENT TO TARGET #################

mkdir tmp_target_alignment

STAR \
  --genomeDir STARrefs/human \
  --readFilesIn test.fastq.gz \
  --readFilesCommand zcat \
  --outSAMmultNmax 1 \
  --outFilterMultimapNmax 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix tmp_target_alignment/test

# add the xf tag with ht-seq count
htseq-count --quiet --stranded=no -a 0 -i gene_name \
  -o tmp_target_alignment/test.counted.bam -p bam \
  tmp_target_alignment/testAligned.sortedByCoord.out.bam $TARGET_GTF > \
  tmp_target_alignment/test.counts.txt 
  
# index this file
samtools index tmp_target_alignment/test.counted.bam


############ LIFT TO QUERY GENOME AND ALIGN TO QUERY #################

mkdir tmp_query_alignment/

# run liftover to query genome
crossfilt-lift \
  -i tmp_target_alignment/test.counted.bam \
  -o tmp_query_alignment/test.lifted \
  -c $TARGET2QUERY \
  -t $TARGET_FA \
  -q $QUERY_FA
  
# convert back to fastq
samtools fastq tmp_query_alignment/test.lifted.bam | gzip > tmp_query_alignment/query.fa.gz

# realign to query genome
STAR \
  --genomeDir STARrefs/chimp \
  --readFilesIn tmp_query_alignment/query.fa.gz \
  --readFilesCommand zcat \
  --outSAMmultNmax 1 \
  --outFilterMultimapNmax 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix tmp_query_alignment/query 
  
# add the xf tag with ht-seq count
htseq-count --quiet --stranded=no -a 0 -i gene_name \
  -o tmp_query_alignment/test.counted.bam -p bam \
  tmp_query_alignment/queryAligned.sortedByCoord.out.bam $QUERY_GTF > \
  tmp_query_alignment/test.counts.txt 

# index this file
samtools index tmp_query_alignment/test.counted.bam

# filter for reads that align to the same spot in the query. Note that
# we provide the counted bam first, so xf tag will carry over to the output. This
# tag is not in the lifted bam file. We will need this xf tag when we lift back
# to the human (target) genome.
crossfilt-filter \
  tmp_query_alignment/test.counted.bam \
  tmp_query_alignment/test.lifted.bam > \
  tmp_query_alignment/test.passed.bam


############ LIFT BACK TO TARGET #################

crossfilt-lift \
  -i tmp_query_alignment/test.passed.bam \
  -o tmp_target_alignment/test.liftedback \
  -c $QUERY2TARGET \
  -t $QUERY_FA \
  -q $TARGET_FA
  
# find reads that lifted back and mapped to the same gene.
crossfilt-filter --xf \
  tmp_target_alignment/test.liftedback.bam \
  tmp_target_alignment/test.counted.bam > \
  test.filtered.bam
  
# the final file here is the test.filtered.bam file. It is now ready for counting
# with htseq-count to get unbiassed counts.

rm -r tmp_query_alignment
rm -r tmp_target_alignment


