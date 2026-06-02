#!/bin/bash 

# set up the conda environment with:
# conda create -n crossfilt bioconda::crossfilt bioconda::star bioconda::samtools bioconda::htseq 'python>=3.10'
# conda activate crossfilt

# This will run a short test of the CrossFilt pipeline. This example problem
# has a small set of 50bp single-end reads. For efficiency the provided
# reads, chains, and genomes only include chr22. In this example I am filtering
# the reads for those that do not show alignment or annotation bias in a 
# comparisson between humans and chimpanzees.

set -euo pipefail
SECONDS=0

############ INPUT FILES #################

TARGET_FA=genomes/human.fa.gz
QUERY_FA=genomes/chimp.fa.gz

TARGET_GTF=gtfs/human.gtf.gz
QUERY_GTF=gtfs/chimp.gtf.gz

TARGET2QUERY=chains/human2chimp.chain.gz
QUERY2TARGET=chains/chimp2human.chain.gz

THREADS=24

############ BUILD STAR REFERENCES #################

# In this example I am aligning with STAR. First we need to build the
# references.

mkdir -p STARrefs/human STARrefs/chimp

#STAR cant use compressed files, so uncompress first
gunzip < "$TARGET_FA" > tmp.fa
gunzip < "$TARGET_GTF" > tmp.gtf

STAR \
  --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir STARrefs/human \
  --genomeFastaFiles tmp.fa \
  --sjdbGTFfile tmp.gtf \
  --sjdbOverhang 49 \
  --sjdbGTFtagExonParentGene gene_name \
  --genomeSAindexNbases 11 \
  --readFilesCommand zcat

gunzip < "$QUERY_FA" > tmp.fa
gunzip < "$QUERY_GTF" > tmp.gtf

STAR \
  --runThreadN "$THREADS" \
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

mkdir -p tmp_target_alignment

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
  tmp_target_alignment/testAligned.sortedByCoord.out.bam "$TARGET_GTF" > \
  tmp_target_alignment/test.counts.txt 
  
# index this file
samtools index tmp_target_alignment/test.counted.bam


############ LIFT TO QUERY GENOME AND ALIGN TO QUERY #################

mkdir -p tmp_query_alignment/

# run liftover to query genome
crossfilt-lift \
  -i tmp_target_alignment/test.counted.bam \
  -o tmp_query_alignment/test.lifted \
  -c "$TARGET2QUERY" \
  -t "$TARGET_FA" \
  -q "$QUERY_FA"
  
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
  tmp_query_alignment/queryAligned.sortedByCoord.out.bam "$QUERY_GTF" > \
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
  -c "$QUERY2TARGET" \
  -t "$QUERY_FA" \
  -q "$TARGET_FA"
  
# find reads that lifted back and mapped to the same gene.
crossfilt-filter --xf \
  tmp_target_alignment/test.liftedback.bam \
  tmp_target_alignment/test.counted.bam > \
  test.filtered.bam
  
# the final file here is the test.filtered.bam file. It is now ready for counting
# with htseq-count to get unbiassed counts.

rm -r tmp_query_alignment
rm -r tmp_target_alignment

############ VERIFY SINGLE-END OUTPUT #################

expected_se="85a8d3c6ed3cdd9be659c4aa6bfb6f4c"
actual_se=$(samtools view test.filtered.bam | md5sum | awk '{print $1}')
if [ "$actual_se" != "$expected_se" ]; then
  echo "FAIL: test.filtered.bam checksum mismatch" >&2
  echo "  expected: $expected_se" >&2
  echo "  actual:   $actual_se" >&2
  exit 1
fi
echo "PASS: test.filtered.bam" >&2


############ REPEAT PROCESS TO TEST PAIRED ENDS #################


############ RUN INITIAL ALIGNMENT TO TARGET #################

mkdir -p tmp_target_alignment

STAR \
  --genomeDir STARrefs/human \
  --readFilesIn test_1.fastq.gz test_2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMmultNmax 1 \
  --outFilterMultimapNmax 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix tmp_target_alignment/test

# add the xf tag with ht-seq count
htseq-count --quiet --stranded=no -a 0 -i gene_name \
  -o tmp_target_alignment/test.counted.bam -p bam \
  tmp_target_alignment/testAligned.sortedByCoord.out.bam "$TARGET_GTF" > \
  tmp_target_alignment/test.counts.txt 
  
# index this file
samtools sort tmp_target_alignment/test.counted.bam > \
  tmp_target_alignment/test.counted.sorted.bam
samtools index tmp_target_alignment/test.counted.sorted.bam

############ LIFT TO QUERY GENOME AND ALIGN TO QUERY #################

mkdir -p tmp_query_alignment/

# run liftover to query genome
crossfilt-lift \
  -i tmp_target_alignment/test.counted.sorted.bam \
  -o tmp_query_alignment/test.lifted \
  -c "$TARGET2QUERY" \
  -t "$TARGET_FA" \
  -q "$QUERY_FA" \
  --paired
  
# convert back to fastq
samtools fastq -1 tmp_query_alignment/query_1.fa \
               -2 tmp_query_alignment/query_2.fa \
               tmp_query_alignment/test.lifted.bam > /dev/null
gzip tmp_query_alignment/query_1.fa
gzip tmp_query_alignment/query_2.fa

# realign to query genome
STAR \
  --genomeDir STARrefs/chimp \
  --readFilesIn tmp_query_alignment/query_1.fa.gz tmp_query_alignment/query_2.fa.gz \
  --readFilesCommand zcat \
  --outSAMmultNmax 1 \
  --outFilterMultimapNmax 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix tmp_query_alignment/query 
  
# add the xf tag with ht-seq count
htseq-count --quiet --stranded=no -a 0 -i gene_name \
  -o tmp_query_alignment/test.counted.bam -p bam \
  tmp_query_alignment/queryAligned.sortedByCoord.out.bam "$QUERY_GTF" > \
  tmp_query_alignment/test.counts.txt 

# index this file
samtools sort tmp_query_alignment/test.counted.bam > \
  tmp_query_alignment/test.counted.sorted.bam
samtools index tmp_query_alignment/test.counted.sorted.bam

# filter for reads that align to the same spot in the query. Note that
# we provide the counted bam first, so xf tag will carry over to the output. This
# tag is not in the lifted bam file. We will need this xf tag when we lift back
# to the human (target) genome.
crossfilt-filter \
  tmp_query_alignment/test.counted.sorted.bam \
  tmp_query_alignment/test.lifted.bam > \
  tmp_query_alignment/test.passed.bam


############ LIFT BACK TO TARGET #################

crossfilt-lift \
  -i tmp_query_alignment/test.passed.bam \
  -o tmp_target_alignment/test.liftedback \
  -c "$QUERY2TARGET" \
  -t "$QUERY_FA" \
  -q "$TARGET_FA" \
  --paired
  
# find reads that lifted back and mapped to the same gene.
crossfilt-filter --xf \
  tmp_target_alignment/test.liftedback.bam \
  tmp_target_alignment/test.counted.sorted.bam > \
  test_paired.filtered.bam
  
# the final file here is the test_paired.filtered.bam file. It is now ready for counting
# with htseq-count to get unbiassed counts.

rm -r tmp_query_alignment
rm -r tmp_target_alignment

############ VERIFY PAIRED-END OUTPUT #################

expected_pe="491558323726da049d99b10cb82937ab"
actual_pe=$(samtools view test_paired.filtered.bam | md5sum | awk '{print $1}')
if [ "$actual_pe" != "$expected_pe" ]; then
  echo "FAIL: test_paired.filtered.bam checksum mismatch" >&2
  echo "  expected: $expected_pe" >&2
  echo "  actual:   $actual_pe" >&2
  exit 1
fi
echo "PASS: test_paired.filtered.bam" >&2

############ RECORD TIMING #################

elapsed=$SECONDS
commit=$(git -C "$(dirname "$0")/.." rev-parse --short HEAD 2>/dev/null || echo "unknown")
log="$(dirname "$0")/timing_log.txt"
printf "%-20s  commit=%-8s  time=%ds\n" "$(date '+%Y-%m-%d %H:%M')" "$commit" "$elapsed" >> "$log"
echo "Test completed in ${elapsed}s  (log: $log)" >&2
