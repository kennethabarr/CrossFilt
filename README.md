# CrossFilt

CrossFilt is a tool developed to identify reads that cause alignment or annotation bias in Cross-species genomic comparisons. We have tested it on RNA-seq and ATAC-seq, but it should be widely applicable to other genomic technologies. This tool works by lifting bam alignments from one species to another. This tool converts any sequence that matches the genome to that of the other species. Then we realign these reads in the other species. Finally, we lift the realigned reads back to the original genome and check which reads return the original coordinates. We only consider these reciprocally mapping reads in genomic comparisons

## Installation

This tool includes four scripts used for implementing this method: liftover_bam.py, identical_reads.py, identical_reads_xf.py and split_bam.py. These are all python scripts that require pysam. For convenience, we have provided a conda environment file. To create and load this envioment run the folloing from the base conda environment:

```
conda env create -f env/pysam.yml
conda activate pysam
```

We have provided a test script to verify that the tool is working properly. This script requires samtools, which is provided by the pysam environment. To run the test, run 
```
bash test.sh
```
This script will lift a set ~500k reads to and then from the chimpanzee genome, then check if they return the same original coordinates. On our system it takes about 2 minutes to run on a single thread. If things work properly it should finish with the message "460826 (100.0%) successfully matched"

## Tools

### liftover_bam.py

```
usage: liftover_bam.py [-h] -i INPUT -o OUTPUT -c CHAIN -t TARGET_FASTA -q QUERY_FASTA [-p] [-b]

Converts genome coordinates and nucleotide sequence for othologous segments in a BAM file

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input BAM file to convert
  -o OUTPUT, --output OUTPUT
                        Prefix for the output files
  -c CHAIN, --chain CHAIN
                        The UCSC chain file
  -t TARGET_FASTA, --target-fasta TARGET_FASTA
                        The genomic sequence of the target (the species we are converting from)
  -q QUERY_FASTA, --query-fasta QUERY_FASTA
                        The genomic sequence of the query (the species we are converting to)
  -p, --paired          Add this flag if the reads are paired
  -b, --best            Only attempt to lift using the best chain
```

This tool will lift reads from the target genome to the query genome using the provided chain file and genomes. It must be run on sorted and indexed bam files, so if the file is not sorted please do so using `samtools sort` and `samtools index`. It is compatible with single and paired end reads, which can be specified by the `--paired` flag. The output is written to a file specified by the output prefix flag. This fill will have '.query.sorted.bam' appended to the prefix to speficy that these reads are sorted and in query genome coordinates. For simple RNA-seq experiments these reads can then be converted back to fastq for realignment using `samtools fastq`. We have also used this on 10x genomics single-cell data using the 'bamtofastq' provided by 10x genomics. 

By default, if a read fails to lift on the best chain, this tool will proceed to the next best chain and try again. It will continue trying for all chains. A user can override this behavior with the `--best` flag, in which case the tool will only attempt to lift using the best chain. In our experience with primates this decreases the number of reads that successfully lift by about 5%, while decreasing the time it takes to run the tool by about 10%. 

In our hands, this tool takes about 2-3 minutes per 1M reads and for most human chain files it requires about 3GB of RAM. For large experiments this may be computationally expensive and we reccomend splitting the bam into smaller peices. The program will only store chains for chromosomes present in the bam file, so the memory requirements will decrease significantly when the bam file is split. For single-end reads you may split the bam file any way you like, but for paired-end reads it is essential that both ends are present in the same file. For that reason we have provided a tool split_bam.py that will split a file into equal sized peices. 

### liftover_bam.py 

```
usage: split_bam.py [-h] -i INPUT -o OUTPUT -c CHUNK_SIZE [-n NCPU] [-p]

Splits a bam file into equal sized chunks, keeping paired reads together

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input BAM file to split
  -o OUTPUT, --output OUTPUT
                        Prefix for the output files
  -c CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        The number of reads per file
  -n NCPU, --ncpu NCPU  The number of CPU cores to use
  -p, --paired          Add this flag if the reads are paired
```

To decrease run-time we reccomend splitting input bam files into smaller peices. This tool will split a bam file into chunks of CHUNK_SIZE reads. If reads are paired, it will ensure that both ends are kept in the same file. the number of CPUs passed to pysam for I/O and sorting can be changed with NCPU. 

### identical_reads.py and identical_reads_xf.py

```
Usage: identical_reads.py bam1 bam2 > reads.txt
Usage: identical_reads_xf.py bam1 bam2 > reads.txt
```

These tools will check whether the reads in two files are identical according to their chromosome, start position, and CIGAR string. Additionally, if running `identical_reads_xf.py`, it will check if the XF tag is identical in two files. The XF tag in a bam file is used by tools like STAR, htseq-count, and 10x cellranger to assign the feature that a read counts towards. 

To run these tools the files must be filtered to include the exact same seq of read names, and the files must be sorted by name. In our case, bam2 is created from lifting or realigning bam1, so bam2 necessarily contains a subset of bam1. To ensure an identical set of reads we only have to filter bam1 for reads in bam2. For example we can prepare these files using the following commands in samtools:

```
samtools view bam2 | cut -f1 > keep.txt
samtools view -b -N keep.txt bam1 > bam1.filtered
samtools sort -n bam1.filtered > bam1.filtered.byname
samtools sort -n bam2 > bam2.byname
identical_reads.py bam1.filtered.byname bam2.byname > reads.txt
```

The file reads.txt then contains read names that are identical in bam1 and bam2. 





