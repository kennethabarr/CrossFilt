# CrossFilt

CrossFilt is a tool developed to filter reads that cause alignment or annotation bias in cross-species genomic comparisons. We have tested it on RNA-seq and ATAC-seq, but it should be widely applicable to other genomic technologies. This tool works by lifting bam alignments from one species to another. This tool converts any sequence that matches the genome to that of the other species. Then we realign these reads in the other species. Finally, we lift the realigned reads back to the original genome and check which reads return the original coordinates. We only consider these reciprocally mapping reads in genomic comparisons.

## Installation

Installation can be through pypi or conda/mamba (reccomended). 

Install through pypi with 

```
pip install crossfilt
```

or conda with

```
conda install bioconda::crossfilt
```

This will create three scripts for implementing our method: crossfilt-lift, crossfilt-filter and crossfilt-split. 

We have included a test script and input files to verify that your installation is working corretly. This also serves as an example of how to run this pipeline to get filtered, unbiased reads for cross-species comparisons. This test will require STAR, htseq-count, and samtools. To run the test, clone this repository, navigate to the test directory and run

```
conda create -n crossfilt bioconda::crossfilt bioconda::star bioconda::samtools bioconda::htseq
conda activate crossfilt
bash test.sh
```

This script will lift a set ~500k human chr22 reads to and then from the chimpanzee genome, then check if they return the same original coordinates and gene tag.

## Tools

### crossfilt-lift

```
usage: crossfilt-lift [-h] -i INPUT -o OUTPUT -c CHAIN -t TARGET_FASTA -q QUERY_FASTA [-p] [-b] [--version]

Converts genome coordinates and nucleotide sequence for othologous segments in a BAM file

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input BAM file to convert
  -o OUTPUT, --output OUTPUT
                        Name prefix for the output file
  -c CHAIN, --chain CHAIN
                        The UCSC chain file
  -t TARGET_FASTA, --target-fasta TARGET_FASTA
                        The genomic sequence of the target (the species we are converting from)
  -q QUERY_FASTA, --query-fasta QUERY_FASTA
                        The genomic sequence of the query (the species we are converting to)
  -p, --paired          Add this flag if the reads are paired
  -b, --best            Only attempt to lift using the best chain
  --version             show program's version number and exit
```

This tool will lift reads from the target genome to the query genome using the provided chain file and genomes. It must be run on sorted and indexed bam files, so if the file is not sorted please do so using `samtools sort` and `samtools index`. It is compatible with single and paired end reads, which can be specified by the `--paired` flag. The output is written to a bam file specified by the output prefix flag. For simple RNA-seq experiments these reads can then be converted back to fastq for realignment using `samtools fastq`. We have also used this on 10x genomics single-cell data using the 'bamtofastq' script provided by 10x genomics. 

By default, if a read fails to lift on the best chain, this tool will proceed to the next best chain and try again. It will continue trying for all chains. A user can override this behavior with the `--best` flag, in which case the tool will only attempt to lift using the best chain. In our experience with primates this decreases the number of reads that successfully lift by about 5%, while decreasing the time it takes to run the tool by about 10%. 

In our hands, this tool takes about 2-3 minutes per 1M reads and for most human chain files it requires about 3GB of RAM. For large experiments this may be computationally expensive and we reccomend splitting the bam into smaller peices. The program will only store chains for chromosomes present in the bam file, so the memory requirements will decrease significantly when the bam file is split. For single-end reads you may split the bam file any way you like, but for paired-end reads it is essential that both ends are present in the same file. For that reason we have provided a tool split_bam.py that will split a file into equal sized peices. 

### crossfilt-split

```
usage: crossfilt-split [-h] -i INPUT -o OUTPUT [-n NCPU] [-p] (-f NFILES | -s FILE_SIZE) [--version]

Splits a bam file into equal sized chunks, keeping paired reads together. This may return fewer files than expected if
many reads are missing a pair.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input BAM file to split
  -o OUTPUT, --output OUTPUT
                        Prefix for the output files
  -n NCPU, --ncpu NCPU  The number of CPU cores to use
  -p, --paired          Add this flag if the reads are paired
  -f NFILES, --nfiles NFILES
                        The number of files to split this into
  -s FILE_SIZE, --file-size FILE_SIZE
                        The number of reads per file
  --version             show program's version number and exit
```

To decrease run-time we reccomend splitting input bam files into smaller peices. The user can specify either the number of reads per file with FILE_SIZE or the number of files to split into with NFILES. If reads are paired, it will ensure that both ends are kept in the same file. The tool will compute the number of files needed based on the total reads present in the index, but if reads are paired and many reads dont have a mate present in the file then it is possible that it will produce fewer files than specified. The number of CPUs passed to pysam for I/O and sorting can be changed with NCPU. 

### crossfilt-filter

```
usage: crossfilt-filter [-h] [-x] [--version] bam1 bam2

Outputs reads from bam1 that that have identical contig, position, CIGAR string, and XF tag (optional) in bam2

positional arguments:
  bam1        Input bam files.
  bam2        Input bam files.

options:
  -h, --help  show this help message and exit
  -x, --xf    Require identical XF tag
  --version   show program's version number and exit
```

Thise tool will check whether the reads in two files are identical according to their chromosome, start position, and CIGAR string. Additionally, if the optional xf flag is included, it will check if the XF tag is identical in two files. The XF tag in a bam file is used by tools like STAR, htseq-count, and 10x cellranger to assign the feature that a read counts towards. 

This tool will run on either position sorted and indexed files or on filtered and name sorted files. If an index file is not provided the tool will proceed under the assumption that reads appear in the exact same order in each file (i.e. both files contain the exact same set of reads and reads are sorted by read name).

This tool will output the bam1 reads that have perfect matches in bam2. 

If bam1 and bam2 are significantly different in size, this tool will be slightly more efficient if bam1 is the larger file.




