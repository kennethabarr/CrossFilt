# CrossFilt

CrossFilt is a tool developed to filter reads that cause alignment or annotation bias in cross-species genomic comparisons. We have tested it on RNA-seq and ATAC-seq, but it should be widely applicable to other genomic technologies. This tool works by lifting bam alignments from one species to another. This tool converts any sequence that matches the genome to that of the other species. Then we realign these reads in the other species. Finally, we lift the realigned reads back to the original genome and check which reads return the original coordinates. We only consider these reciprocally mapping reads in genomic comparisons.

## Changes

The change from v0.2.1 to v0.2.2 includes documentation improvements and minor fixes. Added NumPy-style docstrings throughout the codebase. Fixed several typos in help text and the README (orthologous, decompression, recommend, pieces). Standardized the thread flag to `-@/--threads` across all three commands (crossfilt-split previously used `-n/--ncpu`). Added a `.gitignore`, reference checksums for the test outputs, and timing logging to the test script.

The change from v0.1.5 to v0.2.0 includes major changes to the efficiency of crossfilt-filter. I have observed speedups of about 5x in some of our problems. I have also included two new flags in crossfilt-filter. The first is the --tag option. This will increase the flexibility of this tool, allowing users to choose different aligners that assign features in a different tag than htseq-count. I have also included an option to run on multiple threads, though note that these are compression/decompression threads and there will be minimal benefit beyond a few threads.

## Installation

Installation can be through pypi or conda/mamba (recommended). 

Install through [pypi](https://pypi.org/project/crossfilt/) with 

```
pip install crossfilt
```

or conda with

```
conda install bioconda::crossfilt
```

This will create four scripts for implementing our method: crossfilt-lift, crossfilt-filter, crossfilt-split, and crossfilt-slm. 

We have included a test script and input files to verify that your installation is working correctly. This also serves as an example of how to run this pipeline to get filtered, unbiased reads for cross-species comparisons. This test will require STAR, htseq-count, and samtools. To run the test, clone this repository, navigate to the test directory and run

```
conda create -n crossfilt bioconda::crossfilt bioconda::star bioconda::samtools bioconda::htseq
conda activate crossfilt
bash test.sh
```

This script will lift a set ~500k human chr22 reads to and then from the chimpanzee genome, then check if they return the same original coordinates and gene tag.

## Tools

### crossfilt-lift

```
usage: crossfilt-lift [-h] -i INPUT -o OUTPUT -c CHAIN -t TARGET_FASTA -q QUERY_FASTA [-p] [-b] [-@ THREADS] [--version]

Converts genome coordinates and nucleotide sequence for orthologous segments in a BAM file

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
  -@ THREADS, --threads THREADS
                        Number of compression/decompression threads when reading/writing bam files.
  --version             show program's version number and exit
```

This tool will lift reads from the target genome to the query genome using the provided chain file and genomes. It must be run on sorted and indexed bam files, so if the file is not sorted please do so using `samtools sort` and `samtools index`. It is compatible with single and paired end reads, which can be specified by the `--paired` flag. The output is written to a bam file specified by the output prefix flag. For simple RNA-seq experiments these reads can then be converted back to fastq for realignment using `samtools fastq`. We have also used this on 10x genomics single-cell data using the 'bamtofastq' script provided by 10x genomics. 

By default, if a read fails to lift on the best chain, this tool will proceed to the next best chain and try again. It will continue trying for all chains. A user can override this behavior with the `--best` flag, in which case the tool will only attempt to lift using the best chain. In our experience with primates this decreases the number of reads that successfully lift by about 5%, while decreasing the time it takes to run the tool by about 10%. 

In our hands, this tool takes about 2-3 minutes per 1M reads and for most human chain files it requires about 3GB of RAM. For large experiments this may be computationally expensive and we recommend splitting the bam into smaller pieces. The program will only store chains for chromosomes present in the bam file, so the memory requirements will decrease significantly when the bam file is split. For single-end reads you may split the bam file any way you like, but for paired-end reads it is essential that both ends are present in the same file. For that reason we have provided a tool `crossfilt-split` that will split a file into equal sized pieces. 

### crossfilt-split

```
usage: crossfilt-split [-h] -i INPUT -o OUTPUT [-@ THREADS] [-p] (-f NFILES | -s FILE_SIZE) [--version]

Splits a bam file into equal sized chunks, keeping paired reads together. This may return fewer files than expected if
many reads are missing a pair.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input BAM file to split
  -o OUTPUT, --output OUTPUT
                        Prefix for the output files
  -@ THREADS, --threads THREADS
                        Number of compression/decompression threads when reading/writing bam files.
  -p, --paired          Add this flag if the reads are paired
  -f NFILES, --nfiles NFILES
                        The number of files to split this into
  -s FILE_SIZE, --file-size FILE_SIZE
                        The number of reads per file
  --version             show program's version number and exit
```

To decrease run-time we recommend splitting input bam files into smaller pieces. The user can specify either the number of reads per file with FILE_SIZE or the number of files to split into with NFILES. If reads are paired, it will ensure that both ends are kept in the same file. The tool will compute the number of files needed based on the total reads present in the index, but if reads are paired and many reads don't have a mate present in the file then it is possible that it will produce fewer files than specified. The number of compression/decompression threads passed to pysam for I/O and sorting can be changed with `--threads`. 

### crossfilt-filter

```
usage: crossfilt-filter [-h] [-t TAG] [-x] [-@ THREADS] [--version] bam1 bam2

Outputs reads from bam1 that have identical contig, position, CIGAR string, and tag values (optional) in bam2

positional arguments:
  bam1                  Input bam file 1.
  bam2                  Input bam file 2.

options:
  -h, --help            show this help message and exit
  -t TAG, --tag TAG     Tag values to compare. Can be specified multiple times to compare multiple tags.
  -x, --xf              Compare the XF tag. Equivalent to --tag XF
  -@ THREADS, --threads THREADS
                        Number of compression/decompression threads when reading/writing bam files.
  --version             show program's version number and exit
```

This tool will check whether the reads in two files are identical according to their chromosome, start position, and CIGAR string. Additionally, it will check whether the values of tags are identical in two files. For instance, alignments with STAR or counts with htseq-count store the gene that the read counts towards in the XF tag. Cellranger stores this in the GN tag. The --tag option was added in v0.2.0 and the --xf tag was added to preserve reverse-compatibility. 

This tool will run on either position sorted and indexed files or on filtered and name sorted files. If an index file is not provided the tool will proceed under the assumption that reads appear in the exact same order in each file (i.e. both files contain the exact same set of reads and reads are sorted by read name).

This tool will output the bam1 reads that have perfect matches in bam2. 

If bam1 and bam2 are significantly different in size, this tool will be slightly more efficient if bam1 is the larger file.

### crossfilt-slm

```
usage: crossfilt-slm [-h] -i INPUT -o OUTPUT -c CHAIN -t TARGET_FASTA -q QUERY_FASTA -n NFILES [-p] [-b] [-@ THREADS] [--no-seq] [--keep-chunks] [--version]

Split a BAM by chromosome bins, lift each bin in parallel, and merge the results. Uses the BAM index for fast
chromosome extraction; each worker reads only the chain entries for its chromosomes.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input BAM file (coordinate-sorted; indexed or indexable)
  -o OUTPUT, --output OUTPUT
                        Output BAM prefix (result written to <prefix>.bam)
  -c CHAIN, --chain CHAIN
                        UCSC chain file for the liftover
  -t TARGET_FASTA, --target-fasta TARGET_FASTA
                        FASTA for the target genome (reads are currently aligned to this)
  -q QUERY_FASTA, --query-fasta QUERY_FASTA
                        FASTA for the query genome (reads will be lifted to this)
  -n NFILES, --nfiles NFILES
                        Number of parallel chromosome chunks
  -p, --paired          Input contains paired-end reads
  -b, --best            Only use the highest-scoring chain for each read
  -@ THREADS, --threads THREADS
                        Total threads available. During the parallel lift phase each worker uses 1 BAM I/O thread
                        (liftover is CPU-bound); sort threads per worker = threads // n_chunks. The final merge
                        uses all threads.
  --no-seq              Skip sequence conversion; only update coordinates and CIGAR
  --keep-chunks         Keep intermediate per-chunk BAM files after merging
  --version             show program's version number and exit
```

`crossfilt-slm` is a faster alternative to running `crossfilt-split` and `crossfilt-lift` separately. It combines splitting, lifting, and merging into a single command, and lifts all chromosome chunks in parallel.

Splitting is done by reading the BAM index to assign whole chromosomes to each chunk, rather than iterating over individual reads. Chromosomes are packed into N balanced bins using the LPT (Longest Processing Time First) greedy algorithm, which minimises the imbalance between bins without solving the full bin-packing problem. Each worker then lifts its own chromosomes and reads only the portion of the chain file it needs. Because paired reads are always on the same chromosome, pairs are guaranteed to stay together without any special handling.

Thread allocation is managed automatically from the single `-@` budget. During the parallel lift phase, each worker uses 1 BAM I/O thread (liftover is CPU-bound Python; extra I/O threads would compete across workers). The remaining threads are divided evenly across workers for the per-chunk sort. The final merge uses the full thread budget. For example, with `-n 4 -@ 20`, the lift phase runs 4 workers × 1 I/O thread, the sort phase runs 4 workers × 5 sort threads, and the merge uses 20 threads.

As a rough guide, set `-n` to the number of cores available and `-@` to the same value. For large experiments on a 24-core node, `-n 20 -@ 20` is a reasonable starting point.




