# Croft-seq analysis scripts

## Prerequisites

* Julia (tested on v1.3.1);
* Julia packages:
  * ArgParse
  * BioAlignments
  * BioSequences
  * CSVFiles
  * DataFrames
  * Formatting
  * Gumbo
  * Printf
  * PyCall
  * Statistics

## Usage

Finding potential cleavage sites (single-read):

```bash
./scripts/find-cleavage-patterns.jl \
    --ignore-R2 \
    --bg-include \
    --cleavage-ratio 0.95 \
    --reading-frame 4 \
    --bam sample.bam \
    --depth-change 20.0 \
    --include chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --reference reference.fasta \
    --seq GGAATCCCTTCTGCAGCACC \
    > results.bed
```

Finding potential cleavage sites (pair-end):

```bash
./scripts/find-cleavage-patterns.jl \
    --bg-include \
    --cleavage-ratio 0.95 \
    --reading-frame 4 \
    --bam sample.bam \
    --depth-change 20.0 \
    --include chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --reference reference.fasta \
    --seq GGAATCCCTTCTGCAGCACC \
    > results.bed
```

Combining replication data:

```
```

Visualising cleavage sites:

```
```