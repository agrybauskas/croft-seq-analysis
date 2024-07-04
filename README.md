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
    --bam sample-1.bam \
    --depth-change 20.0 \
    --include chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --reference reference.fasta \
    --seq GGAATCCCTTCTGCAGCACC \
    > sample-1.bed
```

Combining replication cleavage data (with AND operator):

```bash
./scripts/combine-cleavage-patterns.jl \
    --bed sample-1.bed,sample-2.bed,sample-3.bed \
    --operator AND \
    > results-sample-combined.bed
```

Combining replication cleavage data (with OR operator):

```bash
./scripts/combine-cleavage-patterns.jl \
    --bed results-control-1.bed,results-control-2.bed,results-control-3.bed \
    --operator OR \
    > results-control-combined.bed
```

Excluding one cleavage data from another (with NOT operator):

```bash
./scripts/combine-cleavage-patterns.jl \
    --bed results-sample-combined.bed,results-control-combined.bed \
    --operator NOT \
    > results-sample-excluded-control.bed
```

Visualising cleavage sites (PAM sequence can be added by using `--next-seq NGG` option):

```
./scripts/visualize-cleavage-patterns.jl \
    --no-alt-sequence \
    --show-gaps \
    --bed results-sample-excluded-control.bed \
    --top 1000 \
    --include chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --reference reference.fasta \
    --seq GGAATCCCTTCTGCAGCACC \
    --output-format HTML \
    > results-sample-excluded-control.html
```

For the information of all available options for the scripts use `-h` or `--help` options.