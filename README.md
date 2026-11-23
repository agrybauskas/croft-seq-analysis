# Croft-seq analysis scripts

## Prerequisites

* Julia (v1.8.2) packages:
  * ArgParse
  * BioAlignments
  * BioSequences
  * CodecZlib
  * CSVFiles
  * DataFrames
  * FASTX
  * Formatting
  * Gumbo
  * PyCall
  * XAM
* Python (v3.10.12) library:
  * Biopython

## Usage

Finding potential cleavage sites (single-read):

```bash
./scripts/find-cleavage-patterns.jl \
    --ignore-R2 \
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
    > sample-combined.bed
```

Combining replication cleavage data (with OR operator):

```bash
./scripts/combine-cleavage-patterns.jl \
    --bed control-1.bed,control-2.bed,control-3.bed \
    --operator OR \
    > control-combined.bed
```

Excluding one cleavage data from another (with NOT operator):

```bash
./scripts/combine-cleavage-patterns.jl \
    --bed sample-combined.bed,control-combined.bed \
    --operator NOT \
    > sample-combined-excluded-control-combined.bed
```

Visualising cleavage sites (PAM sequence can be added by using `--show-next-seq NGG` option):

```
./scripts/visualize-cleavage-patterns.jl \
    --no-alt-sequence \
    --show-gaps \
    --bed sample-combined-excluded-control-combined.bed \
    --top 1000 \
    --include chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --reference reference.fasta \
    --seq GGAATCCCTTCTGCAGCACC \
    --output-format HTML \
    > sample-combined-excluded-control-combined.html
```

For the information of all available options for the scripts use `-h` or `--help` options.

## Citing

If you use *Croft-seq* in your research, please cite:
* Paulius Toliusis, Algirdas Grybauskas, Tomas Sinkunas, Tautvydas Karvelis,
Giedrius Sasnauskas, Mindaugas Zaremba, Off-target detection of CRISPR-Cas9
nuclease in vitro with CROFT-Seq, Synthetic Biology, Volume 11, Issue 1, 2026,
ysag006. https://doi.org/10.1093/synbio/ysag006
