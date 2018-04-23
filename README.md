# Assembler Components

Components of genome sequence assembly tools

# File formats

- FASTQ reads
- FASTA contigs
- [GFA2](http://gfa-spec.github.io/GFA-spec/) assembly graph
- [SAM/BAM](http://samtools.github.io/hts-specs/) or [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) alignments of reads to draft assembly

# Stages of genome assembly

## Unitig

Assemble unitigs by de Bruijn graph assembly.

FASTQ &rarr; GFA (SE)

- Count *k*-mers
- Filter *k*-mers by abundance
- Compact the graph

## Denoise

Remove sequences and edges caused by sequencing error.

GFA (SE) &rarr; GFA (SE)

- Prune tips
- Collapse bulges due to sequencing errors

## Identify/collapse heterozygous sequences

Identify bulges caused by heterozygosity.

GFA (SE) &rarr; GFA (SE)

- Identify bulges: GFA (SE) &rarr; GFA (SEU)
- Collapse bulges: GFA (SEU) &rarr; GFA (SE)
