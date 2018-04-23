# Assembler Components

Components of genome sequence assembly tools

# File formats

- FASTQ reads
- FASTA contigs
- [GFA2](http://gfa-spec.github.io/GFA-spec/) assembly graph
- [SAM/BAM](http://samtools.github.io/hts-specs/) or [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) alignments of reads to draft assembly

## GFA record types

- GFA (S): sequence segments
- GFA (E): overlap edges
- GFA (G): gap edges
- GFA (U): unordered groups of sequence segments
- GFA (O): ordered paths of sequence segments

## GFA sequence segment attributes

- GFA (S[RC]): read counts
- GFA (S[CN]): copy number estimate

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

Collapsed bulges may select a single sequence or may compute a consensus sequence, using IUPAC ambiguity codes.

## Thread reads

Align reads to the assembly graph.

FASTQ + GFA (SE) + FASTQ &rarr; FASTQ + GFA (SE) + PAF

## Map reads

Map reads to their single best position in the draft genome.

FASTQ + FASTA &rarr; FASTQ + FASTA + BAM

## Estimate copy number

Estimate the copy number of each unitig.

FASTA/GFA(S) + BAM/PAF &rarr; GFA (S[CN])

- Count reads per unitig: GFA (S) + BAM/PAF &rarr; GFA (S[RC])
- Estimate copy number: GFA (S[RC]) &rarr; GFA (S[CN])

## Link unitigs

Identify pairs of unitigs that are proximal. Estimate their relative order, orientation, and distance.

GFA (SE) + BAM/PAF &rarr; GFA (SEG)

## Order and orient

Order and orient paths of unitigs. Contigs are created by contiguous paths of sequence segments. Scaffolds are created by discontiguous paths of sequence segments.

GFA (SEG) &rarr; GFA (SEO)

## Contract paths

Glue vertices of paths and replace each path with a single sequence segment.

GFA (SEO) &rarr; GFA (SE)
