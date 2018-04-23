# Assembler Components

Components of genome sequence assembly tools

# File formats

- FASTQ reads
- FASTA contigs
- [GFA](http://gfa-spec.github.io/GFA-spec/) assembly graph
- [SAM/BAM](http://samtools.github.io/hts-specs/) or [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) alignments of reads to draft assembly

**FASTA/FASTQ**: Optional attributes are found in the comment field, and formatted as in SAM, `XX:x:xxxx`. The comment field follows a space in the header.

**FASTQ**: Paired-end reads are interleaved and may be compressed. Linked read barcodes may be indicated with the `BX` tag.

**FASTA**: The sequence should not be line wrapped. The FASTA file should be indexed (FAI).

**GFA**: GFA2 is preferred over GFA1. The sequence fields may be empty (`*`). The sequences may be stored in an adjacent FASTA file, named for example `assembly.gfa` and `assembly.fa`.

**SAM/BAM**: SAM/BAM files are sorted by position and index, unless otherwise stated. Different sequencing libraries may be indicated with the read group `RG` attribute.

## GFA record types

- GFA (S): sequence segments
- GFA (E): overlap edges
- GFA (G): gap edges
- GFA (U): unordered groups of sequence segments
- GFA (O): ordered paths of sequence segments

## GFA sequence segment attributes

- GFA (S[RC]): read counts
- GFA (S[CN]): copy number estimate

# Notation

type1(record[attributes],&hellip;) + &hellip; &rarr; type2(record[attributes],&hellip;) + &hellip;

This stage requires a file of type1 and produces a file of type2. For example, estimate copy number of unitigs. A GFA file of sequence segments and edges and a BAM or PAF file of mapped reads produces a GFA file with estimated copy numbers of unitigs.

GFA (SE) + BAM/PAF &rarr; GFA (S[CN],E)

# Stages of genome assembly

## Preprocess reads

Remove sequencing artifacts specific to each sequencing technology.

FASTQ &rarr; FASTQ

- Trim adapter sequences
- Split chimeric reads
- Merge overlapping paired-end reads
- Extract barcode sequences

## Correct reads

Correct sequencing errors in reads.

FASTQ &rarr; FASTQ

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

GFA (SE) + BAM/PAF &rarr; GFA (S[CN],E)

- Count reads per unitig: GFA (SE) + BAM/PAF &rarr; GFA (S[RC],E)
- Estimate copy number: GFA (S[RC],E) &rarr; GFA (S[CN],E)

## Link unitigs

Identify pairs of unitigs that are proximal. Estimate their relative order, orientation, and distance.

GFA (SE) + BAM/PAF &rarr; GFA (SEG)

## Order and orient

Order and orient paths of unitigs. Contigs are created by contiguous paths of sequence segments. Scaffolds are created by discontiguous paths of sequence segments.

GFA (SEG) &rarr; GFA (SEO)

## Contract paths

Glue vertices of paths and replace each path with a single sequence segment.

GFA (SEO) &rarr; GFA (SE)

## Scaffold

Scaffolding is the combination of the three stages of link unitigs, order and orient, and contract paths.

FASTA/GFA(S) + BAM/PAF &rarr; FASTA/GFA(S)

# Tools

- [ABySS](https://github.com/bcgsc/abyss#readme)
- [BCALM2](https://github.com/GATB/bcalm#readme)
- [BCOOL](https://github.com/Malfoy/BCOOL#readme)
- [BFC](https://github.com/lh3/bfc#readme)
- [lh3/gfa1](https://github.com/lh3/gfa1#readme)
- [SGA](https://github.com/jts/sga#readme)
- [Unicycler](https://github.com/rrwick/Unicycler#readme)

A tool may combine multiple assembly stages in a single tool.

## Preprocess reads

Tools are specific to each sequencing technology and numerous and so will not be listed here.

## Correct reads

- BFC `bfc`
- BCOOL `Bcool.py`
- SGA `sga index | sga correct`

## Unitig

- ABySS `ABYSS` or `ABYSS-P` or `abyss-bloom-dbg` then `AdjList` or `abyss-overlap`
- BCALM2  `bcalm | convertToGFA.py`
- SGA `sga index | sga filter | sga overlap | sga assemble`

## Denoise

- ABySS `abyss-filtergraph`
- lh3/gfa1 `gfaview -t`

## Identify/collapse heterozygous sequences

- ABySS `PopBubbles`
- lh3/gfa1 `gfaview -b`

## Thread reads

- Minimap2 `minimap2`
- Unicycler `unicycler_align`

## Map reads

- BWA `bwa mem`
- Minimap2 `minimap2 -xsr` for short reads

## Estimate copy number

- SGA `sga-astat.py`

## Link unitigs

- ABySS `DistanceEst` for paired-end and mate-pair reads
- ABySS `abyss-longseqdist` for long reads
- ARCS `arcs` for linked reads

## Order and orient

- ABySS `abyss-scaffold` or `SimpleGraph | MergePaths`
- SGA `sga scaffold`

## Contract paths

- ABySS `MergeContigs`
- SGA `sga scaffold2fasta`
