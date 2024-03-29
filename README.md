# Assembler Components

Components of genome sequence assembly tools

# Rationale

Genome, metagenome and transcriptome assemblers range from fully integrated to fully modular. Fully modular assembly has a number of benefits. This repository is ongoing work to define some important checkpoints in a modular assembly pipeline, along with standard input/output formats. For now we have a bias towards Illumina-type sequencing data (single reads, paired reads, mate-pairs, 10x), but we aim to make the components also compatible with 3rd generation reads.

Feel free to contribute via pull-requests.

# File formats

- FASTQ reads
- FASTA contigs
- [GFA](http://gfa-spec.github.io/GFA-spec/) assembly graph
- [SAM/BAM](http://samtools.github.io/hts-specs/) or [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) alignments of reads to draft assembly
- TSV (tab-separated values) tabular data

**FASTA/FASTQ**: Optional attributes are found in the comment field, and formatted as in SAM, `XX:x:xxxx`. The comment field follows a space in the header.

**FASTQ**: Paired-end reads are interleaved and may be compressed. Linked read barcodes may be indicated with the `BX` tag.

**FASTA**: The sequence should not be line wrapped. The FASTA file should be indexed (FAI).

**GFA**: GFA2 is preferred over GFA1. The sequence fields may be empty (`*`). The sequences may be stored in an adjacent FASTA file, named for example `assembly.gfa` and `assembly.fa`.

**SAM/BAM**: SAM/BAM files are sorted by position and indexed, unless otherwise stated. Different sequencing libraries may be indicated with the read group `RG` attribute.

## GFA record types

- GFA (S): sequence segments
- GFA (E): overlap edges
- GFA (G): gap edges
- GFA (U): unordered groups of sequence segments
- GFA (O): ordered paths of sequence segments

## GFA sequence segment attributes

- GFA (S[RC]): read counts
- GFA (S[DP]): depth of coverage
- GFA (S[CN]): copy number estimate

## File names

An assembly is contained in a single directory. The files are named according to the pattern `[0-9]+_[a-z]+\.[a-z.]+`. The numeric prefixes are zero-padded and identical in length, and they indicate the stage of the assembly. A descriptive name and file type extension follow. The files may be compressed.

### Example

```
0_pe.fq.gz
1_unitig.gfa
2_denoise.gfa
3_debulge.gfa 3_debulge.bam 3_debulge.bam.bai
4_link.gfa
5_scaffold.gfa
6_assembly.gfa 6_assembly.fa 6_assembly.fa.fai 6_assembly.bam 6_assembly.bam.bai
```

# Notation

type1(record[attributes],&hellip;) + &hellip; &rarr; type2(record[attributes],&hellip;) + &hellip;

This stage requires a file of type1 and produces a file of type2. For example, estimate copy number of unitigs. A GFA file of sequence segments and edges and a BAM or PAF file of mapped reads produces a GFA file with estimated copy numbers of unitigs.

GFA (SE) + BAM/PAF &rarr; GFA (S[CN],E)

# Stages of genome assembly

## Preprocess reads

Remove sequencing artifacts specific to each sequencing technology. Improve the quality of the input reads with minimal loss of information, for example heterozygous variants.

FASTQ &rarr; FASTQ

- Trim adapter sequences
- Extract barcode sequences
- Merge overlapping paired-end reads
- Split chimeric reads

## Quality control

Assess the quality of the sequencing, and estimate parameters of the genome.

FASTQ &rarr; TSV

- Assess the quality of the reads
- Estimate sequencing depth
- Estimate parameters of the genome, size as size, heterozygosity, and repetitiveness
- Predict assembly parameters, such as *k*-mer size and minimum *k*-mer abundance

## Correct reads

Correct sequencing errors in reads.

FASTQ &rarr; FASTQ

## Unitig

Assemble unitigs by de Bruijn graph (dBG) assembly or overlap, layout, consensus (OLC) assembly.

FASTQ &rarr; GFA (SE)

### de Bruijn Graph (dBG)

- Count *k*-mers
- Filter *k*-mers by abundance
- Compact the graph
- Compute the sequences of the unitigs

### Overlap, layout, consensus (OLC)

- Find all pairwise overlaps of reads
- Determine the order and orientation of reads
- Compute the consensus sequences of the unitigs

## Denoise

Remove sequencing errors from the assembly graph. Retain variants.

GFA (SE) &rarr; GFA (SE)

- Prune tips
- Remove bulges due to sequencing errors

## Collapse variants

Identify and/or remove variants from the graph.

GFA (SE) &rarr; GFA (SE)

- | **Identify variants**
  | Identify bulges and create unordered groups of sequence segments.
  | GFA (SE) &rarr; GFA (SEU)
- | **Collapse variants**
  | Collapse bulges, possibly creating new sequence segments.
  | GFA (SEU) &rarr; GFA (SE)

A single sequence or path through the bulge may be selected, or the bulge may be replaced by a consensus sequence, possibly using IUPAC ambiguity codes to represent the consensus.

## Map reads

Map reads to the assembly sequences.

FASTQ + GFA (S)/FASTA &rarr; BAM

## Estimate copy number

Estimate the copy number of each sequence segment.

GFA (SE) + BAM/PAF &rarr; GFA (S[CN],E)

- | **Calculate depth of each sequence segment**
  | Count mapped reads and calculate depth of coverage of each sequence segment.
  | GFA (SE) + BAM/PAF &rarr; GFA (S[RC,DP],E)
- | **Estimate the copy number of each sequence segment**
  | GFA (S[RC,DP],E) &rarr; GFA (S[RC,DP,CN],E)

Note that the median depth of coverage is more robust than the mean depth of coverage to the alignment artifacts caused by collapsed repeats, misaligned reads, and other issues.

## Resolve repeats and scaffold

Expand repeats, and order and orient sequence segments into contigs and scaffolds.

FASTA/GFA (SE) + BAM/PAF &rarr; FASTA/GFA (SE)

Contigs are contiguous sequences with no gaps. Creating contigs requires expanding the repetitive sequence found between the unique contigs. Contigs are derived from contiguous paths of sequence segments without any gaps. Scaffolds are derived from discontiguous paths of sequence segments with gaps between the segments.

A tool may implement scaffolding as a single stage of assembly. Scaffolding however may be viewed as composed of the three distinct stages: link unitigs, order and orient unitigs to construct paths, and contract paths to create new sequence segments.

- | **Link unitigs**
  | Identify pairs of unitigs that are proximal. Estimate their relative order, orientation, and distance.
  | GFA (SE) + BAM/PAF &rarr; GFA (SEG)
- | **Order and orient**
  | Order and orient paths of unitigs.
  | GFA (SEG) &rarr; GFA (SEO)
- | **Contract paths**
  | Glue vertices of paths and replace each path with a single sequence segment.
  | GFA (SEO) &rarr; GFA (SE)

## Fill gaps

Assemble the sequence found in the scaffold gaps between adjacent contigs.

FASTQ + FASTA/GFA (S) &rarr; FASTA/GFA (S)

## Polish

Map the reads to the assembly and correct assembly errors.

FASTQ + FASTA/GFA (SE) + BAM/PAF &rarr; FASTA/GFA (SE)

## Visualize the assembly graph

GFA (SE) &rarr; PNG/SVG

## Assess assembly quality

Assess the contiguity and correctness of the assembly.

FASTA/GFA (S) &rarr; TSV

- Compute assembly metrics, such as N50 and NG50
- Align the assembly to the reference genome
- Compute assembly metrics, such as NGA50 and number of misassemblies

# Tools

- [ABySS](https://github.com/bcgsc/abyss#readme)
- [Bandage](https://github.com/rrwick/Bandage#readme)
- [BCALM2](https://github.com/GATB/bcalm#readme)
- [BCOOL](https://github.com/Malfoy/BCOOL#readme)
- [BFC](https://github.com/lh3/bfc#readme)
- [BGREAT2](https://github.com/Malfoy/BGREAT2#readme)
- [EMA](http://ema.csail.mit.edu)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [GenomeScope](https://github.com/schatzlab/genomescope#readme)
- [KmerGenie](http://kmergenie.bx.psu.edu/)
- [lh3/gfa1](https://github.com/lh3/gfa1#readme)
- [Long Ranger](https://10xgenomics.com)
- [Nanopolish](https://github.com/jts/nanopolish#readme)
- [ntCard](https://github.com/bcgsc/ntcard#readme)
- [NxTrim](https://github.com/sequencing/NxTrim#readme)
- [Pilon](https://github.com/broadinstitute/pilon#readme)
- [Porechop](https://github.com/rrwick/Porechop#readme)
- [QUAST](https://github.com/ablab/quast#readme)
- [Racon](https://github.com/isovic/racon#readme)
- [SGA](https://github.com/jts/sga#readme)
- [Tigmint](https://github.com/bcgsc/tigmint#readme)
- [Trimadap](https://github.com/lh3/trimadap#readme)
- [Unicycler](https://github.com/rrwick/Unicycler#readme)

A tool may combine multiple assembly stages in a single tool.

## Preprocess reads

### Illumina mate-pair

- NxTrim

### Illumina paired-end

- Trimadap

### Linked reads

- EMA `ema preprocess`
- Long Ranger `longranger basic`

### Nanopore

- Porechop

## Quality control

- FastQC
- GenomeScope
- KmerGenie
- ntCard
- SGA `sga preqc`

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

## Collapse variants

- ABySS `PopBubbles | MergeContigs`
- lh3/gfa1 `gfaview -b`

## Map reads

- ABySS `abyss-map`
- BGREAT2 `bgreat`
- BWA `bwa mem`
- EMA `ema align` for linked reads
- Long Ranger `longranger align` for linked reads
- Minimap2 `minimap2`
- Unicycler `unicycler_align`

## Estimate copy number

- SGA `sga-astat.py`

## Link unitigs

- ABySS `abyss-fixmate | DistanceEst` for paired-end and mate-pair reads
- ABySS `abyss-longseqdist` for long reads
- ARCS `arcs` for linked reads

## Order and orient

- ABySS `abyss-scaffold` or `SimpleGraph | MergePaths`
- SGA `sga scaffold`

## Contract paths

- ABySS `MergeContigs`
- SGA `sga scaffold2fasta`

## Fill gaps

- ABySS `abyss-sealer`

## Polish

- Nanopolish for Nanopore reads
- Pilon for short reads
- Racon for long reads
- Tigmint to correct large-scale misassemblies with linked reads

# Visualize the assembly graph

- Bandage

## Assess assembly quality

- ABySS `abyss-fac` and `abyss-samtobreak`
- QUAST

# Pipelines

The data is from Unicycler: "These are synthetic reads from plasmids A, B and E from the Shigella sonnei 53G genome assembly". [Shigella sonnei plasmids (synthetic reads)](https://github.com/rrwick/Unicycler/tree/master/sample_data#shigella-sonnei-plasmids-synthetic-reads), [short_reads_1.fastq.gz](https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_1.fastq.gz), [short_reads_2.fastq.gz](https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_2.fastq.gz)

## ABySS

Assemble paired-end reads using ABySS. This ABySS pipeline is a minimal subset of tools run by the complete `abyss-pe` pipeline.

```sh
cd components
k=99
./download_test_data.sh
./abyss_unitigs.sh 0_pe.fq.gz $k
./abyss_contigs_from_unitigs.sh 1_unitig.gfa2 1_unitig.fa $k
./abyss_scaffolding.sh 6_contigs.fa 6_contigs.gfa $k
# Visualize the assembly graph
Bandage load 9_assembly.gfa1 &
```

## BCALM+ABySS

This assembles reads using BCALM for unitigs, then uses the rest of the ABySS pipeline from the previous example.

```sh
cd components
k=99
./download_test_data.sh
./bcalm.sh 0_pe.fq.gz $k
./abyss_contigs_from_unitigs.sh 1_unitig.gfa2 1_unitig.fa $k
./abyss_scaffolding.sh 6_contigs.fa 6_contigs.gfa $k
# Visualize the assembly graph
Bandage load 9_assembly.gfa1 &
```

# Components

## Download test data
 
<!-- markdown-exec(cmd:echo " components/download_test_data.sh" && echo "\`\`\`sh" && cat components/download_test_data.sh && echo "\`\`\`" && echo " --"  ) -->components/download_test_data.sh
```sh
# Download Unicycler test data
# input: nothing
# output: 0_pe.fq.gz

#install dependencies
brew install seqtk curl

seqtk mergepe <(curl -Ls https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_1.fastq.gz) <(curl -Ls https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_2.fastq.gz) | gzip >0_pe.fq.gz
```
 --<!-- /markdown-exec -->

## From reads to unitigs with BCALM

<!-- markdown-exec(cmd:echo "components/bcalm.sh" && echo "\`\`\`sh" && cat components/bcalm.sh && echo "\`\`\`"  && echo " --" ) -->components/bcalm.sh
```sh
# input: [reads] [k]
# e.g. 0_pe.fq.gz 99
#
# output: 1_unitig.gfa2
# contains unitigs created by BCALM

# make sure bcalm and convertToGFA.py are in your path (will later be automatically done by 'conda install bcalm' when someone makes bcalm availaon conda)

reads=$1
k=$2

# Install the dependencies
pip install gfapy
# Unitig with BCALM
bcalm -in $1 -out 1_bcalm -kmer-size $k -abundance-min 1 -verbose 0
mv 1_bcalm.unitigs.fa 1_unitig.fa
convertToGFA.py 1_unitig.fa 1_unitig.gfa $k
# convert bcalm output to gfa2
(printf "H\tVN:Z:1.0\n"; tail -n +2 1_unitig.gfa) >1_unitig.gfa1
gfapy-convert 1_unitig.gfa1 > 1_unitig.gfa2

# cleanup
rm -f 1_bcalm*glue* 1_bcalm.h5
```
 --<!-- /markdown-exec -->

## From reads to unitigs with ABySS

<!-- markdown-exec(cmd:echo "components/abyss_unitigs.sh" && echo "\`\`\`sh" && cat components/abyss_unitigs.sh && echo "\`\`\`"  && echo " --" ) -->components/abyss_unitigs.sh
```sh
# input: [reads] [k]
# e.g. 0_pe.fq.gz 99
#
# output: 1_unitig.gfa2 1_unitig.fa
# contains unitigs created by abyss

reads=$1
k=$2

# setting up
brew install abyss

# Unitig
gunzip -c $reads | ABYSS -k$k -t0 -c0 -b0 -o 1_unitig.fa -
AdjList --gfa2 -k$k 1_unitig.fa >1_unitig.gfa 
mv 1_unitig.gfa 1_unitig.gfa2
```
 --<!-- /markdown-exec -->


## From unitigs to contigs with ABySS 

<!-- markdown-exec(cmd:echo "components/abyss_contigs_from_unitigs.sh"  && echo "\`\`\`sh" && cat components/abyss_contigs_from_unitigs.sh && echo "\`\`\`" && echo " --"   ) -->components/abyss_contigs_from_unitigs.sh
```sh
# input: [unitigs.gfa2] [unitigs.fa] [k]
# e.g. 1_unitig.gfa2 1_unitigs.fa 100
#
# output: 6_contigs.gfa2 and 6_contigs.fa
# contigs produced by ABySS

unitigs_gfa2=$1
unitigs_fa=$2
k=$3

# install dependencies
brew install abyss pigz samtools

# Denoise
abyss-filtergraph --gfa2 -k$k -t200 -c3 -g 2_denoise.gfa $unitigs_gfa2
# Collapse variants
PopBubbles --gfa2 -k$k -p0.99 -g 3_debulge.gfa $unitigs_fa 2_denoise.gfa >3_debulge.path
MergeContigs --gfa2 -k$k -g 3_debulge.gfa -o 3_debulge.fa $unitigs_fa 2_denoise.gfa 3_debulge.path
# Map reads
gunzip -c 0_pe.fq.gz | abyss-map - 3_debulge.fa | pigz >3_debulge.sam.gz
# Link unitigs
gunzip -c 3_debulge.sam.gz | abyss-fixmate -h 4_link.tsv | samtools sort -Osam | DistanceEst --dist -k$k -s500 -n1 4_link.tsv >4_link.dist
# Resolve repeats
samtools faidx 3_debulge.fa
SimpleGraph -k$k -s500 -n5 -o 5_resolverepeats-1.path 3_debulge.gfa 4_link.dist
MergePaths -k$k -o 5_resolverepeats.path 3_debulge.fa.fai 5_resolverepeats-1.path
# Contract paths
MergeContigs --gfa2 -k$k -g 6_contigs.gfa2 -o 6_contigs.fa 3_debulge.fa 3_debulge.gfa 5_resolverepeats.path

# cleanup (comment to keep files)
rm -f 5_resolverepeats.path 3_debulge.gfa 3_debulge.fa 3_debulge.fa.fai 3_debulge.path 5_resolverepeats-1.path 4_link.dist 4_link.tsv 2_denoise.gfa 3_debulge.sam.gz 
```
 --<!-- /markdown-exec -->

## From unitigs to contigs with gfaview

<!-- markdown-exec(cmd:echo "components/gfaview_contigs_from_unitigs.sh"  && echo "\`\`\`sh" && cat components/gfaview_contigs_from_unitigs.sh && echo "\`\`\`"  ) -->components/gfaview_contigs_from_unitigs.sh
```sh
# input: [unitigs_gfa1] [unitigs_gfa2] [unitigs_fa] [k] 
# e.g. 1_unitig.gfa1

# output: 6_contigs.gfa 6_contigs.fa
# 'contigs' generated by popping bubbles and removing tips using gfa1 tool
# https://github.com/lh3/gfa1 

unitigs_gfa1=$1
unitigs_gfa2=$2
unitigs_fa=$3
k=$4

# remove tips (2 rounds), misc trimming, pop bubbles, 
# but then cannot use -u to build unitigs of simplified graph (i.e. contigs), because it doesn't output sequences
gfaview -t -m -t -b $unitigs_gfa1 > 2_denoised.gfa1

# make the output of gfaview properly gfa1 and gfa2
sed -i 's/L1/l1/g' 2_denoised.gfa1 # gfapy wants lower case flags
sed -i 's/L2/l2/g' 2_denoised.gfa1
sort -n -k 2 2_denoised.gfa1 > 2_denoised.gfa1.sorted # because ABySS GFA parser wants S lines in same order as in FASTA file
head -n 1 $unitigs_gfa2 > 2_denoised.gfa2 # copy header as gfaview doesnt 
gfapy-convert 2_denoised.gfa1.sorted >> 2_denoised.gfa2

# would like to use use rgfa to do compaction but itdoesnt work on the Singella example
#./rgfa-mergelinear 2_denoised.gfa1 > 6_contigs.gfa

# convert gfa to fasta
awk '$1 ~/S/ {print ">"$2"\n"$3}' 2_denoised.gfa1.sorted > 2_denoised.fa

# just contract paths in the 2_denoised.gfa2, using ABySS as wasn't able to do it with neither gfaview nor rgfa-mergelinear
abyss-filtergraph --gfa2 -k$k --assemble 2_denoised.gfa2 2_denoised.fa -g 3_denoised.gfa2 > 3_denoised.paths
MergeContigs --gfa2 -k$k 2_denoised.fa 3_denoised.gfa2 3_denoised.paths -o 6_contigs.fa -g 6_contigs.gfa2 

# cleanup (comment to keep files)
rm -f 2_denoised.gfa1.sorted 2_denoise.gfa1 2_denoised.gfa1 3_denoised.paths 3_denoised.gfa2
```
<!-- /markdown-exec -->

## From contigs to scaffolds with ABySS

<!-- markdown-exec(cmd:echo "components/abyss_scaffolding.sh"  && echo "\`\`\`sh" && cat components/abyss_scaffolding.sh && echo "\`\`\`"  ) -->components/abyss_scaffolding.sh
```sh
# input: [contigs.fa] [contigs.gfa] [k]
# e.g. 6_contigs.fa 6_contigs.gfa 99
# 
# output: 9_assembly.gfa and 9_assembly.gfa1
# scaffolds assembly using abyss scaffolder in GFA2 and GFA1 format

contigs_fa=$1
contigs_gfa=$2
k=$3

# Map reads
gunzip -c 0_pe.fq.gz | abyss-map - $contigs_fa | pigz >6_contigs.sam.gz
# Link unitigs
gunzip -c 6_contigs.sam.gz | abyss-fixmate -h 7_link.tsv | samtools sort -Osam | DistanceEst --dot -k$k -s500 -n1 7_link.tsv >7_link.gv
# Order and orient
abyss-scaffold -k$k -s500-1000 -n5-10 $contigs_gfa 7_link.gv >8_scaffold.path
# Contract paths
MergeContigs --gfa2 -k$k -g 9_assembly.gfa -o 9_assembly.fa $contigs_fa $contigs_gfa 8_scaffold.path
# Compute assembly metrics
abyss-fac 9_assembly.fa
# Convert GFA2 to GFA1
abyss-todot --gfa1 9_assembly.gfa >9_assembly.gfa1

# cleanup (comment to remove)
rm -f 6_contigs.sam.gz 7_link.gv 7_link.tsv 8_scaffold.path 8_scaffold.path
```
<!-- /markdown-exec -->
