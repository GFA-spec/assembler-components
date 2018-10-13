# input: [unitigs.gfa2] [unitigs.fa] [k]
# e.g. 1_unitig.gfa2 1_unitigs.fa 100
#
# output: 6_contigs.gfa and 6_contigs.fa
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
MergeContigs --gfa2 -k$k -g 6_contigs.gfa -o 6_contigs.fa 3_debulge.fa 3_debulge.gfa 5_resolverepeats.path

# cleanup (comment to keep files)
rm -f 5_resolverepeats.path 3_debulge.gfa 3_debulge.fa 3_debulge.fa.fai 3_debulge.path 5_resolverepeats-1.path 4_link.dist 4_link.tsv 2_denoise.gfa 3_debulge.sam.gz 
