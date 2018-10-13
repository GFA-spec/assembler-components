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
