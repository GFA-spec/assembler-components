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
