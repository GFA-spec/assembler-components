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
