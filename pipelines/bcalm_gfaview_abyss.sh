cd components
k=99
./download_test_data.sh
./bcalm.sh 0_pe.fq.gz $k
./gfaview_contigs_from_unitigs.sh 1_unitig.gfa1 1_unitig.gfa2 1_unitig.fa $k
./abyss_scaffolding.sh 6_contigs.fa 6_contigs.gfa2 $k

