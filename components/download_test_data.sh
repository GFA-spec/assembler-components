# Download Unicycler test data
# input: nothing
# output: 0_pe.fq.gz

#install dependencies
brew install seqtk curl

seqtk mergepe <(curl -Ls https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_1.fastq.gz) <(curl -Ls https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_2.fastq.gz) | gzip >0_pe.fq.gz
