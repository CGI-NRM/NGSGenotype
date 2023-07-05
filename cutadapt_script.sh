for file in rawdata/*.fq.gz
do
	sample=$(basename -s .fq.gz "$file") 
	cutadapt --discard-untrimmed -e 0.1 --no-indels -g ^file:primers.fa -o filtered2/{name}_$sample.fastq.gz $file
done
