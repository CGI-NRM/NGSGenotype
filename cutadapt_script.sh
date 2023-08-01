for file in Raw_data/*.fq.gz
do
	sample=$(basename -s .fq.gz "$file") 
	cutadapt --discard-untrimmed -e 0.1 --no-indels -g ^file:primers.fa -o Filtered_data/{name}_$sample.fastq.gz $file
done
