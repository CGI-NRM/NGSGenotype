# Specify primer file:
primer="Primers/primer_file.fa"

# Trim primers and create a file for each locus:
for file in Merged_data/*.fq.gz
do
	sample=$(basename -s .fq.gz "$file")
	cutadapt --discard-untrimmed -e 0.1 --no-indels -g ^file:$primer -o Filtered_data/{name}_$sample.fastq.gz $file
done
