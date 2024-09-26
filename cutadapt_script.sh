# Specify primer file:
primer=$1 # run script with the primer file as an argument

# Trim primers and create a file for each locus:
for file in Raw_data/*.fastq.gz
do
	sample=$(basename -s .fastq.gz "$file")
	cutadapt --discard-untrimmed -e 0.1 --no-indels --error-rate 0.1 -a file:$primer -o Filtered_data/{name}_$sample.fastq.gz $file
done
