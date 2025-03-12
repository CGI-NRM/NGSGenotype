mkdir -p Filtered_data/SNP_filtered

filter_snps() {
  while read locus
  do
    read sequence
    locus_name=${locus/'>'/''}
    for fastq in Filtered_data/*$locus_name*
    do
      cutadapt -j 0 --discard-untrimmed -e 0.1 --no-indels --error-rate 0.1 --maximum-length 1 -a $sequence -o ${fastq/'Filtered_data/'/'Filtered_data/SNP_filtered/'} $fastq
    done
  done < Primers/Bear_SNP_surrounding_bases.fa
}

filter_snps | tee cutadapt_snp_trim_log.out
