# Usage:
# bash filter_out_SNPs.sh
# bash filter_out_SNPs.sh chunk_[n]

mkdir -p Filtered_data/$1/SNP_filtered

filter_snps() {
  while read locus
  do
    read sequence
    locus_name=${locus/'>'/''}
    for fastq in Filtered_data/$1/*$locus_name*
    do
      cutadapt -j 0 --discard-untrimmed -e 0.1 --no-indels --error-rate 0.1 --maximum-length 1 -a $sequence -o ${fastq/'Filtered_data/'$1/'Filtered_data/'$1'/SNP_filtered/'} $fastq
    done
  done < ../make_snp_primers/Bear_SNP_surrounding_bases_current_v1.fa
}

filter_snps $1 | tee cutadapt_snp_trim_log_${1//'/'/''}.out