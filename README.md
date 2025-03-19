# NGSGenotype
Analysis pipeline for genotyping individuals based on either microsatellites or SNPs (Single Nucleotide Polymorphisms).

## Usage
Clone down this repository and move into it:
```bash
git clone https://github.com/CGI-NRM/NGSGenotype
cd NGSGenotype/
```

Move your data into the Raw_data-folder:
```bash
mv [folder with your data in it]/*.fastq.gz ./Raw_data/
```

If you have paired end data, merge reads with (log will be saved to 'bbmerge_log.out'):
```bash
cd Merged_data/
bash merge_all_pairs.sh
```

Trim all primers, choose primer file as well as if data was merged of single end (log will be saved to 'cutadapt_log.out'):
```bash
cd ..
bash cutadapt_script.sh merged Primers/Bear_SNP_all_loci.fa # for merged data with bear primers

# or:
bash cutadapt_script.sh single Primers/L_vulgaris.fa # for single end data with Lissotriton primers
```

## Microsatellites
Continue analysis in R using the rmarkdown file 'ngsgenotype.rmd'.

## Single Nucleotide Polymorphisms
First trim the superflous sequences around the SNPs:
```bash
bash filter_out_SNPs.sh # currently hardcoded for bear SNPs
```

Then, either continue in R using the rmarkdown 'ngsgenotype_snp.rmd', or first genotype the SNPs by:
```bash
cd Genotyped_SNPs/
python snpotypewriter.py ../Filtered_data/SNP_filtered/ > genotyped_SNP_data.csv

# or, if data is divided into chunks (where chunk folder names start with "chunk_"):
bash loop_over_chunks.sh ../path/to/folder/with/chunks/
```

And then continue in R.
