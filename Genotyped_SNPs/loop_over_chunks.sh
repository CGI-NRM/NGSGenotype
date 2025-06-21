# Usage (note, it assumes that subfolders are prefixed with "chunk_"):
# bash loop_over_chunks.sh ../Filtered_data/

for cur_folder in $1/chunk*
do
  chunk_name=${cur_folder/*\/}
  cur_time=`python -c "import datetime; cur_time = datetime.datetime.now(); print(f'{cur_time.hour}:{cur_time.minute}:{cur_time.second}')"`
  echo "Parsing "$chunk_name". ("$cur_time")"
  snp_folder=$cur_folder"/SNP_filtered/"
  if [[ -d $snp_folder ]]
  then
    echo " - Genotyping files in "$snp_folder"."
    python ./snpotypewriter.py $snp_folder > $chunk_name"_genotypes.csv"
  else
    echo " - No SNP folder found."
  fi
  echo " - Done."
done
cur_time=`python -c "import datetime; cur_time = datetime.datetime.now(); print(f'{cur_time.hour}:{cur_time.minute}:{cur_time.second}')"`
echo "All parsing done. ("$cur_time")"
