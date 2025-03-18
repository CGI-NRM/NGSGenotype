# Usage (note, it assumes that subfolders are prefixed with "chunk_"):
# bash loop_over_chunks.sh [folder containing chunk folders]

for cur_folder in $1/chunk_*
do
  chunk_name=${cur_folder/*\/}
  cur_time=`python -c "import datetime; cur_time = datetime.datetime.now(); print(f'{cur_time.hour}:{cur_time.minute}:{cur_time.second}')"`
  echo "Parsing "$chunk_name". ("$cur_time")"
  python ./snpotypewriter.py $cur_folder > $chunk_name"_genotypes.csv"
  echo " - Done."
done
cur_time=`python -c "import datetime; cur_time = datetime.datetime.now(); print(f'{cur_time.hour}:{cur_time.minute}:{cur_time.second}')"`
echo "All parsing done. ("$cur_time")"
