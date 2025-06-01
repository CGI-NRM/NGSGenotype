filt_folder="../../Filtered_data/SNP_filtered/"
end_of_name="-CQ*"

minimum_seqs=$1
verbose=$2 # if second argument is "verbose", additional information is printed

if ! [[ -n $minimum_seqs ]]
then
  minimum_seqs=1 # if no threshold is given, set it at 1
fi

echo "Number of loci containing at least $minimum_seqs sequences:"

while read parameter
do
  echo ""
  echo "- "$parameter":"
  samples=`for i in $filt_folder*$parameter*.gz ; do i=${i/*_} ; echo ${i/$end_of_name} ; done | sort | uniq`
  for sample in $samples
  do
    counter=0
    for locus in $filt_folder*$sample*
    do
      n_seqs=`zgrep "^@" $locus | wc -l`
      if [[ $n_seqs -ge $minimum_seqs ]]
      then
	current=1
      else
	current=0
	if [[ $verbose == "verbose" ]]
	then
	  echo "- ${locus/"$filt_folder"/""} only has $n_seqs sequences."
	fi
      fi
      counter=$((counter + current))
    done
    echo $counter"/"`ls $filt_folder*$sample* | wc -l` $sample
  done

done < ./filtering_groups.txt
