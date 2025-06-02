# bash seqs_per_locus.sh [folder] [minimum_seqs] [filter_string] > file.txt

loci=`for i in $1/*.gz ; do i=${i/$1/""} ; i=${i/"/"/""} ; i=${i/_*} ; echo $i ; done | sort | uniq`

for locus in $loci
do
  target_files=`ls $1/*$3* | grep $locus`
  n_seqs=`zgrep "^@" $target_files | wc -l`
  if [[ $n_seqs -gt $2 ]]
  then
    echo "$locus $n_seqs"
  fi
done
