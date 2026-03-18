# Usage:
# bash merge_all_pairs.sh ../Raw_data/chunk1_batch1/

raw_dir="../Raw_data/"
unm_dir="../Unmerged_reads"
fwd_end="_1.fq.gz"
rev_end="_2.fq.gz"

# addition for chunks:
out_dir=${1/$raw_dir/""} # edited for chunks
mkdir -p $out_dir

merge_all ()
{
  for forward_in in $raw_dir/$out_dir/*$fwd_end # edited for chunks
  do
    reverse_in=${forward_in/"$fwd_end"/"$rev_end"}
    merged_out=${forward_in/"$raw_dir/"/""}
    output=${merged_out/"$fwd_end"/""}
    echo $output
#     ~/tools/BBTools/bbmap/bbmerge.sh in1=$forward_in in2=$reverse_in out=$output".merged.fq.gz" outu1=$unm_dir"/"$output".unmerged_1.fq.gz" outu2=$unm_dir"/"$output".unmerged_2.fq.gz"
    ~/tools/BBTools/bbmap/bbmerge.sh in1=$forward_in in2=$reverse_in out="./"$output".merged.fq.gz" outu1=$unm_dir"/"$output".unmerged_1.fq.gz" outu2=$unm_dir"/"$output".unmerged_2.fq.gz" # edited for chunks
  done
}

# Write stout and sterr to log file:
merge_all 2>&1 | tee bbmerge_log.out
