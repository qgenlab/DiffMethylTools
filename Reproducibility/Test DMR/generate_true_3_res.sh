#!/bin/bash

tools=("dss." "methylkit." "methylSig." "DiffMethylTools." "dsseq.")
search_dir=$1  

for dataset in "${tools[@]}"; do
  echo "Processing: $dataset"
  files_to_merge=()
  for t1 in "${tools[@]}"; do
    for t2 in "${tools[@]}"; do
      for t3 in "${tools[@]}"; do
      if [[ "$t1" != "$t2" && "$t1" != "$t3" && "$t2" != "$t3" ]]; then
        if [[ "$t1" != "$dataset" && "$t2" != "$dataset" && "$t3" != "$dataset" ]]; then
          filename="${t1}_${t2}_${t3}.bed"
          fullpath=$(find "$search_dir" -type f -name "$filename")
          if [[ -n "$fullpath" ]]; then
            files_to_merge+=("$fullpath")
          fi
        fi
      fi
    done
  done
done
  if [ ${#files_to_merge[@]} -gt 0 ]; then
    echo "Merging ${#files_to_merge[@]} files into ${dataset}.test.bed"
    cat "${files_to_merge[@]}" | bedtools sort -i - | bedtools merge -i - > "${search_dir}/${dataset}.test.bed"
  else
    echo "No files to merge for dataset $dataset"
  fi
done
