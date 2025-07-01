#!/bin/bash

tools="dss. methylkit. methylSig. generate_DMR. dsseq."

for file in triple_*/*.CpG.not.null.merged.CpG.N.100.trimmed.bed; do
  filename=$(basename "$file")
  dir=$(dirname "$file")

  tool1=""
  tool2=""
  tool3=""
  for tool in $tools; do
    if [[ "$filename" == *"$tool"* ]]; then
      if [ -z "$tool1" ]; then
        tool1="$tool"
      elif [ -z "$tool2" ]; then
	tool2="$tool"
      else
        tool3="$tool"
        break
      fi
    fi
  done

  tool1=${tool1/generate_DMR/DiffMethylTools}
  tool2=${tool2/generate_DMR/DiffMethylTools}
  tool3=${tool3/generate_DMR/DiffMethylTools}
  newname="${tool1}_${tool2}_${tool3}.bed"
  newpath="${dir}/${newname}"

  echo "Copying '$file' -> '$newpath'"
  cp "$file" "$newpath"
done
