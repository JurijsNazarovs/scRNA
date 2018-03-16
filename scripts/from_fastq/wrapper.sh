#!/bin/bash

homePath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #cur. script locat.
source "$HOME/software/bash_scripts/funcList.sh"

argsFile=${1:-"scRna.args"}
posArgs=("dataPath")
ReadArgs "$argsFile" 1 "" ${#posArgs[@]} "${posArgs[@]}"
PrintArgs "" "${posArgs[@]}"

## Get all input files from directory
readarray -t inpFiles <<< "$(find $dataPath -maxdepth 1 -type f -name "*gz")"
inpFiles=($(echo ${inpFiles[@]%%_L[0-9]*} | tr ' ' '\n' | sort -u))

for i in ${inpFiles[@]}; do
  filePrefix="$(basename $i)"
  echo "$filePrefix - screen created"
  readarray -t inpReadEnd1 <<< "$(find $dataPath -maxdepth 1 -type f -name "$filePrefix*_R1*gz")"
  inpReadEnd1=$(JoinToStr "," "${inpReadEnd1[@]}")

  readarray -t inpReadEnd2 <<< "$(find $dataPath -maxdepth 1 -type f -name "$filePrefix*_R2*gz")"
  inpReadEnd2=$(JoinToStr "," "${inpReadEnd2[@]}")

    # screen -dmS scRna_$filePrefix bash -c \
    #       "$HOME/private/scRna/preprocess_umitools.sh $inpReadEnd1 $inpReadEnd2 $argsFile &> $filePrefix.log"

  bash $HOME/private/scRna/preprocess_umitools.sh "$inpReadEnd1" "$inpReadEnd2" "$argsFile"
done
