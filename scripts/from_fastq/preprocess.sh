#!/bin/bash

shopt -s nullglob #allows create an empty array
homePath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #cur. script locat.
source "$HOME/software/bash_scripts/funcList.sh"

curScrName=${0##*/} #delete all before last backSlash
EchoLineBold
echo "[Start] $curScrName"
EchoLineSh
lenStr=${#curScrName}
lenStr=$((25 + lenStr))
printf "%-${lenStr}s %s\n"\
       "The location of $curScrName:"\
       "$homePath"
printf "%-${lenStr}s %s\n"\
       "The $curScrName is executed from:"\
       "$PWD"
EchoLineSh

# inpPairEnd - join different lanes using ,
inpPairEnd1=${1:-"A_S1_L001_R1_001.fastq.gz,A_S1_L002_R1_001.fastq.gz"}
inpPairEnd2=${2:-"A_S1_L001_R2_001.fastq.gz,A_S1_L002_R2_001.fastq.gz"}
argsFile=${3:-"scRna.args"}

posArgs=("dataPath" "outPath" "alignIndDir" "geneGtf" "geneBed" "geneFa"
         "chromSize" "nCpu" "stages")
ReadArgs "$argsFile" 1 "" ${#posArgs[@]} "${posArgs[@]}"
PrintArgs "" "${posArgs[@]} inpPairEnd1 inpPairEnd2"

# if [[ -n "$dataPath" ]]; then
#     if [[ -n "$inpPairEnd1" ]]; then
#         inpPairEnd1="$dataPath/$inpPairEnd1"
#     fi

#     if [[ -n "$inpPairEnd1" ]]; then
#         inpPairEnd2="$dataPath/$inpPairEnd2"
#     fi
# fi

stageQC=1
stageAlign=2
stageSort=3
stageAlignQC=4
stageExpress=5

mkdir -p "$outPath"
experName=${inpPairEnd1%%_L[0-9]*}
experName=${experName##*/}

tmpDir=$(mktemp -q -u /store01/tmp/scRNAXXXXXX)
## Quality Control
# Can use FASTQC, but seems useless
if [[ $(IsInInt "$stages" $stageQC) -eq 1 ]]; then   
    printf "QC ... \n"
    qcDir="$outPath/QC"
    mkdir -p "$qcDir"
    fastqc -o "$qcDir" \
           $(echo "$inpPairEnd1" | tr ',' ' ') \
           $(echo "$inpPairEnd2" | tr ',' ' ')

    if [[ $? -ne 0 ]]; then
        WarnMsg "QC failed"
    else
      printf "done\n"
    fi
fi


## Alignment
alignDir="$outPath/align"
alignOut="$alignDir/${experName}_"

if [[ $(IsInInt "$stages" $stageAlign) -eq 1 ]]; then
    printf "Alignment ... \n"
    mkdir -p "$alignDir"

    STAR --runThreadN $nCpu \
         --runMode alignReads \
         --readFilesIn "$inpPairEnd1" "$inpPairEnd2" \
         --readFilesCommand zcat \
         --genomeDir "$alignIndDir" \
         --outFileNamePrefix "$alignOut" \
         --outSAMtype BAM SortedByCoordinate \
         --outTmpDir "$tmpDir"

    #--outSAMtype BAM SortedByCoordinate \

    if [[ $? -ne 0 ]]; then
        ErrMsg "Alignment failed" "$?"
    else
      printf "done\n"
      rm -rf "$alignTmpDir"
    fi
fi

alignOut="${alignOut}Aligned.sortedByCoord.out.bam"

## Transforming to bam and sorting with respect to coordinates
# if [[ $(IsInInt "$stages" $stageSort) -eq 1 ]]; then
#     alignIn="$alignOut"
#     alignOut="${alignIn%.sam*}.srt.bam"
#     echo "$alignIn $alignOut"
#     samtools sort "$alignIn" -o "$alignOut" -O bam -T "$tmpDir"
# fi


## Quality Control of Mapped reads
if [[ $(IsInInt "$stages" $stageAlignQC) -eq 1 ]]; then
samtools index "$alignOut"
    printf "Gene body coverage ... "
    alignDirQC="$outPath/alignQC"
    mkdir -p "$alignDirQC"

    # Gene body coverage
    python2 ~/software/RSeQC-2.6.4/scripts/geneBody_coverage.py -i "$alignOut"\
            -r "$geneBed"  -o "$alignDirQC/$experName"

    if [[ $? -ne 0 ]]; then
        WarnMsg "QC of mapped reads failed"
    else
      printf "done\n"
    fi
    

    # # Stat
    # printf "QC of mapped reads ... "
    # statQC="$alignDirQC/$experName.stat"
    # python2 ~/software/RSeQC-2.6.4/scripts/bam_stat.py -i "$alignOut" > "$statQC"

    # if [[ $? -ne 0 ]]; then
    #     WarnMsg "QC of mapped reads failed"
    # else
    #   printf "done\n"
    # fi
fi


## Expression matrix
if [[ $(IsInInt "$stages" $stageExpress) -eq 1 ]]; then
    printf "Expression matrix ... "
    expresDir="$outPath/express"
    mkdir -p "$expresDir"
    expressOut="$expresDir/$experName.txt"

    # htseq-count --nonunique all \
    #             -r pos \
    #             -f bam \
    #             "$alignOut"\
    #             "$geneGtf" > "$expressOut"

    featureCounts -O -M -Q 30 -p -a "$geneGtf" -o "$expressOut" "$alignOut"

    if [[ $? -ne 0 ]]; then
        ErrMsg "Quantification failed" "$?"
    else
      print "done\n"
    fi
fi
