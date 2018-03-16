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

stageQC=0
stageCat=1
stageWhitelist=2
stageExtractBC=3
stageAlign=4
stageExpress=5
stageCountMolec=6
stageAlignQC=7

bcPattern='(?P<cell_1>.{16})(?P<umi_1>.{10})'

mkdir -p "$outPath"
experName=${inpPairEnd1%%_L[0-9]*}
experName=${experName##*/}
outPath="$outPath/$experName"

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

## Concatenation
inpPairEnd1Tmp="$outPath/cat/${experName}_read1.fastq.gz"
inpPairEnd2Tmp="$outPath/cat/${experName}_read2.fastq.gz"
if [[ $(IsInInt "$stages" $stageCat) -eq 1 ]]; then
    mkdir -p "$(dirname $inpPairEnd1Tmp)"
    printf "Concatenation of lanes ... "
    readarray -t inpPairEnd1 <<< "$(echo "$inpPairEnd1" | tr ',' '\n')"
    readarray -t inpPairEnd2 <<< "$(echo "$inpPairEnd2" | tr ',' '\n')"

    cat "${inpPairEnd1[@]}" > $inpPairEnd1Tmp
    cat "${inpPairEnd2[@]}" > $inpPairEnd2Tmp

    if [[ $? -ne 0 ]]; then
        ErrMsg "One of concatenation failed" "$?"
    else
      printf "done\n"
      rm -rf "$alignTmpDir"
    fi
fi
inpPairEnd1=("$inpPairEnd1Tmp")
inpPairEnd2=("$inpPairEnd2Tmp")

## Whitelist
whiteList="$outPath/stat/${experName}_whitelist.txt"
if [[ $(IsInInt "$stages" $stageWhitelist) -eq 1 ]]; then
    mkdir -p "$(dirname $whiteList)"
    whiteListPlot="${whiteList%.*}"
    printf "Create a whitelist ... "
    # Umi tool predict number of cells
    # use --expect-cells= - if i know the expectedn umber of cells
    # or --set-cell-number = if it is exact
    umi_tools whitelist\
              --stdin "$inpPairEnd1" \
              --extract-method regex \
              --bc-pattern "$bcPattern" \
              --plot-prefix "$whiteListPlot" \
              -L "${whiteList%.*}.log" \
              > "$whiteList"
    if [[ $? -ne 0 ]]; then
        ErrMsg "Whitelist failed" "$?"
    else
      printf "done\n"
    fi
    whiteListPlot="${whiteListPlot}_cell_barcode_counts_density.png"
    echo "Check: $whiteListPlot"
fi

## Extract barcode and filter the reads
inpPairEnd1Tmp="$outPath/final_reads/${experName}_read1.fastq.gz"
inpPairEnd2Tmp="$outPath/final_reads/${experName}_read2.fastq.gz"
if [[ $(IsInInt "$stages" $stageExtractBC) -eq 1 ]]; then
    printf "Extract barcode ... "
    mkdir -p "$(dirname $inpPairEnd1Tmp)"
    umi_tools extract\
              --extract-method regex \
              --bc-pattern "$bcPattern" \
              --stdin "$inpPairEnd1" \
              --stdout "$inpPairEnd1Tmp" \
              --read2-in "$inpPairEnd2" \
              --read2-out "$inpPairEnd2Tmp" \
              --filter-cell-barcode \
              --whitelist "$whiteList" \
              -L "${whiteList%.*}.extracted.log"
    if [[ $? -ne 0 ]]; then
        ErrMsg "Extract barcode failed" "$?"
    else
      #rm -rf /tmp/*
      printf "done\n"
    fi
fi
inpPairEnd1=("$inpPairEnd1Tmp")
inpPairEnd2=("$inpPairEnd2Tmp")


## Alignment
alignOut="$outPath/align/${experName}_"
if [[ $(IsInInt "$stages" $stageAlign) -eq 1 ]]; then
    printf "Alignment ... \n"
    mkdir -p "$(dirname $alignOut)"

    STAR --runThreadN $nCpu \
         --runMode alignReads \
         --readFilesIn "$inpPairEnd2" \
         --readFilesCommand zcat \
         --genomeDir "$alignIndDir" \
         --outFileNamePrefix "$alignOut" \
         --outFilterMultimapNmax 1 \
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


## Expression matrix
expressOut="$outPath/expres/$experName.txt"
expressOutSort="${expressOut%.*}.srt.bam"
if [[ $(IsInInt "$stages" $stageExpress) -eq 1 ]]; then
    printf "Expression matrix ... "
    mkdir -p "$(dirname $expressOut)"
    #featureCounts -O -M -Q 30 -p -a "$geneGtf" -o "$expressOut" "$alignOut"
    featureCounts -R BAM -T "$nCpu" -a "$geneGtf" -o "$expressOut" "$alignOut"

    if [[ $? -ne 0 ]]; then
        ErrMsg "Quantification failed" "$?"
    else
      print "done\n"
    fi
fi


## Counting molecules
expressCountOut="${expressOut%.*}.counts.tsv.gz"
expressCountMatrix="${expressOut%.*}.matrix.tsv.gz"
if [[ $(IsInInt "$stages" $stageCountMolec) -eq 1 ]]; then
    printf "Sorting ..."
    countInp="$(dirname $expressOut)/$(basename $alignOut).featureCounts.bam"
    countOut="${countInp%.*}.sort.bam"
    samtools sort "$countInp" -o "$countOut"
    samtools index "$countOut"
    if [[ $? -ne 0 ]]; then
        ErrMsg "Sorting failed" "$?"
    else
      print "done\n"
    fi

    printf "Count molecules ..."
    umi_tools count --per-gene \
              --gene-tag=XT \
              --per-cell \
              -I "$countOut" \
              -S "$expressCountOut" \
              -L "${expressCountOut%.*}.log"
    if [[ $? -ne 0 ]]; then
        ErrMsg "Count molecules failed" "$?"
    else
      print "done\n"
    fi


    printf "Count molecules matrix ..."
    umi_tools count --per-gene\
              --gene-tag=XT \
              --per-cell \
              --wide-format-cell-counts \
              -I "$countOut" \
              -S "$expressCountMatrix"\
              -L "${expressCountMatrix%.*}.log"
    if [[ $? -ne 0 ]]; then
        ErrMsg "Count molecules matrix failed" "$?"
    else
      print "done\n"
    fi
fi


## Quality Control of Mapped reads
if [[ $(IsInInt "$stages" $stageAlignQC) -eq 1 ]]; then
    printf "Create index ... "
    samtools index "$alignOut"
    if [[ $? -ne 0 ]]; then
        ErrMsg "Creating index failed" "$?"
    else
      print "done\n"
    fi
    
    printf "Gene body coverage ... "
    alignDirQC="$outPath/alignQC"
    mkdir -p "$alignDirQC"

    # Gene body coverage
    python2 ~/software/RSeQC-2.6.4/scripts/geneBody_coverage.py\
            -i "$alignOut" \
            -r "$geneBed" \
            -o "$alignDirQC/$experName"

    if [[ $? -ne 0 ]]; then
        WarnMsg "Genebody of mapped reads failed"
    else
      printf "done\n"
    fi
    
    # Stat
    printf "Stat of mapped reads ... "
    statQC="$alignDirQC/$experName.stat"
    python2 ~/software/RSeQC-2.6.4/scripts/bam_stat.py \
            -i "$alignOut" > "$statQC"

    if [[ $? -ne 0 ]]; then
        WarnMsg "Stat of mapped reads failed"
    else
      printf "done\n"
    fi
fi
