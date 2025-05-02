#!/bin/bash
## Usage function
usage() {
    echo "Usage: HAPLOTYPES -s, --samples SAMPLES_CSV -p, --primers PRIMERS_CSV -o,--output OUTPUT_FOLDER  [-q, --quality QUALITY] 
    [-n, --nreads NUMBER_OF_READS] [-t, --threads NUMBER_OF_THREADS] [-f,--force] [-b,--blast]

    REQUIRED ARGUMENTS:
    -o | --output: Provide the path to the output folder
    -p | --primers: Provide the path to a csv file that contains the IDs, the forward and the reverse sequences of the primers and the expected fragment size. Find an example in the section dedicated to demultiplexing.
    -s | --samples: Provide the path to a csv file that contains the names (without extension) and the path of the fastq files containing the raw reads that will be analysed
    
    OPTIONAL ARGUMENTS:
    -b | --blast: Use BLAST instead of BOLD for species identification by setting this parameter and by providing an email address for using NCBI Entrez (default: disabled)
    -c | --similar_consensus: Change the minimum similarity threshold (in percent) of amplicon_sorter. If the similarity of two clusters are smaller or equal to the threshold, they are considered separartely otherwise they are collapsed prior to consensus reconstruction (default: 96)
    -e | --exclude: Provide a text file with samples and loci to exclude from the analyses. Each row contains one comma-separated sample-ID and locus-ID to be excluded.
    -f | --force: Force overwriting the previous ouput folder (default: Cowardly refusing to overwrite)
    -g | --outgroup: Define an outgroup sample or a comma,separated list of outgroup samples (need to match the IDs in --samples)
    -i | --partition: Use partitionmodel for iqtree with combined dataset (default: disabled)
    -k | --kthreshold: Define the threshold k for the maximum allowed mismatches for the primer alignment during demultiplexing (default: 0.05)
    -m | --minreads: Set the minimum number of reads required for the reconstruction of consensus sequences; If there are not enough reads, the locus will be ignored (default: 100)
    -n | --nreads: Provide the absolute number or the percentage of top quality reads you want the pipeline to take into account for consensus sequences generation and variant calling (default is 500)
    -q | --quality: Provide an integer that represent the minimum quality which the reads will be filtered for (default is 10)
    -r | --sizerange: Define the allowed size buffer around the expected locus length as defined in the primers file (default: 100)
    -t | --threads: Number of threads the program will be using (default: 10)
    -w | --nowatermark: Remove the watermark from the tree figures
    -y | --freqthreshold: Retain consensus sequences for further analyses which are supported by raw reads, whose frequency in the total pool of reads is larger or equal to this threshold (default: 0.1).

  
  Input HAPLOTYPES -h,--help to show the help message"
    exit 1
}

## Initialize variables with default values
quality=7           # -q
similarconsensus=96 # -c
nreads=500          # -n
sizerange=100       # -r
minreads=100        # -m
threads=10          # -t
kthres=0.05         # -k
force="no"          # -f
blast="no"          # -b
partition="no"      # -i
outgroup="no"       # -g
freqthreshold=0.1   # -y


os="$(uname -s)"

## Find basedir
tmp=$(pwd)
wd=$tmp  #${tmp%/*} ######## <- wd is now the base directory below the shell/ folder

# extra check for container execution to account for different paths
if [[ $AMICONTAINER == "yes" ]]; then
	echo "I am in docker"
	ampdir="/app/AmpliPiper"
else
	echo "I am local"
	ampdir=$wd
fi

if
    [[ -d $wd/logs ]] \
        ;
then
    rm -rf $wd/logs
    mkdir $wd/logs
fi

## Test if software already installed if not, start installation
if
    [[ ! -d $ampdir/envs ]] \
        ;
then
    bash $ampdir/shell/setup.sh
fi

## Test if dependencies correctly installed
if
    [[ -f $wd/logs/dep.err ]] \
        ;
then
    printf "\n######################\nError during installation:\n\n"

    # print errors, if there are any and quit
    while read -r line; do
        echo $line
    done <$wd/logs/dep.err
    printf "######################\n\n\n"
    usage
    exit 1
fi

## Test if software correctly installed
if
    [[ -f $wd/envs/logs/setup.err ]] \
        ;
then
    printf "\n######################\nError during installation of the following software packages:\n\n"

    # print errors, if there are any and quit
    while read -r line; do
        echo $line
    done <$ampdir/envs/logs/setup.err
    printf "######################\n\n\n"
    usage
    exit 1
fi

## whereis conda <- initiate conda to load software
eval "$(conda shell.bash hook)"

# Loop through the commandline arguments
while
    [[ $# -gt 0 ]] \
        ;
do
    case "$1" in
    -h | --help)
        usage
        ;;
    -b | --blast)
        blast="$2"
        shift 2
        ;;
    -c | --similar_consensus)
        similarconsensus="$2"
        shift 2
        ;;
    -e | --exclude)
        exclude="$2"
        shift 2
        ;;
    -f | --force)
        force="yes"
        shift 1
        ;;
    -g | --outgroup)
        outgroup="$2"
        shift 2
        ;;
    -i | --partition)
        partition="yes"
        shift 1
        ;;
    -k | --kthreshold)
        kthres="$2"
        shift 2
        ;;
    -m | --minreads)
        minreads="$2"
        shift 2
        ;;
    -n | --nreads)
        nreads="$2"
        shift 2
        ;;
    -o | --output)
        output="$2"
        shift 2
        ;;
    -p | --primers)
        primers="$2"
        shift 2
        ;;
    -q | --quality)
        quality="$2"
        shift 2
        ;;
    -r | --sizerange)
        sizerange="$2"
        shift 2
        ;;
    -s | --samples)
        samples="$2"
        shift 2
        ;;
    -t | --threads)
        threads="$2"
        shift 2
        ;;
    -y | --freqthreshold)
        freqthreshold="$2"
        shift 2
        ;;
    -w | --nowatermark)
        nowatermark="yes"
        shift 2
        ;;
    *)
        echo "Unknown option: $1"
        usage
        ;;
    esac
done


#check if output path is absolute or relative
if [[ "$output" = /* ]]
then
   echo "Output direcory: $output"
else
   output="$wd/$output"
   echo "New Output directory: $output"
fi

## Is the sample input file provided?
if
    [[ -z "${samples}" || ! -f "${primers}" ]] \
        ;
then
    printf "\n######################\nMissing required argument: INPUT FILE with names and paths to raw data in csv format\n######################\n\n\n"
    usage
fi

## Is the primer input file provided?
if
    [[ -z "${primers}"  || ! -f "${primers}" ]] \
        ;
then
    printf "\n######################\nMissing required argument: PRIMERS_TABLE\n######################\n\n\n"
    usage
fi

## Is the path to the outputfolder defined?
if
    [[ -z "${output}" ]] \
        ;
then
    printf "\n######################\nMissing required argument: path to OUTPUT folder\n######################\n\n\n"
    usage
fi

## force overwrite pre-existing output folder?
if
    [[ -d "${output}" ]] &&
        [[ ${force} == "yes" ]] \
        ;
then
    rm -rf ${output}
    rm -rf ${output}

fi

if
    [[ -d "${output}" ]] \
        ;
then
    printf "\n######################\nOutput directory already exists, use -f to override\n######################\n\n\n"
    usage
fi

## remove watermark in phylogenetic trees if this is REALLY wanted
if
    [[ ${nowatermark} == "yes" ]] \
        ;
then
    WM="NO"
else
    WM="YES"
fi

###
parametersettings="${quality},${similarconsensus},${nreads},${sizerange},${minreads},${threads},${kthres},${force},${blast},${partition},${outgroup}"
###

## make sure the EOL character is \n and NOT \r\n in the samples and primers files ;-)
${ampdir}/envs/python_dependencies/bin/sed -i 's/\r$//' ${samples}
${ampdir}/envs/python_dependencies/bin/sed -i 's/\r$//' ${primers}

## remove spaces
${ampdir}/envs/python_dependencies/bin/sed -i 's/ //g' ${samples}
${ampdir}/envs/python_dependencies/bin/sed -i 's/ //g' ${primers}

## remove empty lines
${ampdir}/envs/python_dependencies/bin/sed -i '/^$/d' ${samples}
${ampdir}/envs/python_dependencies/bin/sed -i '/^$/d' ${primers}

# Test if FASTQ input files in samples file do exist
while IFS=$"," read -r samplename file; do

    ## skip header
    if
        [[ ${samplename} == "ID" ]] \
            ;
    then
        continue
    fi

    if
        [[ ! -f ${file} ]] \
            ;
    then

        printf "\n######################\n${file} does not exist, quitting\n######################\n\n\n"
        usage
    fi

done <${samples}

## fill array for samples/loci to exclude
declare -a EXCLUDE=()
if
    [[ ! -z "${exclude}" ]] \
        ;
then
    mapfile -t EXCLUDE < <(awk -F',' '{print $1$2}' "${exclude}")
fi

## set SampleID Search engine:

if [[ ${blast} != "no" ]]; then
    SE="BLAST"
else
    SE="BOLD"
fi

#################### Start executing if all parameters set and input files available
now=$(date +%Y-%m-%d' '%H:%M:%S)
echo "+++++++++ Program started at ${now} +++++++++"

## test similiarity among primers
echo "***** Test primer sequence similarity *****"

## make output folder structure
mkdir -p ${output}/results/summary/primers
mkdir -p ${output}/data/raw
mkdir -p ${output}/log/demulti
mkdir ${output}/data/demultiplexed
mkdir -p ${output}/shell/demult1
mkdir -p ${output}/log/summary
mkdir ${output}/data/filtered

conda activate ${ampdir}/envs/python_dependencies

python3 ${wd}/scripts/CompPrimers.py \
    -p ${primers} \
    -o ${output}/results/summary/primers/primers_dist.csv \
    >${output}/results/summary/min_edit.dist 2>${output}/log/summary/primerdist.log

conda deactivate

echo "finished"
## count number of loci in primers file:
LOCI=$(awk '!/^ID,/' ${primers} | wc -l)

## loop through all input files in raw data folder and store filtered files in new folder
echo "***** Copying files, starting filtering and demultiplexing by locus *****"

while IFS=$"," read -r samplename file; do

    ## skip header
    if
        [[ ${samplename} == "ID" ]] \
            ;
    then
        continue
    fi

    echo """
    set -e
    
    eval \"\$(conda shell.bash hook)\"

    ## copy raw data to new output folder

    if [[ ${file} == "\*.gz$" ]]; then
        cp -n ${file} ${output}/data/raw/${samplename}.fastq.gz
    else
        conda activate ${ampdir}/envs/chopper
        pigz -c ${file} > ${output}/data/raw/${samplename}.fastq.gz
        conda deactivate
    fi

    ## filter raw sequences with chopper
    conda activate ${ampdir}/envs/chopper
    pigz -dc ${output}/data/raw/${samplename}.fastq.gz |
        chopper -q ${quality} 2>> ${output}/log/demulti/${samplename}_demulti.log | pigz \
        >${output}/data/filtered/${samplename}-filt.fastq.gz
        
    conda deactivate

    ## demultiplex per locus and sample
    mkdir ${output}/data/demultiplexed/${samplename}

    conda activate ${ampdir}/envs/python_dependencies
    python3 ${wd}/scripts/DemultFastq.py \
        -i ${output}/data/filtered/${samplename}-filt.fastq.gz \
        -p $primers \
        -o ${output}/data/demultiplexed/${samplename} \
        -th ${kthres} \
        -sr ${sizerange} \
        -mr ${minreads} \
        -rp ${nreads} >> ${output}/log/demulti/${samplename}_demulti.log 2>&1
    conda deactivate

    #echo '${samplename} finished'

    """ >${output}/shell/demult1/${samplename}.sh

done <${samples}

## Execute shell scripts in parallel
conda activate ${ampdir}/envs/parallel
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j${threads} bash ::: ${output}/shell/demult1/*.sh
conda deactivate

## Consensus sequence reconstruction
echo "***** Files demultiplexed, starting consensus haplotype reconstruction *****"

mkdir -p ${output}/results/consensus_seqs
mkdir -p ${output}/log/ampliconsorter
mkdir -p ${output}/shell/demult2

while IFS=$"," read -r samplename file; do

    ## skip header
    if [[ ${samplename} == "ID" ]]; then
        continue
    fi

    while IFS=$"," read -r primername fwd rev size; do

        ## skip samples/loci in exclude list
        if [[ " ${EXCLUDE[*]} " == *" $samplename$primername "* ]]; then
            echo "skipping "$samplename": "$primername
            continue
        fi

        ## skip header
        if
            [[ ${primername} == "ID" ]] \
                ;
        then
            continue
        fi

        echo """
	set -e

        ## source conda
        eval \"\$(conda shell.bash hook)\"

        ## make shortcut for path and create directories
        SamPrim=${output}/results/consensus_seqs/${samplename}/${primername}
        mkdir -p ${output}/results/consensus_seqs/${samplename}
        cd ${output}/results/consensus_seqs/${samplename}

        ##AmpliconSorter consensus
        conda activate ${ampdir}/envs/python_dependencies
        ## test if OS is Mac
        if
            [[ ${os} == "Darwin" ]] \
                ;
        then
            amplicon_sorter.py \
                -i ${output}/data/demultiplexed/${samplename}/${primername}.fastq \
                -np 1 \
                -mac \
                --similar_consensus ${similarconsensus} \
                -maxr ${nreads} \
                -o \${SamPrim} \
            >> ${output}/log/ampliconsorter/${samplename}_${primername}_AS.log 2>&1
        else
            amplicon_sorter.py \
                -i ${output}/data/demultiplexed/${samplename}/${primername}.fastq \
                -np 1 \
                --similar_consensus ${similarconsensus} \
                -maxr ${nreads} \
                -o \${SamPrim} \
            >> ${output}/log/ampliconsorter/${samplename}_${primername}_AS.log 2>&1
        fi
        """ >${output}/shell/demult2/${samplename}_${primername}.sh

    done <${primers}

done <${samples}

## Execute shell scripts in parallel
conda activate ${ampdir}/envs/parallel
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j${threads} bash ::: ${output}/shell/demult2/*.sh
conda deactivate

## make CSV summary
echo "***** Summarize Ampliconsorter output *****"

${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/ParseSummary.py \
    --path ${output} \
    --primer ${primers} \
    --samples ${samples} \
    >${output}/results/summary/summary.csv

conda activate ${ampdir}/envs/python_dependencies

## choose haplotypes
${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/ChooseCons.py \
    --input ${output}/results/summary/summary.csv \
    --path ${output}/results/consensus_seqs \
    --output ${output}/results/haplotypes \
    --FreqTH ${freqthreshold} \
    >>${output}/log/ampliconsorter/Summary.log 2>&1

conda deactivate

mv ${output}/results/summary/summary.csv.ploidy ${output}/results/summary/summary.csv

## Plot Heatmap with missing data
conda activate ${ampdir}/envs/R
Rscript ${ampdir}/scripts/MissingDataHeatmap.r \
    ${output}/results/summary/summary.csv \
    >>${output}/log/ampliconsorter/Summary_AS.log 2>&1
conda deactivate

echo " finished"
## align haplotypes
echo "***** align haplotypes *****"
mkdir -p ${output}/log/tree

while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    ## skip if less than four aligned sequences
    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}.fasta ]] \
            ;
    then
        continue
    fi

    if
        [[ $(grep "^>" ${output}/results/haplotypes/${primername}/${primername}.fasta | wc -l) -lt 4 ]] \
            ;
    then
        continue
    fi

    ## mafft alignment
    conda activate ${ampdir}/envs/mafft
    mafft \
        --adjustdirectionaccurately \
        ${output}/results/haplotypes/${primername}/${primername}.fasta 2>>${output}/log/tree/${primername}_TREE.log |
        ${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/fixIDAfterMafft.py \
            --input ${output}/results/haplotypes/${primername}/${primername}.fasta \
            --Alignment - \
            >${output}/results/haplotypes/${primername}/${primername}_aln.fasta

    conda deactivate

    echo "locus ${primername} finished"

done <${primers}

echo "***** Calculate Genetic Distances *****"

conda activate ${ampdir}/envs/R

${ampdir}/envs/R/bin/Rscript ${ampdir}/scripts/GeneticDist.r \
	${output}/results/haplotypes \
    >>${output}/log/ampliconsorter/GeneticDistance_AS.log 2>&1

echo "finished"

conda deactivate

## Species identification with BOLD
mkdir ${output}/log/SpecID
PRINT=0

while IFS=$"," read -r primername fwd rev size; do

    if
        [[ ${primername} == "COX1" ]] ||
            [[ ${primername} == "ITS" ]] ||
            [[ ${primername} == "MATK_RBCL" ]] \
            ;
    then

        ## Print Header
        if
            [[ ${PRINT} == 0 ]] \
                ;
        then
            echo "***** Species ID from ${SE} *****"
            PRINT=1
        fi

        ## use BOLD/BLAST API for Species identification
        conda activate ${ampdir}/envs/python_dependencies

        if [[ ${blast} != "no" ]]; then
            mkdir -p ${output}/results/SpeciesID/${SE}/${primername}/summarized_outputs
            ${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/BLASTapi.py \
                -i ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
                -e ${blast} \
                -o ${output}/results/SpeciesID/${SE}/${primername}/summarized_outputs \
                >>${output}/log/SpecID/${primername}_SI.log 2>&1
        else
            mkdir -p ${output}/results/SpeciesID/${SE}/${primername}
            ${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/bold_api/BOLDapi.py \
                -i ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
                -p ${primername} \
                -c ${wd}/scripts/style.css \
                -n 10 \
                -o ${output}/results/SpeciesID/${SE}/${primername} \
                >>${output}/log/SpecID/${primername}_SI.log 2>&1
        fi
        conda deactivate

        echo ${primername} "finished"

    fi

done <${primers}

## test if IDs found for samples for any marker, perferably COX1
for locus in COX1 ITS MATK_RBCL; do
    if
        [[ -f ${output}/results/SpeciesID/${SE}/${locus}/summarized_outputs/final.csv ]] \
            ;
    then
        ID=${locus}
        break
    fi
done

if [[ ${LOCI} -gt 1 && $(ls -l ${output}/results/haplotypes/*/*_aln.fasta | wc -l) -gt 1 && ${freqthreshold} != "0" && ${freqthreshold} != "0.0" && ${freqthreshold} != "0.00" ]]; then

    ## Species delineation with ASAP for all concatenated Haplotypes
    echo "***** concatenate all loci *****"

    conda activate ${ampdir}/envs/asap
    mkdir -p ${output}/results/haplotypes/Concatenated_loci

    ${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/MergeAln.py \
        --input "${output}/results/haplotypes/*/*_aln.fasta" \
        --output ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci

    conda deactivate
    echo "finished"
fi

## reconstruct trees if > 3 haplotypes in input
echo "***** reconstruct ML trees per locus *****"
echo "WORKING DIR: $wd"

mkdir -p ${output}/log/tree

## adjust WIDTH to account for longer name if frequencies are also printed at freqthreshold == 0
if [[ ${freqthreshold} == 0 ]]; then
    WIDTH=10
else
    WIDTH=8
fi

while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}.fasta ]] \
            ;
    then
        continue
    fi

    ## skip if less than four aligned sequences
    if
        [[ $(grep "^>" ${output}/results/haplotypes/${primername}/${primername}.fasta | wc -l) -lt 4 ]] \
            ;
    then
        continue
    fi

    ## Phylogeny using IQtree with 100 bootsrapping rounds
    mkdir -p ${output}/results/tree/${primername}/
    cp ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
        ${output}/results/tree/${primername}/${primername}

    conda activate ${ampdir}/envs/iqtree

    iqtree \
        -s ${output}/results/tree/${primername}/${primername} \
        -ntmax ${threads} \
        -B 1000 \
        >>${output}/log/tree/${primername}_TREE.log 2>&1

    conda deactivate

    conda activate ${ampdir}/envs/R

    ## adjust tree height based on samples in dataset
    if
        [[ $(($(grep "^>" ${output}/results/haplotypes/${primername}/${primername}_aln.fasta | wc -l) / 3)) -gt 8 ]] \
            ;
    then
        HEIGHT=$(($(grep "^>" ${output}/results/haplotypes/${primername}/${primername}_aln.fasta | wc -l) / 3))
    else
        HEIGHT=8
    fi

    ## append BOLD/BLAST names if available and adjust x-axis offset to account for longer names
    OFFSET=0.3
    outgroupNew="no"
    if
        [[ ! -z ${ID} ]] \
            ;
    then
        OFFSET=0.7
        outgroupNew=$(${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/RenameTreeLeaves.py \
            --primername ${primername} \
            --input ${output}/results/tree/${primername}/${primername}.treefile \
            --name ${output}/results/SpeciesID/${SE}/${ID}/summarized_outputs/final.csv \
            --outgroup ${outgroup})
    fi

    ## plot trees with ggplot
    ${ampdir}/envs/R/bin/Rscript ${ampdir}/scripts/PlotTree.r \
        ${output}/results/tree/${primername}/${primername}.treefile \
        ${output}/results/tree/${primername}/${primername} \
        ${primername} \
        ${OFFSET} \
        ${WIDTH} \
        ${HEIGHT} \
        ${outgroupNew} \
        $WM \
        >>${output}/log/tree/${primername}_TREE.log 2>&1

    conda deactivate

    echo "ML tree for locus ${primername} finished"

done <${primers}

## repeat for combined FASTA data
if
    [[ -s "${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta" ]] \
        ;
then

    if
        [[ ${partition} == "yes" ]] \
            ;
    then
        ## Phylogeny using IQtree with 100 bootsrapping rounds
        mkdir -p ${output}/results/tree/Concatenated_loci/
        cp ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta \
            ${output}/results/tree/Concatenated_loci/Concatenated_loci.fasta
        cp ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.part \
            ${output}/results/tree/Concatenated_loci/Concatenated_loci

        conda activate ${ampdir}/envs/iqtree

        cd ${output}/results/haplotypes/Concatenated_loci

        iqtree \
            -s ${output}/results/tree/Concatenated_loci/Concatenated_loci.fasta \
            -ntmax ${threads} \
            -B 1000 \
            -p ${output}/results/tree/Concatenated_loci/Concatenated_loci \
            >>${output}/log/tree/Concatenated_loci_TREE.log 2>&1

        conda deactivate

    else

        ## Phylogeny using IQtree with 100 bootsrapping rounds
        mkdir -p ${output}/results/tree/Concatenated_loci/
        cp ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta \
            ${output}/results/tree/Concatenated_loci/Concatenated_loci

        conda activate ${ampdir}/envs/iqtree

        cd ${output}/results/haplotypes/Concatenated_loci

        iqtree \
            -s ${output}/results/tree/Concatenated_loci/Concatenated_loci \
            -ntmax ${threads} \
            -B 1000 \
            >>${output}/log/tree/Concatenated_loci_TREE.log 2>&1

        conda deactivate

    fi

    conda activate ${ampdir}/envs/R

    ## adjust tree height based on samples in dataset
    if
        [[ $(($(grep "^>" ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta | wc -l) / 3)) -gt 8 ]] \
            ;
    then
        HEIGHT=$(($(grep "^>" ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta | wc -l) / 3))
    else
        HEIGHT=8
    fi

    ## append BOLD names if available and adjust x-axis offset to account for longer names
    OFFSET=0.3
    outgroupNew="no"
    if
        [[ ! -z ${ID} ]] \
            ;
    then
        OFFSET=0.7
        outgroupNew=$(${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/RenameTreeLeaves.py \
            --primername combined \
            --input ${output}/results/tree/Concatenated_loci/Concatenated_loci.treefile \
            --name ${output}/results/SpeciesID/${SE}/${ID}/summarized_outputs/final.csv \
            --outgroup ${outgroup})

    fi

    echo ${outgroup} ${outgroupNew}
    ## plot trees with ggplot
    ${ampdir}/envs/R/bin/Rscript ${ampdir}/scripts/PlotTree.r \
        ${output}/results/tree/Concatenated_loci/Concatenated_loci.treefile \
        ${output}/results/tree/Concatenated_loci/Concatenated_loci \
        Concatenated_loci \
        ${OFFSET} \
        ${WIDTH} \
        ${HEIGHT} \
        ${outgroupNew} \
        ${WM} \
        >>${output}/log/tree/Concatenated_loci_TREE.log 2>&1

    conda deactivate

    echo "ML tree for Concatenated_loci dataset finished"

fi

if [[ ${LOCI} -gt 1 && ${freqthreshold} != "0" && ${freqthreshold} != "0.0" && ${freqthreshold} != "0.00" ]]; then
    ## now do ASTRAL concatenated trees
    mkdir -p ${output}/results/tree/ASTRAL

    ## concatenate trees
    cat ${output}/results/tree/*/*.treefile >${output}/results/astral_input.tree
    mv ${output}/results/astral_input.tree ${output}/results/tree/ASTRAL

    ## test if ANY trees available in file
    if [[ ! -z $(grep '[^[:space:]]' ${output}/results/tree/ASTRAL/astral_input.tree) ]]; then

        echo "***** reconstruct ASTRAL tree across all loci *****"

        if
            [[ $(($(cat ${samples} | wc -l) / 3)) -gt 8 ]] \
                ;
        then
            HEIGHT= $(($(cat ${samples} | wc -l) / 3))
        else
            HEIGHT=8
        fi

        if
            [[ ! -z ${ID} ]] \
                ;
        then
            OFFSET=0.5
        else
            OFFSET=0.1
        fi

        ## reconstruct astral consensus tree
        conda activate ${ampdir}/envs/aster

        wastral \
            -i ${output}/results/tree/ASTRAL/astral_input.tree \
            -o ${output}/results/tree/ASTRAL/ASTRAL.tree \
            >>${output}/log/tree/Astral.log 2>&1

        conda deactivate

        conda activate ${ampdir}/envs/R

        ## plot tree
        ${ampdir}/envs/R/bin/Rscript ${ampdir}/scripts/PlotTree_astral.r \
            ${output}/results/tree/ASTRAL/ASTRAL.tree \
            ${output}/results/tree/ASTRAL/ASTRAL \
            ASTRAL \
            ${OFFSET} \
            ${WIDTH} \
            ${HEIGHT} \
            ${outgroupNew} \
            ${WM} \
            >>${output}/log/tree/Astral.log 2>&1

        conda deactivate

        echo "finished"

    fi
fi
## Species delineation with ASAP
mkdir ${output}/log/SpecDelim

echo "***** Species delimitation with ASAP for each locus *****"
echo "WORKING DIR: $wd"
conda activate ${ampdir}/envs/asap

while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}_aln.fasta ]] \
            ;
    then
        continue
    fi

    mkdir -p ${output}/results/SpeciesDelim/${primername}
    cd ${output}/results/SpeciesDelim/${primername}

    asap \
        -a ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
        -o ${output}/results/SpeciesDelim/${primername} \
        >>${output}/log/SpecDelim/${primername}_SD.log 2>&1

    # fi

    echo "locus ${primername} finished"

done <${wd}/${primers}

## run ASAP on concatenated FASTA if file not empty
if
    [[ -s "${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta" ]] \
        ;
then

    mkdir ${output}/results/SpeciesDelim/Concatenated_loci

    asap \
        -a ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta \
        -o ${output}/results/SpeciesDelim/Concatenated_loci \
        >>${output}/log/SpecDelim/Concatenated_loci_SD.log 2>&1
    #fi
    conda deactivate

fi
conda deactivate

cd $wd
## make HTML summary

mkdir -p ${output}/Output

echo "***** Summarize in HTML file *****"

conda activate ${ampdir}/envs/python_dependencies

mkdir ${output}/results/html
mkdir ${output}/log/html

## visualize MSA
unset SESSION_MANAGER

${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/msa_to_html.py \
    -hap ${output}/results/haplotypes/ \
    -py "${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/msaviz.py" \
    -hf ${output}/results/html/MSA_analysis.html \
    -sf ${wd}/scripts/scripts/style.css \
    >>${output}/log/html/MSA.log 2>&1

unset SESSION_MANAGER

samplestba=$(awk -F',' '{print$1}' "${samples}")
primerstba=$(awk -F',' 'NR>1 {print $1}' "${primers}")

## make complete HTML output

mkdir ${output}/.logos
cp ${wd}/imgs/tettris.png ${output}/.logos/tettris.png
cp ${wd}/imgs/nhm.svg.png ${output}/.logos/nhm.svg.png

${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/altersvg.py \
    -f ${output}/results/SpeciesDelim/ \
    >>${output}/log/html/HTML.log 2>&1

${ampdir}/envs/python_dependencies/bin/python3 ${ampdir}/scripts/DisplayOutput.py \
    -p ${parametersettings} \
    -r ${output}/results \
    -out ${output}/Output \
    -loci ${primerstba} \
    -samples ${samplestba} \
    >>${output}/log/html/HTML.log 2>&1

conda deactivate

if [[ -d ${output}/results/SpeciesID ]]; then
    for PA in ${output}/results/SpeciesID/${SE}/*; do
        IDlocus=${PA##*/}
        cp ${PA}/summarized_outputs/final.csv ${output}/Output/summary/SpeciesID_${IDlocus}.csv
    done
fi

## copy haplotypes to Output folder
while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}_aln.fasta ]] \
            ;
    then
        continue
    fi

    cp ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
        ${output}/Output/haplotypes/${primername}

done <${primers}

rm -f ${output}/results/haplotypes/Rplots.pdf
rm -f ${output}/Output/haplotypes/Rplots.pdf

## Finished
now=$(date +%Y-%m-%d' '%H:%M:%S)
echo " +++++++++ Program finished at ${now} +++++++++ "
