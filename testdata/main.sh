## Run testdataset

# (1) Define path to AmpliPiper folder
# replace '<path to AmpliPiper folder>' with real path
WD='<path to AmpliPiper folder>'
## e.g., WD=/media/inter/pipelines/AmpliPiper

# (2) generate samples.csv file

## print header for samples.csv
printf "ID,PATH\n" >${WD}/testdata/data/samples.csv

## loop through input FASTQ files
for Filepath in ${WD}/testdata/reads/*fastq.gz; do

    ## get Filename
    Filename=${Filepath##*/}

    ## get File ID
    ID=${Filename%.fastq.gz*}

    ## print to samples.csv
    echo ${ID},${Filepath} >>${WD}/testdata/data/samples.csv
done

# (3) run AmpliPiper

bash ${WD}/shell/AmpliPiper.sh \
    --samples ${WD}/testdata/data/samples.csv \
    --primers ${WD}/testdata/data/primers.csv \
    --output ${WD}/testdata/results/demo2 \
    --quality 10 \
    --nreads 200 \
    --blast capoony@gmail.com \
    --similar_consensus 97.5 \
    --threads 200 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41 \
    --force
