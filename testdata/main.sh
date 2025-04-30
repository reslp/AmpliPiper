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

# (3a) run AmpliPiper

bash ${WD}/shell/AmpliPiper.sh \
    --samples ${WD}/testdata/data/samples.csv \
    --primers ${WD}/testdata/data/primers.csv \
    --output ${WD}/testdata/results/demo \
    --quality 10 \
    --nreads 1000 \
    --blast your@email.com \
    --similar_consensus 97 \
    --threads 10 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41 \
    --force

# (3b) run AmpliPiper and keep all consensus sequences

bash ${WD}/shell/AmpliPiper.sh \
    --samples ${WD}/testdata/data/samples.csv \
    --primers ${WD}/testdata/data/primers.csv \
    --output ${WD}/testdata/results/demo_allCons \
    --quality 10 \
    --nreads 1000 \
    --blast your@email.com \
    --similar_consensus 97 \
    --threads 100 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41 \
    --freqthreshold 0 \
    --force
