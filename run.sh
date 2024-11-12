# define the path to your test folder
WD='/app/test'

# generate samples.csv file

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

# run AmpliPiper
bash /app/shell/AmpliPiper.sh \
    -s /app/test/testdata/data/samples.csv \
    -p /app/test/testdata/data/primers.csv \
    -o /app/test/testdata/results/demo \
    --quality 10 \
    --nreads 1000 \
    --blast your@email.com \
    --similar_consensus 97 \
    --threads 200 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41 \
    --force