## created by cjx
## Wed May 22 16:54:43 CST 2024


## fastqc
echo "beginning...fastqc..."
mkdir fastqc-results

for i in `find . -name "*fq.gz"`;
do
    fastqc \
    --noextract \
    --format fastq \
    --threads 8 \
    --outdir ./fastqc-results \
    $i
done
echo "done"

cd fastqc-results && multiqc .

## Optional trim

# ##single
# for i in `ls *fq.gz`
# do
# trim_galore --quality 30 --length 20 $i  -gzip  -o trim  --cores 8  --fastqc 
# done

##paired
for i in `ls *_1.clean.fq.gz`
do
sample=${i%%_1.clean*}
trim_galore --cores 8 --paired --quality 30 --length 20 -o trim ${sample}_1.clean.fq.gz ${sample}_2.clean.fq.gz 
done

## alignment:/paired-end/  |  /single-end/ need to change bowtie&&sambamba's argument

# path=/mnt/san3/usr/cjx/neurospora_crassa/
# fqpath=/mnt/san3/usr/cjx/neurospora_crassa/fqclean/
# genome=/home/cjx/Data/data/neurospora/

mkdir $path/bam

for i in `find ${fqpath} -name "*_1.clean.fq.gz"`
do

    sampleid=`basename $i`
    sample=${sampleid%%_1.clean*}
    species=Nc2489
    
    #mapping and sort
    bowtie2 \
    -x ${genome}/$species \
    -p 8 \
    -1 ${fqpath}/${sample}_1.clean.fq.gz \
    -2 ${fqpath}/${sample}_2.clean.fq.gz \
    | samtools sort \
    --threads 8 \
    -O bam \
    -o $path/bam/${sample}_sort.bam \
    >> $path/bam/${sample}.log 2>&1

    #dedup
    java -jar picard.jar MarkDuplicates \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=true \
        INPUT=$path/bam/${sample}_sort.bam \
        OUTPUT=$path/bam/${sample}_sort.rmdup.bam \
        METRICS_FILE=$path/bam/${sample}_sort.rmdup.metrics \
        > $path/bam/${sample}.log 2>&1 
    
    samtools index $path/bam/${sample}_sort.rmdup.bam && rm $path/bam/${sample}_sort.bam

done

## call peak , need ip-sample's and input sample's ID begginning with "ip-"/"input-"
## alignment: "--broad" used when chip was enriched into broad peaks, like H3K27me3/H3K9me3, not like H3K4me3.
cd $path/bam
mkdir macs2_results
for i in `find $path -name "ip*bam"`
do
    ipsample=`basename $i`
    name=${ipsample#*ip-}
    control=input-$name
    out=${name%%_sort*}
    macs2 callpeak \
        -f BAM \
        -t $path/bam/$ipsample \
        -c $path/bam/$control \
        -n $out \
        -g 40996763 \
        --outdir $path/bam/macs2_results \
        --qvalue 0.05 --nomodel --SPMR --bdg --broad \
    2> $path/bam/macs2_results/${out}.log

    macs2 bdgcmp \
        -t $path/bam/macs2_results/${out}_treat_pileup.bdg \
        -c $path/bam/macs2_results/${out}_control_lambda.bdg \
        -o $path/bam/macs2_results/${out}_FE.bdg \
        -m FE

    bedGraphToBigWig \
        $path/bam/macs2_results/${out}_FE.bdg \
        /home/cjx/Data/data/neurospora/Nc2489.chrom.size \
        $path/bam/macs2_results/${out}_FE.bw
done

