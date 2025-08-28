#！bin/bash
methy_dir=/mnt/san3/usr/cjx/neurospora_crassa/cxn_temp/bismark_results
genome_dir=/home/cjx/Data/data/neurospora/N2489-genome
fq_dir=/mnt/san3/usr/cjx/neurospora_crassa/cxn_temp/fq

# cd /mnt/san3/usr/cjx/neurospora_crassa/Eugenecsr/genome/changed_genome
# /mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark_genome_preparation .

for file in `find ${fq_dir} -name "*1.fq.gz" | sed 's/_1.fq.gz$//'`
do
    sample=`basename ${file}`

    cd ${methy_dir} &&  mkdir ${sample}

    /mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark \
    --genome ${genome_dir} \
    -o ${sample} \
    -basename ${sample} \
    -bowtie2 \
    -p 8\
    -1 ${fq_dir}/${sample}_1.fq.gz\
    -2 ${fq_dir}/${sample}_2.fq.gz\
    &> ${sample}/${sample}.log

    /mnt/san1/usr/thw/tansoft/bismark_v0.22.1/deduplicate_bismark \
    --bam \
    --output_dir ${sample} ${sample}/${sample}_pe.bam \
    &>> ${sample}/${sample}.log

    /mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark_methylation_extractor \
    --gzip \
    --buffer_size 80G \
    --comprehensive \
    --bedGraph \
    --CX \
    --cytosine_report \
    --genome_folder ${genome_dir} \
    -o ${sample} ${sample}/${sample}_pe.deduplicated.bam \
    &>> ${sample}/${sample}.log

    cd ${sample} && /mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark2report

    /mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark2summary ${sample}_pe.bam \
    -o ${sample}.summary_report

    # merge window methylation from site

    less ${sample}/${sample}_pe.deduplicated.CX_report.txt.gz \
    |awk -v OFS='\t' '$4+$5>5{print $1,$2-1,$2,$4/($4+$5)}' |LC_ALL=C sort -S 80% -k 1,1 -k 2n - \
    -o ${sample}/${sample}.CX.bed

    LC_ALL=C sort -S 80% -k 1,1 -k 2n /home/cjx/Data/data/neurospora/Nc2489_200bin.bed \
    |bedtools map -a - \
    -b ${sample}/${sample}.CX.bed \
    -c 4 -o mean | awk '$4!="."{print $0}' \
    > ${sample}/${sample}.CX.200w.bed

    /home/cjx/Data/set_up_backage/UCSC_Kent/bedGraphToBigWig ${sample}/${sample}.CX.200w.bed \
    /home/cjx/Data/data/neurospora/Nc2489.chrom.size \
    ${sample}/${sample}.CX.200w.bw

    gzip ${sample}/${sample}.CX.bed

    ## Calculation area methylation

    echo -n -e "${sample}\t" && bedtools intersect -a ${sample}/${sample}.CX.200w.bed \
    -b /home/cjx/Data/data/neurospora/dup/duplication_3kinds.s.bed \
    -wo \
    |awk -v OFS='\t' '{print $8,$9,$4,$1,$2,$3}' \
    |sort \
    |awk '$2==200' \
    |groupBy -g 1 -c 3 -o mean \
    |xargs > ${sample}/${sample}.dupkind.methylation.txt

done 
