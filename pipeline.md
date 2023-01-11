## pipeline

```shell
####### merge mRNA fastq files
cat mLpo2_L2_1001453.R1.fastq.gz mLpo2_L4_1001453.R1.fastq.gz >mLpo2_merge_1001453.R1.fastq.gz
cat mLpo2_L2_1001453.R2.fastq.gz mLpo2_L4_1001453.R2.fastq.gz >mLpo2_merge_1001453.R2.fastq.gz
cat mLpo5_L2_1002905.R1.fastq.gz mLpo5_L4_1002905.R1.fastq.gz >mLpo5_merge_1002905.R1.fastq.gz
cat mLpo5_L2_1002905.R2.fastq.gz mLpo5_L4_1002905.R2.fastq.gz >mLpo5_merge_1002905.R2.fastq.gz

####### fastqc
qcdir=~/project/3.qc/raw_qc/mRNA/merge_mRNA
fqdir=~/project/1.raw_fq/mRNA/merge_mRNA
fastqc -t 16 -o $qcdir $fqdir/*.fastq.gz

####### trim_galore
rawdata=~/project/1.raw_fq/mRNA/merge_mRNA
cleandata=~/project/2.clean_fq/mRNA/merge_mRNA
cat ~/project/1.raw_fq/mRNA/merge_mRNA/config | while read id
do
        echo "trim_galore  --cores 32 --phred33 -q 20 --length 36 --stringency 3 --fastqc --paired --max_n 3 -o ${cleandata} ${rawdata}/${id}.R1.fastq.gz ${rawdata}/${id}.R2.fastq.gz"
done >trim_galore.sh

sh trim_galore.sh

####### star
index=~/database/index/star_index
inputdir=~/project/2.clean_fq/mRNA/merge_mRNA
outdir==~/project/4.mapping/mRNA/merge_mRNA
cat ~/project/1.raw_fq/mRNA/merge_mRNA/config | while read id
do
        echo "STAR --runThreadN 32 \
            --readFilesCommand zcat \
            --genomeDir ${index}  \
            --readFilesIn ${inputdir}/${id}.R1_val_1.fq.gz  ${inputdir}/${id}.R2_val_2.fq.gz  \
            --outFileNamePrefix ${id} \
            --outSAMtype BAM Unsorted >${id}.log"
done >STAR.sh
sh STAR.sh >STAR.log

####### samtools
outdir=~/project/4.mapping/mRNA/merge_mRNA
cat ~/project/1.raw_fq/mRNA/merge_mRNA/config | while read id
do
        echo "samtools view -@ 16 -hb -q 255 ${outdir}/${id}Aligned.out.bam | \
        samtools sort -@ 16 -o ${outdir}/${id}_uniq_sort.bam
        samtools index -@ 16 ${outdir}/${id}_uniq_sort.bam   ${outdir}/${id}_uniq_sort.bam.bai"
done >samtools_sort.sh
sh samtools_sort.sh >samtools_sort.log

####### featurecount
gtf=~/database/genome/Human/Ensembl/GRCh38_release102/Homo_sapiens.GRCh38.102.chr.gtf
inputdir=~/project/4.mapping/mRNA/merge_mRNA
outdir=~/project/5.Expression/featureCount/mRNA/count
featureCounts -T 16 -a $gtf -p -t exon  -g gene_id  -o sle.all.id.txt $inputdir/*_uniq_sort.bam
```

