#!/bin/bash
conda activate rnaseq
~/get_a_worker_node.sh ##Don't forget to get a worker node first!!! Do not run this on the main node

IN_SS=SampleSheet.txt
OUT_DIR=RUN2
IN_DIR=NEO-fastqs
TRIM_FASTQ_DIR=RUN1

#For Kallisto
KALLISTO_IDX=~/projects/bulk_rna_seq/kallisto/genome/homo_sapiens/transcriptome.idx

##Edit block above to point to appropriate locations
##-----------------------------------------------------------------------------------------------------------------------

echo "Input dir: "$IN_DIR
echo "Output dir: "$OUT_DIR
date

while IFS=$'\t' read line
do
	set $line
    SAMPLE=$1
    COND=$2
    echo "Started working on: "$SAMPLE $COND

#Merge Lanes
cat $IN_DIR/$SAMPLE*R1* > ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz
cat $IN_DIR/$SAMPLE*R2* > ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz

#Run fastp
fastp -h ${OUT_DIR}/${SAMPLE}.${COND}.fastp.html -j ${OUT_DIR}/${SAMPLE}.${COND}.fastp.json \
-i ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz -I ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz \
-o ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz -O ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz \
--dedup \
--cut_front \
--cut_front_window_size \
--cut_tail_window_size \
--cut_tail \
--n_base_limit \
--length_required 40 \
--low_complexity_filter \
--overrepresentation_analysis
#-y -3 --cut_tail_window_size 4 -5 --cut_front_window_size 4 --length_required 40

#Run fastp
fastqc $OUT_DIR

#Run Kallisto
mkdir -p ${OUT_DIR}/QUANT/${SAMPLE}.${COND}
kallisto quant -i $KALLISTO_IDX -o ${OUT_DIR}/QUANT/${SAMPLE}.${COND} --pseudobam --bias -b 100 -t 8 --fusion --rf-stranded ${TRIM_FASTQ_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${TRIM_FASTQ_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz &> ${OUT_DIR}/QUANT/${SAMPLE}.${COND}/${SAMPLE}.${COND}.kallisto.log

#Run MultiQC
multiqc -f $OUT_DIR ${OUT_DIR}/QUANT/

done < <(cat $IN_SS)
date
