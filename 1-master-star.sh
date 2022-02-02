#!/bin/bash
conda activate rnaseq
~/get_a_worker_node.sh ##Don't forget to get a worker node first!!! Do not run this on the main node

IN_SS=SampleSheet.txt
OUT_DIR=RUN2
IN_DIR=NEO-fastqs
TRIM_FASTQ_DIR=RUN1

#For STAR
GTF="~/projects/bulk_rna_seq/HPVDetection/HPV-Hybrid/human_hpv.gtf"
STAR_HPV_HYBRID_IDX="~/projects/bulk_rna_seq/HPVDetection/HPV-Hybrid/index"

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

#Run STAR
mkdir -p ${OUT_DIR}/STAR/${SAMPLE}.${COND}

STAR --genomeDir $STAR_HPV_HYBRID_IDX \
--readFilesCommand zcat \
--runThreadN 6 \
--readFilesIn ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFileNamePrefix ${SAMPLE}.${COND}

#Alignment QC
qualimap rnaseq -outdir $OUT_DIR \
-a proportional \
-bam ${OUT_DIR}/STAR/${SAMPLE}.${COND}.bam \
-p strand-specific-reverse \
-gtf $GTF \
--java-mem-size=8G

#Run MultiQC
multiqc -f $OUT_DIR ${OUT_DIR}/STAR/

done < <(cat $IN_SS)
date



