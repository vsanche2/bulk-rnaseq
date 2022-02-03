#!/bin/bash
conda activate rnaseq
~/get_a_worker_node.sh ##Don't forget to get a worker node first!!! Do not run this on the main node

IN_SS=SampleSheet.txt
OUT_DIR=RUN2
OUT_QC=RUN2/QC
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
fastp -h ${OUT_QC}/${SAMPLE}.${COND}.fastp.html -j ${OUT_QC}/${SAMPLE}.${COND}.fastp.json \
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

#Run fastQC https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help
fastqc ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz --outdir $OUT_QC
fastqc ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz --outdir $OUT_QC

#Run STAR
mkdir -p ${OUT_DIR}/STAR/${SAMPLE}.${COND}

STAR --runMode alignReads \
--genomeDir $STAR_HPV_HYBRID_IDX \
--sjdbGTFfile $gtfreference $GTF \
--readFilesCommand zcat \
--runThreadN 6 \
--quantMode TranscriptomeSAM GeneCounts \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMultimapNmax 1 \
--outFileNamePrefix ${SAMPLE}.${COND}. \
--readFilesIn ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${OUT_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz

#--outFilterMatchNmin 35 \



#Alignment QC http://qualimap.conesalab.org/doc_html/command_line.html
qualimap rnaseq -outdir $OUT_QC \
-a proportional \
-bam ${OUT_DIR}/STAR/${SAMPLE}.${COND}/${SAMPLE}.${COND}.Aligned.sortedByCoord.out.bam \
-gtf $GTF \
--java-mem-size=8G
#-p strand-specific-reverse \

#RNASEQC https://github.com/getzlab/rnaseqc
rnaseqc --coverage $GTF ${OUT_DIR}/STAR/${SAMPLE}.${COND}/${SAMPLE}.${COND}.Aligned.sortedByCoord.out.bam $OUT_QC

#Run MultiQC
multiqc -f $OUT_QC ${OUT_DIR}/STAR/

done < <(cat $IN_SS)
date



