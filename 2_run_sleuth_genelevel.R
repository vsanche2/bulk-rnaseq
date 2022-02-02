
sample_id = dir(file.path("~/projects/bulk_rna_seq/RUN1", "QUANT"))
kal_dirs = file.path("~/projects/bulk_rna_seq/RUN1", "QUANT", sample_id)

sample_sheet = "~/projects/bulk_rna_seq/SampleSheet.txt"
so_object_name = "~/projects/bulk_rna_seq/RUN1/slth_obj_norm_TPMtrans"
tpm_file = "~/projects/bulk_rna_seq/RUN1/abundance_genelevel_TPMtrans.txt"
### Edit above block to appropriate paths/filenames
###-----------------------------------------------------------------------------------------

sessionInfo()
library(sleuth)
library(dplyr)

sample_id
kal_dirs

s2c <- read.table(sample_sheet, header = FALSE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = V1, condition = V2)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c

library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
head(t2g,100)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = t2g, read_bootstrap_tpm=TRUE,
                  aggregation_column="ext_gene", gene_mode=TRUE, normalize=TRUE)
#,transformation_function = function(x) log2(x+0.5))
#transform_fun_tpm = function(x) log2(x+0.5))

sleuth_save(so,so_object_name)

pdf("PCA.pdf")
plot_pca(so, color_by = 'sample')
dev.off()

pdf("Densities.pdf")
plot_group_density(so, use_filtered = TRUE, units = "tpm",trans = "log2", grouping = "sample", offset = 1)  
dev.off()

#"scaled_reads_per_base" for gene-level analyses #scaled_reads_per_base
sleuth_matrix_genelevel <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
write.table(sleuth_matrix_genelevel, tpm_file, row.names = T, sep = "\t", quote = F, append = F)

## Need replicates for finding within condition variance, so the below will FAIL.
#so <- sleuth_fit(so, ~condition, 'full')
#so <- sleuth_fit(so, ~1, 'reduced')
#so <- sleuth_lrt(so, 'reduced', 'full')
#sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
#head(sleuth_significant, 20)

