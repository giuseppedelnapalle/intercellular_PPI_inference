#!/usr/bin/env Rscript
# $title preparte gene expression data for CIBERSORTx
# $description keep genes from LRdb_man_combined PPIs
# $dataset LUAD
# $input gene_expression_matrix
# $author giuseppe
# $date Apr 2020


# 0 setting up R session --------------------------------------------------

dtset <- "LUAD_lung_adenocarcinoma"

library(dplyr)
#library(TCGAbiolinks)
library(SummarizedExperiment)

# turn off stringsAsFactors
options(stringsAsFactors = F)

# directory of filtered_TPM
dir_tpm <- "/home/nikola/Project_Data/R_data/cancer_research/intercellular_ppi_ca/proof_concept/LUAD/output/"

# directory of genes of PPI_db
dir_gene <- "/home/nikola/Project_Data/R_data/PPI_db/merge_LRdb_man/output/"

# # directory of output
dir_output <- "./output/"


# 1 load data -------------------------------------------------------------

load(file = paste0(dir_tpm, "mat_TPM_filt2.RData"))
load(file = paste0(dir_gene, "genes_LRdb_man.RData"))
length(genes)


# 2 filter gene expression data -------------------------------------------

mat_TPM_filt2_comb <- mat_TPM_filt2[rownames(mat_TPM_filt2) %in% genes,]
# dim(mat_TPM_filt2_comb)
mat_TPM_filt2_comb <- cbind(as.data.frame(rownames(mat_TPM_filt2_comb)), as.data.frame(mat_TPM_filt2_comb))
colnames(mat_TPM_filt2_comb)[1] <- "gene"

write.table(mat_TPM_filt2_comb, file = paste0(dir_output, "mat_TPM_filt_LRdb_man_", 
                                              sapply(strsplit(dtset,"_"), "[", 1), ".txt"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# save image --------------------------------------------------------------

rm(mat_TPM_filt2, mat_TPM_filt2_comb)
save.image("~/Project_Data/R_data/cancer_research/intercellular_ppi_ca/proof_concept/LUAD/filter_TPM_by_LRdb_man_ppi_genes.RData")
