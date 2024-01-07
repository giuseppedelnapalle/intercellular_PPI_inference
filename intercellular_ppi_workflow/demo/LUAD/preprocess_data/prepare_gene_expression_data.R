#!/usr/bin/env Rscript
# $title preparte gene expression data for CIBERSORTx
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

# directory of input
dir_input <- "/home/nikola/Project_Data/R_data/cancer_research/TCGA_datasets"
dir_input <- paste(dir_input, dtset, "input/", sep = "/")

# directory of imported functions
dir_func <- "/home/nikola/Project_Data/R_data/functions_scripts/gene_id_conversion/"
dir_func2 <- "/home/nikola/Project_Data/R_data/functions_scripts/data_cleaning/"

# directory of genes of PPI_db
dir_gene <- "/home/nikola/Project_Data/R_data/PPI_db/merge_db/output/"
dir_gene2 <- "/home/nikola/Project_Data/R_data/PPI_db/LRdb_prep/"

# # directory of output
dir_output <- "./output/"
dir.create(dir_output)


# 1 load data -------------------------------------------------------------

load(file = paste0(dir_input, "rna.hg38.RData"))
load(file = paste0(dir_input, "mat_TPM.RData"))
load(file = paste0(dir_gene, "g_comb.RData"))
load(file = paste0(dir_gene2, "genes_LRdb.rda"))

# 2 convert ensembl_ID to gene_symbol -------------------------------------

source(file = paste0(dir_func, "ensembl_to_geneSymbol.R"))
mat_TPM_gs <- ensembl_to_geneSym(mat_TPM)
dim(mat_TPM)
dim(mat_TPM_gs)
sum(duplicated(rownames(mat_TPM_gs)))

# 3 remove lowly-expressed genes ------------------------------------------

# remove genes if TPM < thresh_tpm in over thresh_prop of samples
source(file = paste0(dir_func2, "filter_genes.R"))

summary(apply(mat_TPM_gs, 1, median))
summary(apply(mat_TPM_gs, 2, median))
ncol(mat_TPM_gs)

thresh_tpm <- 1
thresh_prop <- .98

mat_TPM_filt <- filter_genes(mat_TPM_gs, thresh_tpm, thresh_prop)

# test --------------------------------------------------------------------

tmp <- filter_genes(mat_TPM_gs, thresh_tpm, thresh_prop)

dim(tmp)
sum(rownames(tmp) %in% g_comb)
length(g_comb)

summary(apply(tmp, 1, median))
summary(apply(tmp, 2, median))
rm(tmp)

# 4 remove outlying samples ------------------------------------------------

# phenotype data
phenoDT <- colData(rna.hg38)
rm(rna.hg38)
# verify order of sample matched
sum(match(colnames(mat_TPM_filt),rownames(phenoDT))==c(1:ncol(mat_TPM_filt)))==ncol(mat_TPM_filt)

# convert TPM to log2 space
log2_TPM <- log2(mat_TPM_filt+1)

# 4.1 PCA -----------------------------------------------------------------

source(file = paste0(dir_func2, "PCA.R"))

# centered ----------------------------------------------------------------

# transform mean of genes to 0
y <- log2_TPM - rowMeans(log2_TPM)
summary(apply(y, 1, mean))

# labeled by pathology
PCA_scatterplot(t(y), meta = phenoDT, batch = phenoDT[,1:3], color = "primary_diagnosis",
         scale = F, height = 10, width = 20, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "pathology_centered")

# labeled by stage
PCA_scatterplot(t(y), meta = phenoDT, batch = phenoDT[,1:3], color = "tumor_stage",
         scale = F, height = 10, width = 16, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "stage_centered")

# labeled by sample_type
PCA_scatterplot(t(y), meta = phenoDT, batch = phenoDT[,1:3], color = "sample_type",
         scale = F, height = 10, width = 16, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "sample_type_centered")

# scaled ------------------------------------------------------------------

# labeled by pathology
PCA_scatterplot(t(log2_TPM), meta = phenoDT, batch = phenoDT[,1:3], color = "primary_diagnosis",
         scale = T, height = 10, width = 20, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "pathology_scaled")

# labeled by stage
PCA_scatterplot(t(log2_TPM), meta = phenoDT, batch = phenoDT[,1:3], color = "tumor_stage",
         scale = T, height = 10, width = 16, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "stage_scaled")

# labeled by sample_type
PCA_scatterplot(t(log2_TPM), meta = phenoDT, batch = phenoDT[,1:3], color = "sample_type",
         scale = T, height = 10, width = 16, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "sample_type_scaled")

# 4.2 outlier detection & removal -----------------------------------------

# reference WGCNA tutorial

# constructing a sample network for outlier detection
source(file = paste0(dir_func2, "remove_outliers.R"))

#table(phenoDT$primary_diagnosis)

log2_TPM_filt <- remove_outliers(log2_TPM, phenoDT, trait = c("sample_type", "primary_diagnosis", "tumor_stage"),
                                thresh_z = -2.5, height = 16, width = 20,
                                directory = paste0(dir_output, "plots"), 
                                file_suffix = "log2TPM_thresh_m2.5")

# labeled by pathology
PCA_scatterplot(t(log2_TPM_filt), meta = phenoDT[rownames(phenoDT) %in% colnames(log2_TPM_filt),], 
         batch = phenoDT[rownames(phenoDT) %in% colnames(log2_TPM_filt),1:3], color = "primary_diagnosis",
         scale = T, height = 10, width = 20, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "pathology_scaled_obs_filt_z2.5")

# labeled by sample_type
PCA_scatterplot(t(log2_TPM_filt), meta = phenoDT[rownames(phenoDT) %in% colnames(log2_TPM_filt),], 
         batch = phenoDT[rownames(phenoDT) %in% colnames(log2_TPM_filt),1:3], color = "sample_type",
         scale = T, height = 10, width = 16, point.size = 1,
         directory = paste0(dir_output, "plots"), file_suffix = "sample_type_scaled_obs_filt_z2.5")


# 5 export mixture data for CIBERSORTx ------------------------------------

# transform to TPM
mat_TPM_filt2 <- 2^log2_TPM_filt-1

save(mat_TPM_filt2, file = paste0(dir_output, "mat_TPM_filt2.RData"))

# keep only PPI_db genes to reduce size of file
mat_TPM_filt2_sub <- mat_TPM_filt2[rownames(mat_TPM_filt2) %in% g_comb,]
mat_TPM_filt2_sub <- cbind(as.data.frame(rownames(mat_TPM_filt2_sub)), as.data.frame(mat_TPM_filt2_sub))
colnames(mat_TPM_filt2_sub)[1] <- "gene"

write.table(mat_TPM_filt2_sub, file = paste0(dir_output, "mat_TPM_filt_sub_", 
                                         sapply(strsplit(dtset,"_"), "[", 1), ".txt"),
            quote = F, sep = "\t", row.names = F, col.names = T)

# keep LRdb genes
mat_TPM_filt2_LRdb <- mat_TPM_filt2[rownames(mat_TPM_filt2) %in% genes_LRdb,]
mat_TPM_filt2_LRdb <- cbind(as.data.frame(rownames(mat_TPM_filt2_LRdb)), as.data.frame(mat_TPM_filt2_LRdb))
colnames(mat_TPM_filt2_LRdb)[1] <- "gene"

write.table(mat_TPM_filt2_LRdb, file = paste0(dir_output, "mat_TPM_filt_LRdb_", 
                                             sapply(strsplit(dtset,"_"), "[", 1), ".txt"),
            quote = F, sep = "\t", row.names = F, col.names = T)

# corresponding phenotype data
phenoDT_filt2 <- phenoDT[rownames(phenoDT) %in% colnames(mat_TPM_filt2),]
sum(match(rownames(phenoDT_filt2), colnames(mat_TPM_filt2))==c(1:nrow(phenoDT_filt2)))==nrow(phenoDT_filt2)
# View(as.data.frame(phenoDT_filt2))
# table(as.data.frame(phenoDT_filt2)[,"sample_type"])

save(phenoDT_filt2, file = paste0(dir_output, "phenoDT_filt2.RData"))

# inspect removed samples
# phenoDT[!rownames(phenoDT) %in% colnames(mat_TPM_filt2),"sample_type"]

# sample_type
stype_filt2 <- phenoDT_filt2$sample_type
stype_filt2 <- gsub(" ", "_", stype_filt2)
stype_filt2 <- ifelse(grepl("Tumor", stype_filt2, fixed = T), "tumour", "control")
write.table(stype_filt2, file = paste0(dir_output, "sample_type_filt2_", 
                                             sapply(strsplit(dtset,"_"), "[", 1), ".txt"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# save image --------------------------------------------------------------

rm(log2_TPM, log2_TPM_filt, mat_TPM, mat_TPM_filt, mat_TPM_filt2, mat_TPM_filt2_sub, mat_TPM_gs, y, mat_TPM_filt2_LRdb)

save.image("~/Project_Data/R_data/cancer_research/intercellular_ppi_ca/proof_concept/LUAD/prepare_gene_expression_data.RData")
