#!/usr/bin/env Rscript
# $title create phenotype_class file
# $description
#   create phenotype_classes file to genereate signature matrix for CIBERSORTx
#   reference https://cibersortx.stanford.edu/tutorial.php
# $author giuseppe
# $date Apr 2020


# generate phenotype_class file for sorted RNA-seq data

# "Each data point should have a value of "0", "1" or "2" (without the double quotes). 
# A value of "1" indicates membership of the reference sample to the class as defined in that row, 
# a value of "2" indicates the class that the sample will be compared against, 
# and a value of "0" indicates that the comparison will be ignored."

pheno_class_RNA_seq <- function(ref_sample = ref_sample,  # data_frame, required by CIBERSORTx
                                skipped_sample = skipped_sample,  # data_frame, labeling if sample to be ignored in each class (col)
                                directory = directory,
                                file_suffix = file_suffix){
  obs_lst <- colnames(ref_sample)[-1]
  pheno_cls <- as.character(levels(factor(obs_lst)))
  
  pheno_lst <- list()
  for (cls in pheno_cls) {
    i <- match(cls, pheno_cls)
    pheno_lst[[i]] <- ifelse(obs_lst == cls, 1, 2)
    
    # override class if specified by skipped_sample as 0
    v_skip <- skipped_sample[,cls]
    pheno_lst[[i]] <- sapply(1:length(obs_lst), function(j){
      if (v_skip[j] == 0) pheno_lst[[i]][j] <- 0
      return(pheno_lst[[i]][j])
    })  # end of function of j
  }  # end of for loop
  
  # transform pheno_lst
  df <- as.data.frame(do.call(cbind, pheno_lst))
  colnames(df) <- pheno_cls
  df <- t(df)
  
  write.table(df, file = paste0(directory, "pheno_class_", file_suffix, ".txt"), 
              quote = F, sep = "\t", row.names = T, col.names = F)
}
