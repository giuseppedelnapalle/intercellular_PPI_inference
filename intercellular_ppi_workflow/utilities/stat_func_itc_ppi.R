#!/usr/bin/env Rscript
# $title statistical functions for intercellular PPI inference projects
# $author giuseppe
# $date Apr 2020

library(coin)
library(multtest)

# required function
dir_tmp <- "/home/nikola/Project_Data/R_data/script_function/statistics/boot_t_test/"
source(file = paste0(dir_tmp, "bootstrap_t_test.R"))
rm(dir_tmp)


# summary of variance in each cell_type_comb
var_summary <- function(df_ppi = df_ppi, 
                        cell_type_combs = cell_type_combs,
                        axis = 1){
  res <- lapply(ctb_trgt, function(x){
    dt <- df_ppi[, grep(x, colnames(df_ppi), fixed = T, value = F)]
    summary(apply(dt, axis, var))
  })
  
  names(res) <- ctb_trgt
  return(res)
}

# summary of median in each cell_type_comb
median_summary <- function(df_ppi = df_ppi, 
                           cell_type_combs = cell_type_combs,
                           axis = 1){
  res <- lapply(ctb_trgt, function(x){
    dt <- df_ppi[, grep(x, colnames(df_ppi), fixed = T, value = F)]
    summary(apply(dt, axis, median))
  })
  
  names(res) <- ctb_trgt
  return(res)
}

# summary of mean in each cell_type_comb
mean_summary <- function(df_ppi = df_ppi, 
                           cell_type_combs = cell_type_combs,
                           axis = 1){
  res <- lapply(ctb_trgt, function(x){
    dt <- df_ppi[, grep(x, colnames(df_ppi), fixed = T, value = F)]
    summary(apply(dt, axis, mean))
  })
  
  names(res) <- ctb_trgt
  return(res)
}

# filter PPIs by median & ppi_thresh
filter_ppi_by_median <- function(df_ppi = df_ppi, 
                                 sample_type = sample_type,  # vector
                                 cell_type_combs = cell_type_combs,
                                 ppi_thresh = .5){
  ppi <- rownames(df_ppi)
  
  st_lvl <- levels(as.factor(sample_type))
  
  if (length(st_lvl) != 2) {
    print("# sample type not equal to 2, please check parameter sample_type.")
    return(-1)
  } else {
    
    ppi_ctb <- lapply(cell_type_combs, function(ctb){
      dt <- df_ppi[, grep(ctb, colnames(df_ppi), fixed = T, value = F)]
      
      keep <- sapply(ppi, function(p){
        md <- sapply(c(1,2), function(i) median(t(dt[p, sample_type == st_lvl[i]])))
        return(md[1] > ppi_thresh | md[2] > ppi_thresh)
      })
      
      return(ppi[keep])
    })
    
    names(ppi_ctb) <- cell_type_combs
    return(ppi_ctb)
  }  # end of else
}

# filter PPIs by proportion of variance
filter_ppi_by_var_prop <- function(df_ppi = df_ppi, 
                                   cell_type_combs = cell_type_combs,
                                   var_prop_thresh = .5){
  ppi <- rownames(df_ppi)
  
  ppi_ctb <- lapply(cell_type_combs, function(ctb){
    dt <- df_ppi[, grep(ctb, colnames(df_ppi), fixed = T, value = F)]
    
    v <- sapply(ppi, function(p) var(t(dt[p,])))
    
    return(ppi[order(v, decreasing = T)][c(1:round(nrow(stat)*var_prop_thresh))])
  })
  
  names(ppi_ctb) <- cell_type_combs
  return(ppi_ctb)
}

# filter PPIs by variance threshold
filter_ppi_by_var_thresh <- function(df_ppi = df_ppi, 
                                     cell_type_combs = cell_type_combs,
                                     var_thresh = .5){
  ppi <- rownames(df_ppi)
  
  ppi_ctb <- lapply(cell_type_combs, function(ctb){
    dt <- df_ppi[, grep(ctb, colnames(df_ppi), fixed = T, value = F)]
    
    v <- sapply(ppi, function(p) var(t(dt[p,])))
    
    return(ppi[v >= var_thresh])
  })
  
  names(ppi_ctb) <- cell_type_combs
  return(ppi_ctb)
}

# merge ppi_filtered data
merge_ppi_filt <- function(ppi_filt1 = ppi_filt1,
                           ppi_filt2 = ppi_filt2,
                           cell_type_combs = cell_type_combs){
  p <- lapply(cell_type_combs, function(x) ppi_filt1[[x]][ppi_filt1[[x]] %in% ppi_filt2[[x]]])
  names(p) <- cell_type_combs
  return(p)
}

# ppi for each cell_type_combination
ppi_no_filter <- function(df_ppi = df_ppi, 
                          cell_type_combs = cell_type_combs){
  ppi <- rownames(df_ppi)
  
  ppi_ctb <- lapply(cell_type_combs, function(x) ppi)
  
  names(ppi_ctb) <- cell_type_combs
  return(ppi_ctb)
}

# permutation test for comparison of PPI scores between groups
perm_test_ppi <- function(df_ppi = df_ppi, 
                          sample_type = sample_type,  # vector
                          cell_type_combs = cell_type_combs,  # vector
                          n_resample = 1e4,
                          seed = 789){
  # verify # sample_type > 1
  if (length(levels(as.factor(sample_type))) == 1) {
    print("Only one sample_type found. Abort the procedure.")
    return(-1)
  }
  
  ppi <- rownames(df_ppi)
  
  #In the coin package, categorical variables and ordinal variables must be coded 
  # as factors and ordered factors, respectively
  sample_type <- as.factor(sample_type)
  
  set.seed(seed)
  
  # loop over cell_type_combs
  res_ctb <- lapply(cell_type_combs, function(ctb){
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi), fixed = T, value = F)]
    
    res_inter <- lapply(ppi, function(inter){
      test <- length(levels(as.factor(as.numeric(dt0[inter,])))) == 1
      
      if (test) {
        pval <- 1.0
      } else {
        fmla <- as.formula(paste0(inter, "~", "sample_type"))
        dt <- data.frame(t(dt0[inter,]), sample_type)
        res0 <- oneway_test(fmla, data=dt, distribution=approximate(nresample=n_resample))
        pval <- pvalue(res0)[[1]]
      }
      
      res = list(inter, pval)
      return(res)
    })
    
    return(res_inter)
  })
  
  names(res_ctb) <- cell_type_combs
  return(res_ctb)
}

# permutation test for comparison of PPI scores between groups [ppi filtered]
perm_test_filt_ppi <- function(df_ppi = df_ppi, 
                          sample_type = sample_type,  # vector
                          cell_type_combs = cell_type_combs,  # vector
                          ppi_filtered = ppi_filtered,  # list returned by filter_ppi_by_median
                          n_resample = 1e4,
                          seed = 789){
  # verify # sample_type > 1
  if (length(levels(as.factor(sample_type))) == 1) {
    print("Only one sample_type found. Abort the procedure.")
    return(-1)
  }
  
  #In the coin package, categorical variables and ordinal variables must be coded 
  # as factors and ordered factors, respectively
  sample_type <- as.factor(sample_type)
  
  set.seed(seed)
  
  # loop over cell_type_combs
  res_ctb <- lapply(cell_type_combs, function(ctb){
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi), fixed = T, value = F)]
    
    ppi <- ppi_filtered[[ctb]]
    
    res_inter <- lapply(ppi, function(p){
      test <- length(levels(as.factor(as.numeric(dt0[p,])))) == 1
      
      if (test) {
        pval <- 1.0
      } else {
        fmla <- as.formula(paste0(p, "~", "sample_type"))
        dt <- data.frame(t(dt0[p,]), sample_type)
        res0 <- oneway_test(fmla, data=dt, distribution=approximate(nresample=n_resample))
        pval <- pvalue(res0)[[1]]
      }
      
      return(list(p, pval))
    })
    
    return(res_inter)
  })
  
  names(res_ctb) <- cell_type_combs
  return(res_ctb)
}

# require multtest package
# adjust pvalue from perm_test_ppi or boot_t_ppi or perm_test_filt_ppi
adjust_pval <- function(res_test,  # result of test
                        proc = "BH", 
                        alpha = .05){
  ct_combs <- names(res_test)
  
  padj_ctb <- lapply(ct_combs, function(ctb){
    ppi <- sapply(c(1:length(res_test[[ctb]])), function(i) res_test[[ctb]][[i]][[1]])
    
    # vector of pvalue
    pval <- sapply(c(1:length(res_test[[ctb]])), function(i) res_test[[ctb]][[i]][[2]])
    
    # proc = c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", "TSBH")
    res_padj <- mt.rawp2adjp(pval, proc = proc, alpha = alpha, na.rm = T)
    
    # append 'contrast' column
    res_padj[[1]] <- as.data.frame(res_padj[[1]])
    res_padj[[1]] <- transform(res_padj[[1]],
                               contrast = ppi[res_padj[[2]]])
    
    return(res_padj)
  })
  
  names(padj_ctb) <- names(res_test)
  return(padj_ctb)
}

# bootstrap t-test for comparison of PPI scores between groups
boot_t_ppi <- function(df_ppi = df_ppi, 
                       sample_type = sample_type,  # vector
                       cell_type_combs = cell_type_combs,  # vector
                       alternative = "two.sided",  # c("two.sided", "less", "greater")
                       reps = 1e3,
                       seed = 789){
  # verify # sample_type > 1
  if (length(levels(as.factor(sample_type))) == 1) {
    print("Only one sample_type found. Abort the procedure.")
    return(-1)
  }
  
  ppi <- rownames(df_ppi)
  
  # sample_type
  st1 <- levels(as.factor(sample_type))[1]
  st2 <- levels(as.factor(sample_type))[2]
  
  set.seed(seed)
  
  # loop over cell_type_combs
  res_ctb <- lapply(cell_type_combs, function(ctb){
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi), fixed = T, value = F)]
    
    res_inter <- lapply(ppi, function(inter){
      test <- length(levels(as.factor(as.numeric(dt0[inter,])))) == 1
      
      if (test) {
        pval <- 1.0
      } else {
        dt <- dt0[inter,]
        x <- as.vector(t(dt[,sample_type == st1]))
        y <- as.vector(t(dt[,sample_type == st2]))
        pval <- boot.t.test(x, y, reps = reps, mu = 0, alternative = alternative)$p.value
      }
      
      return(list(inter, pval))
    })
    
    return(res_inter)
  })
  
  names(res_ctb) <- cell_type_combs
  return(res_ctb)
}

# bootstrap t-test for comparison of PPI scores between groups
boot_t_ppi_filt <- function(df_ppi = df_ppi, 
                       sample_type = sample_type,  # vector
                       cell_type_combs = cell_type_combs,  # vector
                       ppi_filtered = ppi_filtered,
                       alternative = "two.sided",  # c("two.sided", "less", "greater")
                       reps = 1e3,
                       seed = 789){
  # verify # sample_type > 1
  if (length(levels(as.factor(sample_type))) == 1) {
    print("Only one sample_type found. Abort the procedure.")
    return(-1)
  }
  
  # sample_type
  st1 <- levels(as.factor(sample_type))[1]
  st2 <- levels(as.factor(sample_type))[2]
  
  set.seed(seed)
  
  # loop over cell_type_combs
  res_ctb <- lapply(cell_type_combs, function(ctb){
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi), fixed = T, value = F)]
    
    ppi <- ppi_filtered[[ctb]]
    
    res_inter <- lapply(ppi, function(p){
      test <- length(levels(as.factor(as.numeric(dt0[p,])))) == 1
      
      if (test) {
        pval <- 1.0
      } else {
        dt <- dt0[p,]
        x <- as.vector(t(dt[,sample_type == st1]))
        y <- as.vector(t(dt[,sample_type == st2]))
        pval <- boot.t.test(x, y, reps = reps, mu = 0, alternative = alternative)$p.value
      }
      
      return(list(p, pval))
    })
    
    return(res_inter)
  })
  
  names(res_ctb) <- cell_type_combs
  return(res_ctb)
}
