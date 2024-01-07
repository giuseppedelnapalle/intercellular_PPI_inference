#!/usr/bin/env Rscript
# $title handy functions for intercellular PPI inference projects
# $author giuseppe
# $date Apr 2020


# 1 string processing -----------------------------------------------------


# remove leading n fields in a string
remove_n_fileds <- function(x, 
                            sep = "_", 
                            fixed = T){
  str_rm <- sapply(strsplit(x, split = sep, fixed = fixed), "[", 1)
  y = gsub(paste0(str_rm, sep), "", x, fixed = fixed)
  return(y)
}

# convert string containing cell_type_comb to vector
ct_comb_vector <- function(x, 
                           sep, 
                           m, n, 
                           fixed = T){
  elem1 <- sapply(strsplit(x, sep, fixed = fixed), "[", m)
  elem2 <- sapply(strsplit(x, sep, fixed = fixed), "[", n)
  y <- c(elem1, elem2)
  return(y)
}

# convert "." to "_", strip leading & tailing "." in name
modify_name <- function(name_vector = name_vector){
  nm2 <- sapply(name_vector, function(x){
    # convert "." to "_"
    y <- gsub("..", "_", x, fixed = T)
    y <- gsub("__", "_", y, fixed = T)
    # remove leading & tailing "_"
    y <- gsub("^_", "", y, fixed = F)  # regular expression
    y <- gsub("_$", "", y, fixed = F)
    return(y)
  })
}


# 2 filtering -------------------------------------------------------------

# convert values of < thresh to set value
filter_ppi <- function(mat, thresh = .5, min_val = .01){
  mat_filt <- mat
  mat_filt[mat < thresh] <- min_val
  return(mat_filt)
}

# log2FoldChange of PPI score
log2fc <- function(df_ppi = df_ppi,  # PPI as rows, sample_cell_type_pair as columns
                   sample_type = sample_type,  # vector of sample_type
                   ref_sample_type = ref_sample_type,  # str
                   cell_type_combs = cell_type_combs){
  ppi <- rownames(df_ppi)
  
  res_ctb <- lapply(cell_type_combs, function(ctb){
    dt0 <- df_ppi[,grep(ctb, colnames(df_ppi), fixed = T)]
    
    if (ncol(dt0) == length(sample_type)) {
      logfc_ppi <- sapply(ppi, function(p){
        dt <- dt0[p,]
        fc <- mean(t(dt[,sample_type != ref_sample_type])) / 
          mean(t(dt[,sample_type == ref_sample_type]))
        logfc <- log2(fc)
        return(logfc)
      })
      return(logfc_ppi)
      
    } else {
      print(paste0("Wrong number of samples of cell_type_combs ", ctb))
      return(NULL)
    }
  })
  
  names(res_ctb) <- cell_type_combs
  return(res_ctb)
}

# filter log2fc by threshold
log2fc_filter <- function(log2fc_list = log2fc_list,
                          thresh_logfc = .8){
  log_filt <- lapply(names(log2fc_lst), function(nm){
    filt <- abs(log2fc_lst[[nm]]) >= thresh_logfc
    return(log2fc_lst[[nm]][filt])
  })
  names(log_filt) <- names(log2fc_list)
  return(log_filt)
}

# filter PPI by log2fc & padj
ppi_filt_logfc_padj <- function(adj_p = adj_p,
                                log2fc_list = log2fc_list,
                                padj_thresh = .05,
                                logfc_thresh = .8){
  log2fc_filt <- log2fc_filter(log2fc_list, logfc_thresh)
  
  ppi_filt <- lapply(names(adj_p), function(nm){
    x <- adj_p[[nm]][[1]][,3][adj_p[[nm]][[1]][,2] < padj_thresh]
    y <- names(log2fc_filt[[nm]])
    z <- x[x %in% y]
    print(paste0(nm, ": ", length(z)))
    return(z)
  })
  
  names(ppi_filt) <- modify_name(names(adj_p))
  return(ppi_filt)
}


