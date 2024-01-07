#!/usr/bin/env Rscript
# $title plot functions for intercellular PPI inference project
# $author giuseppe
# $date Apr 2020

library(ggplot2)
library(ComplexHeatmap)


# box & scatter plot
# plot PPI score with adjp of < .05
box_scatter_plt_ppi <- function(df_ppi = df_ppi,
                               sample_type = sample_type,
                               res_padj = res_padj,
                               cell_type_combs = cell_type_combs,
                               padj_thresh = .05,
                               ppi_score = ppi_score,
                               xlab = "Sample_type",
                               notch = F,
                               point.size = 2,
                               height = 10, width = 16,
                               theme_ggplot = theme_ggplot,
                               name_plt = name_plt,
                               directory = directory){
  for (ctb in cell_type_combs) {
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi))]
    
    filt <- res_padj[[ctb]][[1]][,2] < padj_thresh
    ppi <- res_padj[[ctb]][[1]][,3][filt]
    dt <- dt0[ppi,]
    
    # modify ctb for output
    ctb2 <- modify_ctb_name(ctb)
    dir_plt <- paste0(directory, ctb2, "/")
    dir.create(dir_plt)
    
    for (y in ppi) {
      dt_plt <- data.frame(t(dt[y,]), sample_type = sample_type)
      
      plot = ggplot(dt_plt, aes_string(x = "sample_type", y = y, color = "sample_type", fill = "sample_type")) +
        # argument "outlier.shape = NA" remove outlier in boxplot
        geom_boxplot(color = "black", alpha = .2, notch = notch, outlier.shape = NA) +
        geom_point(position = "jitter", size = point.size, alpha = .8) +
        labs(title=paste0(ppi_score, " of ", y), x=xlab, y=y) +
        theme_ggplot

      ggsave(filename = paste0(dir_plt,"Boxplot_",tolower(ppi_score),"_",tolower(ctb2),"_",y,"_",name_plt,".pdf"), 
             plot = plot, height = height, width = width, units = "cm")
    }  # end of for y
  }  # end of for ctb
}

# box & scatter plot 2
# ppi provided by ppi_filtered param
box_scatter_plt_ppi2 <- function(df_ppi = df_ppi,
                                 sample_type = sample_type,
                                 cell_type_combs = cell_type_combs,
                                 ppi_filtered = ppi_filtered,  # list returned by ppi_filt_logfc_padj func
                                 ppi_score = ppi_score,
                                 xlab = "Sample_type",
                                 notch = F,
                                 point.size = 2,
                                 height = 10, width = 16,
                                 theme_ggplot = theme_ggplot,
                                 name_plt = name_plt,
                                 directory = directory){
  for (ctb in cell_type_combs) {
    i <- match(ctb, cell_type_combs)
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi))]
    ppi <- ppi_filtered[[i]]
    dt <- dt0[ppi,]
    
    # modify ctb for output
    ctb2 <- modify_ctb_name(ctb)
    dir_plt <- paste0(directory, ctb2, "/")
    dir.create(dir_plt)
    
    for (y in ppi) {
      dt_plt <- data.frame(t(dt[y,]), sample_type = sample_type)
      
      plot = ggplot(dt_plt, aes_string(x = "sample_type", y = y, color = "sample_type", fill = "sample_type")) +
        # argument "outlier.shape = NA" remove outlier in boxplot
        geom_boxplot(color = "black", alpha = .2, notch = notch, outlier.shape = NA) +
        geom_point(position = "jitter", size = point.size, alpha = .8) +
        labs(title=paste0(ppi_score, " of ", y), x=xlab, y=y) +
        theme_ggplot
      
      ggsave(filename = paste0(dir_plt,"Boxplot_",tolower(ppi_score),"_",tolower(ctb2),"_",y,"_",name_plt,".pdf"), 
             plot = plot, height = height, width = width, units = "cm")
    }  # end of for y
  }  # end of for ctb
}


# scatter plot
# plot PPI score with adjp of < .05
scatter_plt_ppi <- function(df_ppi = df_ppi,
                               sample_type = sample_type,
                               res_padj = res_padj,
                               cell_type_combs = cell_type_combs,
                               padj_thresh = .05,
                               ppi_score = ppi_score,
                               xlab = "Sample_type",
                               point.size = 2,
                               height = 10, width = 16,
                               theme_ggplot = theme_ggplot,
                               name_plt = name_plt,
                               directory = directory){
  for (ctb in cell_type_combs) {
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi))]
    
    filt <- res_padj[[ctb]][[1]][,2] < padj_thresh
    ppi <- res_padj[[ctb]][[1]][,3][filt]
    dt <- dt0[ppi,]
    
    # modify ctb for output
    ctb2 <- modify_ctb_name(ctb)
    dir_plt <- paste0(directory, ctb2, "/")
    dir.create(dir_plt)
    
    for (y in ppi) {
      dt_plt <- data.frame(t(dt[y,]), sample_type = sample_type)
      
      plot = ggplot(dt_plt, aes_string(x = "sample_type", y = y, color = "sample_type", fill = "sample_type")) +
        geom_point(position = "jitter", size = point.size, alpha = .8) +
        labs(title=paste0(ppi_score, " of ", y), x=xlab, y=y) +
        theme_ggplot

      ggsave(filename = paste0(dir_plt,"Scatterplot_",tolower(ppi_score),"_",tolower(ctb2),"_",y,"_",name_plt,".pdf"), 
             plot = plot, height = height, width = width, units = "cm")
    }  # end of for y
  }  # end of for ctb
}

# violin plot
violin_plt_ppi <- function(df_ppi = df_ppi,
                           sample_type = sample_type,
                           res_padj = res_padj,
                           cell_type_combs = cell_type_combs,
                           padj_thresh = .05,
                           ppi_score = ppi_score,
                           xlab = "Sample_type",
                           point.size = 2,
                           height = 10, width = 16,
                           theme_ggplot = theme_ggplot,
                           name_plt = name_plt,
                           directory = directory){
  for (ctb in cell_type_combs) {
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi))]
    
    filt <- res_padj[[ctb]][[1]][,2] < padj_thresh
    ppi <- res_padj[[ctb]][[1]][,3][filt]
    dt <- dt0[ppi,]
    
    # modify ctb for output
    ctb2 <- modify_ctb_name(ctb)
    dir_plt <- paste0(directory, ctb2, "/")
    dir.create(dir_plt)
    
    for (y in ppi) {
      dt_plt <- data.frame(t(dt[y,]), sample_type = sample_type)
      
      plot = ggplot(dt_plt, aes_string(x = "sample_type", y = y, color = "sample_type", fill = "sample_type")) +
        geom_point(position = "jitter", size = point.size, alpha = .8) +
        geom_violin(color = "black", alpha = .2) +
        labs(title=paste0(ppi_score, " of ", y), x=xlab, y=y) +
        theme_ggplot
      
      ggsave(filename = paste0(dir_plt,"Violin_plot_",tolower(ppi_score),"_",tolower(ctb2),"_",y,"_",name_plt,".pdf"), 
             plot = plot, height = height, width = width, units = "cm")
    }  # end of for y
  }  # end of for ctb
}

# violin plot 2
# ppi provided by ppi_filtered param
violin_plt_ppi2 <- function(df_ppi = df_ppi,
                                 sample_type = sample_type,
                                 cell_type_combs = cell_type_combs,
                                 ppi_filtered = ppi_filtered,  # list returned by ppi_filt_logfc_padj func
                                 ppi_score = ppi_score,
                                 xlab = "Sample_type",
                                 notch = F,
                                 point.size = 2,
                                 height = 10, width = 16,
                                 theme_ggplot = theme_ggplot,
                                 name_plt = name_plt,
                                 directory = directory){
  for (ctb in cell_type_combs) {
    i <- match(ctb, cell_type_combs)
    dt0 <- df_ppi[, grep(ctb, colnames(df_ppi))]
    ppi <- ppi_filtered[[i]]
    dt <- dt0[ppi,]
    
    # modify ctb for output
    ctb2 <- modify_ctb_name(ctb)
    dir_plt <- paste0(directory, ctb2, "/")
    dir.create(dir_plt)
    
    for (y in ppi) {
      dt_plt <- data.frame(t(dt[y,]), sample_type = sample_type)
      
      plot = ggplot(dt_plt, aes_string(x = "sample_type", y = y, color = "sample_type", fill = "sample_type")) +
        geom_point(position = "jitter", size = point.size, alpha = .8) +
        geom_violin(color = "black", alpha = .2) +
        labs(title=paste0(ppi_score, " of ", y), x=xlab, y=y) +
        theme_ggplot
      
      ggsave(filename = paste0(dir_plt,"Violin_plot_",tolower(ppi_score),"_",tolower(ctb2),"_",y,"_",name_plt,".pdf"), 
             plot = plot, height = height, width = width, units = "cm")
    }  # end of for y
  }  # end of for ctb
}

# heatmap of PPI_score
# require ComplexHeatmap package
heatmap_ppi <- function(df_ppi = df_ppi,
                        sample_type = sample_type,
                        ref_sample_type = ref_sample_type,  # e.g. "control"
                        cell_type_combs = cell_type_combs,
                        adj_pvalue = adj_pvalue,  # list returned by adjust_pval
                        padj_thresh = .05,
                        log2fc_filtered = log2fc_filtered,  # list returned by log2fc_filter
                        name = name,  # name of legend
                        show_row_names = F,
                        show_column_names = F,
                        row_title = character(0),
                        row_title_side = "left",
                        column_title = character(0),
                        column_title_side = "top",
                        annotation_name_side = "left",
                        cluster_rows = T,
                        cluster_columns = T,
                        ppi_score = ppi_score,  # name of ppi_score
                        height = 6, 
                        width = 8,
                        pointsize = 14,  # for pdf function
                        name_plt = name_plt,
                        directory = directory){
  # convert label of sample_type
  sample_type <- ifelse(sample_type == ref_sample_type, "control", "tumour")
  
  for (ctb in cell_type_combs) {
    filt <- adj_pvalue[[ctb]]$adjp[,2] < padj_thresh
    ppi <- adj_pvalue[[ctb]]$adjp[,3][filt]
    ppi2 <- names(log2fc_filtered[[ctb]])
    ppi3 <- ppi[ppi %in% ppi2]
    
    print(paste0("cell_type_comb: ", ctb))
    print(paste0("# of PPI: ", length(ppi3)))
    
    mat <- as.matrix(df_ppi[ppi3, grep(ctb, colnames(df_ppi), fixed = T)])
    
    column_ha = HeatmapAnnotation(sample_type = sample_type,
                                  col = list(sample_type = c("control" = "purple", "tumour" = "turquoise")),
                                  annotation_legend_param = list(sample_type = list(nrow = 1)),
                                  annotation_name_side = annotation_name_side)
    
    ht <- Heatmap(mat, name = name, top_annotation = column_ha,
            show_row_names = show_row_names, show_column_names = show_column_names,
            row_title = row_title, column_title = column_title,
            cluster_rows = cluster_rows, cluster_columns = cluster_columns,
            heatmap_legend_param = list(direction = "horizontal"))
    
    plot <- draw(ht, merge_legend = T, heatmap_legend_side = "top", annotation_legend_side = "top")
    
    ctb2 <- modify_ctb_name(ctb)
    pdf(file = paste0(directory, "Heatmap_",tolower(ppi_score),"_",tolower(ctb2),"_",name_plt,".pdf"),
        width = width, height = height, pointsize = pointsize)
    print(plot)
    dev.off()
  }
}

# modify ctb for output
modify_ctb_name <- function(ctb){
  # convert "." to "_"
  ctb2 <- gsub("..", "_", ctb, fixed = T)
  ctb2 <- gsub("__", "_", ctb2, fixed = T)
  # remove leading & tailing "_"
  ctb2 <- gsub("^_", "", ctb2, fixed = F)  # regular expression
  ctb2 <- gsub("_$", "", ctb2, fixed = F)
  return(ctb2)
}
