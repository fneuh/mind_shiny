## helper script for web-app ##

feature_plot_fun <- function(seurat_obj, gene) {
  if(!gene %in% rownames(seurat_obj)) {
    return(NA)
  }
  FeaturePlot(seurat_obj, features = gene, reduction="umap")
}


feature_plot_by_stage_fun <- function(seurat_obj, gene) {
  if(!gene %in% rownames(seurat_obj)) {
    return(NA)
  }
  FeaturePlot(seurat_obj, features = gene, reduction="umap", split.by = "Stage", ncol = 3)
}


network_plot <- function(eRegulon_md_df, tf, mm10_tfs, only_tfs = TRUE) {
  ## subset md_df:
  eRegulon_sub <- eRegulon_md_df[eRegulon_md_df$TF == tf, ]
  
  if(only_tfs) {
    eRegulon_sub <- eRegulon_sub[eRegulon_sub$Gene %in% mm10_tfs, ]
  }
  
  ## create graph:
  edge_df <- eRegulon_sub[, c("TF", "Gene", "TF2G_regulation", "TF2G_importance_x_rho")]
  colnames(edge_df) <- c("to","from","TF2G_regulation", "TF2G_importance_x_rho")
  dg <- graph_from_data_frame(edge_df, directed = TRUE)
  
  ## simplify graph:
  dg <- igraph::simplify(dg, remove.multiple = T, remove.loops = T, 
                         edge.attr.comb = c(TF2G_importance_x_rho="mean", TF2G_regulation="mean"))
  
  ## graph visuals for plotting:
  E(dg)$width <- 1 + E(dg)$TF2G_importance_x_rho*0.3
  E(dg)$edge.arrow.size <- 1 + E(dg)$TF2G_importance_x_rho*0.3
  E(dg)$edge.arrow.width <- 1+ E(dg)$TF2G_importance_x_rho*0.5
  edge_color_vec <- c("1"= "blue","-1" = "red", "0" = "grey")
  E(dg)$color <- edge_color_vec[as.character(E(dg)$TF2G_regulation)]
  
  ## plot graph:
  plot(dg)
}