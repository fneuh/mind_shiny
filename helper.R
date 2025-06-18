## helper script for web-app ##

umap_plot_ggplot <- function(df_list, dataset, col_attr, split_attr = NA, point_size = 1) {
  if(split_attr == "nothing, show combined") {split_attr <- NA}

  ## select dataset:
  if(dataset == "Embryonic Inhibitory") {
    df <- df_list$INH_EMB
  } else if(dataset == "Embryonic Inhibitory and Excitatory") {
    df <- df_list$EI_EMB
  } else if (dataset == "Neonatal Inhibitory") {
    df <- df_list$INH_PN
  } else {
    print("No valid dataset specified")
    df <- NULL
  }

  ## plot aesthetics:
  if(length(unique(df[, col_attr])) <= 26) {
    col_vec <- alphabet2(n = length(unique(df[, col_attr])))
  } else {
    col_vec <- c(alphabet2(n = 26), alphabet(n = length(unique(df[, col_attr])) - 26))
  }
  
  names(col_vec) <- unique(df[, col_attr])
  
  g <- ggplot(df, aes(x = UMAP2_1, y = UMAP2_2, color = !! sym(col_attr))) +
    geom_point(size = point_size) +
    scale_color_manual(values = col_vec) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.text = element_text(size=15), legend.title=element_text(size=15), plot.title = element_text(size = 15)) +
    ggtitle(col_attr) +
    xlab("UMAP1") + ylab("UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 2)))
  if(!is.na(split_attr)) {
    g <- g +
      facet_wrap(sym(split_attr))
  }
  return(g)
}

feature_plot_ggplot <- function(df_list, mtx_list, dataset, gene_name, split_attr = NA, point_size = 1) {
  if(split_attr == "nothing, show combined") {split_attr <- NA}

  ## choose datasets
  if(dataset == "Embryonic Inhibitory") {
    df <- df_list$INH_EMB
    mtx_file <- mtx_list$INH_EMB
  } else if(dataset == "Embryonic Inhibitory and Excitatory") {
    df <- df_list$EI_EMB
    mtx_file <- mtx_list$EI_EMB
  } else if (dataset == "Neonatal Inhibitory") {
    df <- df_list$INH_PN
    mtx_file <- mtx_list$INH_PN
  } else {
    print("No valid dataset specified")
    df <- NULL
    mtx_file <- NULL
  }
  
  ## add expression to df:
  gene_name_vec <- as.vector(mtx_file$"gene_names/gene_name_mtx")
  if(gene_name %in% gene_name_vec) {
    gene_idx <- which(gene_name_vec == gene_name)
  } else {
    gene_idx <- NULL
    return("Gene not found in this dataset.")
  }
  
  gene_expr_vec <- h5read(file = mtx_file, name = "log_expression/log_expr_mtx", index = list(NULL, gene_idx))
  df[, gene_name] <- gene_expr_vec
  
  g <- ggplot(df, aes(x = UMAP2_1, y = UMAP2_2, color = !! sym(gene_name))) +
    geom_point(size = point_size) +
    scale_color_gradient(low = "#ffffbf", high = "#a50026") +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.title=element_text(size=15), plot.title = element_text(size = 15)) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggtitle(gene_name)
  if(!is.na(split_attr)) {
    g <- g +
      facet_wrap(sym(split_attr))
  }
  return(g)
}



network_plot <- function(eRegulon_md_df, tf1, tf2 = NA, tf3 = NA, mm10_tfs, only_tfs = TRUE) {
  tf_vec <- c(tf1, tf2, tf3)
  
  ## subset md_df:
  eRegulon_sub <- eRegulon_md_df[eRegulon_md_df$TF %in% tf_vec, ]
  
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
  E(dg)$width <- 1 + E(dg)$TF2G_importance_x_rho*0.2
  E(dg)$edge.arrow.size <- 1 + E(dg)$TF2G_importance_x_rho*0.2
  E(dg)$edge.arrow.width <- 1+ E(dg)$TF2G_importance_x_rho*0.4
  edge_color_vec <- c("1"= "blue","-1" = "red", "0" = "grey")
  E(dg)$color <- edge_color_vec[as.character(E(dg)$TF2G_regulation)]
  
  ## plot graph:
  plot(dg)
  legend("topright", legend = c("activating","repressing"), lty = 1, col = c("blue","red"), lwd = 4)
}
