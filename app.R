## web application using shiny ##

## load data and libraries ##
library(shiny)
library(bslib)
#library(periscope)

#library(Seurat)
library(igraph)
library(ggplot2)
library(pals)
library(data.table)
#library(Matrix)

## this needs to go to server.R or something don't really know actually..:
#library(Cairo)
#options(shiny.usecairo=T)

## read helper functions:
source("helper.R")


## read pre-processed tables:

## EI:
EI_umap_embedding <- fread("data/EI_merged_umap2_df.tsv", sep = "\t")
EI_feature_df <- fread("data/EI_merged_log_count_mtx.csv")


EI_umap_embedding <- as.data.frame(EI_umap_embedding)
EI_umap_embedding$cluster <- factor(EI_umap_embedding$cluster, levels = c(
  "Gas1_Ldha","Hes1_Fabp7","Fabp7_Mt3","Hist1h1b_Top2a","Ccnd2_Nudt4","Nkx2-1_Lhx8","Npy_Nxph1","Sst_Maf","Nr2f2_Nr2f1","Isl1_Zfp503","Foxp1_Gucy1a3","Ebf1_Foxp1",
  "Neurog2_Rrm2","Neurog2_Eomes","Neurod2_Neurod6","Neurod6_Mef2c"
))
EI_feature_df <- as.data.frame(EI_feature_df)


## INH:
INH_umap_embedding <- fread("data/inhibitory_datasets_umap2_df.tsv", sep = "\t")
INH_feature_df <- fread("data/inhibitory_datasets_log_count_mtx.csv")

INH_umap_embedding <- as.data.frame(INH_umap_embedding)
INH_umap_embedding$cluster <- factor(INH_umap_embedding$cluster, levels = c(
  "Fabp7","Top2a","Fabp7_Ccnd2","Ube2c","Nkx2_1","Abracl","Npy","Maf_Sst","Snhg11_Lhx8","Snhg11","Tcf4_Nr2f2",
  "Tshz1","Six3_Gucy1a3","Gucy1a3","Ebf1_Isl1"
))
INH_feature_df <- as.data.frame(INH_feature_df)

## combine umap and feature dfs for easier handling:
umap_embedding_list <- list(
  "INH" = INH_umap_embedding,
  "EI" = EI_umap_embedding
)
feature_list <- list(
  "INH" = INH_feature_df,
  "EI" = EI_feature_df
)

## GRN:
eRegulon_md_df <- read.table("/data/mayerlab/neuhaus/dorsal_ventral_comp/cfse_network/results_wArchRPeaks/eRegulon_metadata_filtered.tsv", sep = "\t", h=T)

mm10_tfs <- read.table("/datastore_share/Users/neuhaus/for_Ann/mm_mgi_tfs.txt"); mm10_tfs <- mm10_tfs$V1


## -----------------------------------------------------------------------------

## UI ##

ui <- page_navbar(
  title = "Mouse Inhibitory Neuron Development",
  bg = "#b2abd2",
  
  
  ## ---------------------------------------------------------------------------
  ## HOME ##
  nav_panel(
    title = "Home",
    ## HTML static content
    
    h2("Mouse Inhibitory Neuron Development"),
    
    layout_columns(
      card(
        p("This webserver accompanies this publication: "),
        a(href="https://www.biorxiv.org/content/10.1101/2024.03.18.585524v1", "Kotlyarenko & Bright et al. 2024"),
        h5("Abstract"),
        p("Inhibitory neurons of the telencephalon are generated from progenitors in the ganglionic eminences that mature and differentiate into specialized cell types. Here, we used single cell transcriptomics and single cell chromatin accessibility together with lineage tracing and birthdating techniques to investigate the influence of progenitor competence on the development of GABAergic precursors. We found that the timing of neurogenesis influences the maturation competence of progenitors to develop towards a fully functional state, but not their differentiation competence to evolve into transcriptomically diverse states. The underlying mechanism defining maturation competence was chromatin priming, orchestrated by the transcription factor Nfib in collaboration with regulators of inhibitory neuron development. Finally, transplantation experiments revealed an interplay between both intrinsic and extrinsic cues acting upon maturation competence. These findings identify a mechanism that coordinates inhibitory neuron development by changing its maturation to achieve maximum adaptability to their environment.")
      ),
      
      card(
        img(src = "ge_scheme.png")
      ),
      col_widths = c(7,3)
    ),
    
    p("We provide the following panels to explore our data:"),
    h5("Dataset UMAPs"),
    p("Explore structure of single-cell RNA-seq datasets. Visualize clusters, stages, experiments and studies. The data stems from this study (Kotylarenko et al. 2024) and two additional datasets: Development of inhibitory neurons from Bandler et al. 2022 and development of excitatory neurons in somato-sensory cortex from Di Bella et al. 2012."),
    h5("RNA expression"),
    p("Visualize expression of a gene of interest. UMAP plots can be split by cluster, stage, experiment and study."),
    h5("Tracks"),
    p("Link to UCSC session"),
    h5("Network"),
    p("Explore interactions between up to 3 TFs and their direct target genes.")
  ),
  
  
  ## ---------------------------------------------------------------------------
  ## DATA UMAP ##
  nav_panel(
    title = "Dataset UMAPs",
    
    layout_columns(
      card(
        helpText(
          "Create UMAP plot."
        ),
        selectInput(
          inputId = "umap_dataset",
          label = "Choose dataset:",
          choices = c("Inhibitory", "Inhibitory and Excitatory"),
          selected = "Inhibitory"
        ),
        selectInput(
          inputId = "umap_color_by",
          label = "Annotate cells by:",
          choices = c("stage", "cluster", "experiment", "study"),
          selected = "cluster"
        ),
        selectInput(
          inputId = "umap_split_by",
          label = "Split plots by:",
          choices = c("nothing, show combined", "stage", "cluster", "experiment","study"),
          selected = "nothing, show combined"
        ),
      submitButton("Update View", icon("refresh"))
      ),
      
      card(
        card_header("UMAP plot"),
        plotOutput("umap")
      ),
      
      col_widths = c(3,7)
    )
  ),
  
  ## ---------------------------------------------------------------------------
  ## RNA FEATURE PLOT ##
  nav_panel(
    title = "RNA Expression",
    
    layout_columns(
      card(
        helpText(
          "Create FeaturePlot for gene of interest."
        ),
        selectInput(
          inputId = "feature_dataset",
          label = "Choose dataset:",
          choices = c("Inhibitory", "Inhibitory and Excitatory"),
          selected = "Inhibitory"
        ),
        textInput(
          inputId = "gene",
          label = "Select your gene of interest:",
          value = "Nfib"
        ),
        selectInput(
          inputId = "feature_split_by",
          label = "Split plots by:",
          choices = c("nothing, show combined", "stage", "cluster", "experiment","study"),
          selected = "nothing, show combined"
        ),
        submitButton("Update View", icon("refresh"))
      ),
    
      card(
        card_header("Scaled Expression for Feature"),
        plotOutput("feature")
      ),
    
      col_widths = c(3,7)
    ),
  ),
  
  
  ## ---------------------------------------------------------------------------
  ## Tracks ##
  nav_panel(
    "Tracks", 
    p("First page content.")
  ),
  
  ## Network ##
  nav_panel(
    "Network",
    
    layout_columns(
      
      card(
        helpText(
          "Create subnetwork for TFs of interest."
        ),
        textInput(
          inputId = "tf1",
          label = "Type your TF of interest: ",
          value = "Nfib"
        ),
        textInput(
          inputId = "tf2",
          label = "Type another TF of interest: ",
          value = "not chosen"
        ),
        textInput(
          inputId = "tf3",
          label = "Type a third TF of interest: ",
          value = "not chosen"
        ),
        selectInput(
          inputId = "only_TF",
          label = "Plot only TFs?",
          choices = c("Yes", "No"),
          selected = "Yes"
        ),
        submitButton("Update View", icon("refresh"))
      ),
      
      card(
        plotOutput("network")
      ),
      
      col_widths = c(3,7)
    )
  ),
  
  nav_spacer()
)


## -----------------------------------------------------------------------------
## server logic ##

server <- function(input, output) {
  ## 1: UMAP
  output$umap <- renderPlot({
    umap_plot_ggplot(df_list = umap_embedding_list, dataset = input$umap_dataset, col_attr = input$umap_color_by, split_attr = input$umap_split_by, point_size = 0.2)
  })
  
  ## 2: FEATURE 
  output$feature <- renderPlot({
    feature_plot_ggplot(df_list = umap_embedding_list, mtx_list = feature_list, dataset = input$feature_dataset, gene_name = input$gene, split_attr = input$feature_split_by, point_size = 0.2) 
  })

  ## 3: TRACKS
  
  ## todo
  
  ## 4: NETWORK
  output$network <- renderPlot({
    only_tfs <- switch(input$only_TF,
                       "Yes" = TRUE,
                       "No" = FALSE)
    network_plot(
      eRegulon_md_df = eRegulon_md_df,
      tf1 = input$tf1,
      tf2 = input$tf2,
      tf3 = input$tf3,
      mm10_tfs = mm10_tfs,
      only_tfs = only_tfs
    )
  })
  # output$network_plot_download <- downloadHandler(
  #   filename = function() { "network_plot.pdf" },
  #   content = function(file) {
  #     pdf(file, paper = "default")
  #     plot({
  #       only_tfs <- switch(input$only_TF,
  #                          "Yes" = TRUE,
  #                          "No" = FALSE)
  #       network_plot(
  #         eRegulon_md_df = eRegulon_md_df,
  #         tf1 = input$tf1,
  #         tf2 = input$tf2,
  #         tf3 = input$tf3,
  #         mm10_tfs = mm10_tfs,
  #         only_tfs = only_tfs
  #       )
  #     })
  #     dev.off()
  #   }
  # )
}


## run app ##

shinyApp(ui, server)