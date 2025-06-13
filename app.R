## web application using shiny ##

## load data and libraries ##
library(shiny)
library(bslib)

library(igraph)
library(ggplot2)
library(pals)
library(data.table)
library(rhdf5)

## read helper functions:
source("helper.R")

## read pre-processed tables:
## EI:
EI_umap_embedding <- fread("data/EI_merged_umap2_df.tsv", sep = "\t")
#EI_feature_df <- fread("data/EI_merged_log_count_sub_mtx.csv")
EI_feature_file <- H5Fopen("data/EI_merged_datasets.h5")

EI_umap_embedding <- as.data.frame(EI_umap_embedding)
EI_umap_embedding$cluster <- factor(EI_umap_embedding$cluster, levels = c(
  "Gas1_Ldha","Hes1_Fabp7","Fabp7_Mt3","Hist1h1b_Top2a","Ccnd2_Nudt4","Nkx2-1_Lhx8","Npy_Nxph1","Sst_Maf","Nr2f2_Nr2f1","Isl1_Zfp503","Foxp1_Gucy1a3","Ebf1_Foxp1",
  "Neurog2_Rrm2","Neurog2_Eomes","Neurod2_Neurod6","Neurod6_Mef2c"
))
EI_umap_embedding$class <- factor(EI_umap_embedding$class, levels = "Mitotic", "Inhibitory Neuron Precursor", "Excitatory Neuron Precursor")
EI_umap_embedding$study <- factor(EI_umap_embedding$study, levels = c("Bright et al. 2025","Bandler et al. 2022","Di Bella et al. 2021"))

#EI_feature_df <- as.data.frame(EI_feature_df)


## INH EMBRYONIC:
INH_umap_embedding <- fread("data/inhibitory_datasets_umap2_df.tsv", sep = "\t")
#INH_feature_df <- fread("data/inhibitory_datasets_log_count_sub_mtx.csv")
INH_feature_file <- H5Fopen("data/inhibitory_datasets.h5")

INH_umap_embedding <- as.data.frame(INH_umap_embedding)
INH_umap_embedding$cluster <- factor(INH_umap_embedding$cluster, levels = c(
  "Fabp7","Top2a","Fabp7_Ccnd2","Ube2c","Nkx2_1","Abracl","Npy","Maf_Sst","Snhg11_Lhx8","Snhg11","Tcf4_Nr2f2",
  "Tshz1","Six3_Gucy1a3","Gucy1a3","Ebf1_Isl1"
))
INH_umap_embedding$experiment <- factor(INH_umap_embedding$experiment, levels = c("WT","CFSE","LINEAGE"))
INH_umap_embedding$stage <- factor(INH_umap_embedding$stage, levels = c("E12","E14","E16","P0"))

#INH_feature_df <- as.data.frame(INH_feature_df)

## INH POSTNATAL:
INH_PN_umap_embedding <- fread("data/STICR_umap_df.tsv", sep = "\t")
#INH_PN_feature_df <- fread("data/STICR_log_count_sub_mtx.csv")
INH_PN_feature_file <- H5Fopen("data/STICR_datasets.h5")

INH_PN_umap_embedding <- as.data.frame(INH_PN_umap_embedding)
INH_PN_umap_embedding$stage <- factor(INH_PN_umap_embedding$stage, levels = c("P5","P7","P8","P9","P11","P13/P14","P15"))
INH_PN_umap_embedding$experiment <- factor(INH_PN_umap_embedding$experiment, levels = c(
  "STICR E10.5 - P8","STICR E10.5 - P9","STICR E10.5 - P11",
  "STICR E12.5 - P7","STICR E12.5 - P8","STICR E12.5 - P9","STICR E12.5 - P15",
  "STICR E13.5 - P5","STICR E13.5 - P7","STICR E13.5 - P13/P14",
  "STICR E14.5 - P7"
))
#INH_PN_feature_df <- as.data.frame(INH_PN_feature_df)


## combine umap and feature dfs for easier handling:
umap_embedding_list <- list(
  "INH_EMB" = INH_umap_embedding,
  "EI_EMB" = EI_umap_embedding,
  "INH_PN" = INH_PN_umap_embedding
)
feature_file_list <- list(
  "INH_EMB" = INH_feature_file,
  "EI_EMB" = EI_feature_file,
  "INH_PN" = INH_PN_feature_file
)

all_gene_names <- unique(c(
  as.vector(INH_feature_file$"gene_names/gene_name_mtx"),
  as.vector(EI_feature_file$"gene_names/gene_name_mtx"),
  as.vector(INH_PN_feature_file$"gene_names/gene_name_mtx")
))
all_gene_names <- all_gene_names[order(all_gene_names)]

## GRN:
eRegulon_md_df <- read.table("data/eRegulon_metadata_filtered.tsv", sep = "\t", h=T)
mm10_tfs <- read.table("data/mm_mgi_tfs.txt"); mm10_tfs <- mm10_tfs$V1


## -----------------------------------------------------------------------------

## UI ##

ui <- page_navbar(
  title = "Mouse Inhibitory Neuron Development",
  bg = "#b2abd2",
  tags$head(tags$link(rel = "icon", type = "image/png", href = "CM_2025Logo.png")),
  
  ## ---------------------------------------------------------------------------
  ## HOME ##
  nav_panel(
    title = "Home",
    ## HTML static content
    
    h2("Mouse Inhibitory Neuron Development"),
    
    layout_column_wrap(
      card(
        p("This webserver accompanies this publication: "),
        a(href="https://www.biorxiv.org/content/10.1101/2024.03.18.585524v2", "Bright, Kotlyarenko & Neuhaus et al. 2025"),
        h4("About this resource"),
        p("This browser enables interactive exploration of published single-cell RNA sequencing datasets covering mouse forebrain development, with a focus on inhibitory neuron populations."),
        p("It includes:"),
        
        tags$ul(tags$li("Embryonic datasets from Bright, Kotlyarenko & Neuhaus et al. (2025), Di Bella & Habibi et al. (2021) and Bandler, Vitali & Delgado et al. (2021)"), tags$li("Neonatal dataset from Bandler, Vitali & Delgado et al. (2021)")),
        
        p("Together, these datasets span key stages of inhibitory neuron specification, migration, and maturation, across both dorsal and ventral forebrain structures. The browser provides access to gene expression profiles and gene regulatory network predictions."),
        
      ),
      card(
        img(src = "CM_FN_logo.png"),
      )
    ),
    
    #p("This webserver accompanies this publication: "),
    #a(href="https://www.biorxiv.org/content/10.1101/2024.03.18.585524v2", "Bright, Kotylarenko & Neuhaus et al. 2025"),
    #h4("About this resource"),
    #p("This browser enables interactive exploration of published single-cell RNA sequencing datasets covering mouse forebrain development, with a focus on inhibitory neuron populations."),
    #p("It includes:"),
    #tags$ul(tags$li("Embryonic datasets from Bright, Kotlyarenko & Neuhaus et al. (2025), Di Bella & Habibi et al. (2021) and Bandler, Vitali & Delgado et al. (2021)"), tags$li("Neonatal dataset from Bandler, Vitali & Delgado et al. (2021)")),
    #p("Together, these datasets span key stages of inhibitory neuron specification, migration, and maturation, across both dorsal and ventral forebrain structures. The browser provides access to gene expression profiles and gene regulatory network predictions."),
    
    
    # layout_columns(
    #   card(
    #     p("This webserver accompanies this publication: "),
    #     a(href="https://www.biorxiv.org/content/10.1101/2024.03.18.585524v2", "Bright, Kotylarenko & Neuhaus et al. 2025"),
    #     h5("Abstract"),
    #     p("Diverse types of GABAergic projection neurons and interneurons of the telencephalon derive from progenitors in a ventral germinal zone, called the ganglionic eminence. Using single-cell transcriptomics, chromatin accessibility profiling, lineage tracing, birthdating, heterochronic transplantation, and perturbation sequencing in mouse embryos, we investigated how progenitor competence influences the maturation and differentiation of these neurons. We found that the progression of neurogenesis over developmental time shapes maturation competence in ganglionic eminence progenitors, influencing how they progress into mature states. In contrast, differentiation competence, which defines the ability to produce diverse transcriptomic identities, remains largely unaffected by the stages of neurogenesis. Chromatin remodeling alongside a NFIB-driven regulatory gene module influences maturation competence in late-born neurons. These findings provide key insights into how transcriptional programs and chromatin accessibility govern neuronal maturation and the diversification of GABAergic neuron subtypes during neurodevelopment.")
    #   ),
    # 
    #   card(
    #     img(src = "ge_scheme.png")
    #   ),
    #   col_widths = breakpoints(
    #     sm = c(4,1),
    #     md = c(6,3),
    #     lg = c(8,4)
    #   ),
    #   row_heights = breakpoints(
    #     sm = c(20),
    #     md = c(16),
    #     lg = c(12)
    #   )
    # ),
    
    p("We provide the following panels to explore our data:"),
    
    tags$ul(
      tags$li(tags$b("Dataset UMAPs:"), " Explore structure of single-cell RNA-seq datasets. Visualize cell clusters, cell classes, developmental stages, experiments and studies. The data stems from this study and two additional datasets: Development of inhibitory neurons from Bandler, Vitali & Delgado et al. (2021) and development of excitatory neurons in somatosensory cortex from Di Bella & Habibi et al. (2021)."),
      tags$li(tags$b("RNA expression:"), " Visualize expression of a gene of interest. UMAP plots can be split by cell cluster, cell class, developmental stage, experiment and study. Genes were filtered for genes that are either highly-variable (n=5000) or having log-normalized expression greater than 0.5."),
      tags$li(tags$b("Tracks:"), " Data from scATAC-seq and NFIB CUT&RUN experiments can be explored in UCSC genome browser"),
      tags$li(tags$b("Network:"), " Explore interactions between up to 3 TFs and their direct target genes. Gene-regulatory networks were precdited using SCENIC+: González-Blas & De Winter et al. (2023).")
    ),
    
    #p("Explore structure of single-cell RNA-seq datasets. Visualize cell clusters, cell classes, developmental stages, experiments and studies. The data stems from this study and two additional datasets: Development of inhibitory neurons from Bandler, Vitali & Delgado et al. (2021) and development of excitatory neurons in somatosensory cortex from Di Bella & Habibi et al. (2021)."),
    #p("RNA expression"),
    #p("Visualize expression of a gene of interest. UMAP plots can be split by cell cluster, cell class, developmental stage, experiment and study. Genes were filtered for genes that are either highly-variable (n=5000) or having log-normalized expression greater than 0.5."),
    #p("Tracks"),
    #p("Data from scATAC-seq and NFIB CUT&RUN experiments can be explored in UCSC genome browser"),
    #p("Network"),
    #p("Explore interactions between up to 3 TFs and their direct target genes. Gene-regulatory networks were precdited using SCENIC+: González-Blas & De Winter et al. (2023).")
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
          choices = c("Embryonic Inhibitory", "Embryonic Inhibitory and Excitatory","Neonatal Inhibitory"),
          selected = "Embryonic Inhibitory"
        ),
        selectInput(
          inputId = "umap_color_by",
          label = "Annotate cells by:",
          choices = c("stage", "cluster", "experiment", "study","class"),
          selected = "cluster"
        ),
        selectInput(
          inputId = "umap_split_by",
          label = "Split plots by:",
          choices = c("nothing, show combined", "stage", "cluster", "experiment","study","class"),
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
          choices = c("Embryonic Inhibitory", "Embryonic Inhibitory and Excitatory","Neonatal Inhibitory"),
          selected = "Embryonic Inhibitory"
        ),
        selectInput(
          inputId = "gene",
          label = "Select your gene of interest:",
          choices = all_gene_names,
          selected = "Nfib"
        ),
        selectInput(
          inputId = "feature_split_by",
          label = "Split plots by:",
          choices = c("nothing, show combined", "stage", "cluster", "experiment","study","class"),
          selected = "nothing, show combined"
        ),
        submitButton("Update View", icon("refresh"))
      ),
    
      card(
        card_header("Normalized Expression for Feature"),
        plotOutput("feature")
      ),
    
      col_widths = c(3,7)
    ),
  ),
  
  
  ## ---------------------------------------------------------------------------
  ## Tracks ##
  nav_panel(
    "Tracks", 
    p("Data from scATAC-seq experiments (e12.5 and e16.5) together with data from NFIB CUT&RUN experiments can be explored in UCSC genome browser session:"),
    a(href="https://genome.ucsc.edu/s/annrosebright/MIND%20Trackhub", "LINK to genome browser"),
    p("Available tracks include chromatin accessibility at e12.5 and e16.5 for all cells and split by broad cell states (AP, BP, precursor). For CUT&RUN experiments we show 2 tracks: NFIB and H3K4me3")
  ),
  
  ## Network ##
  nav_panel(
    "Network",
    
    layout_columns(
      
      card(
        helpText(
          "Create subnetwork for TFs of interest."
        ),
        selectInput(
          inputId = "tf1",
          label = "Type your TF of interest: ",
          selected = "Nfib",
          choices = unique(eRegulon_md_df$TF)
        ),
        selectInput(
          inputId = "tf2",
          label = "Type another TF of interest: ",
          selected = "not chosen",
          choices = c("not chosen", unique(eRegulon_md_df$TF))
        ),
        selectInput(
          inputId = "tf3",
          label = "Type a third TF of interest: ",
          selected = "not chosen",
          choices = c("not chosen", unique(eRegulon_md_df$TF))
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
    feature_plot_ggplot(df_list = umap_embedding_list, mtx_list = feature_file_list, dataset = input$feature_dataset, gene_name = input$gene, split_attr = input$feature_split_by, point_size = 0.2) 
  })

  ## 3: TRACKS
  
  ## static content
  
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