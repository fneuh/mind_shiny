## web application using shiny ##

## load data and libraries ##
library(shiny)
library(bslib)

library(Seurat)
library(igraph)

source("helper.R")

#load("/datastore_share/Users/ykotlyarenko/Inhibitory_cells_development/6_All_inhibitory_datasets/Datasets/Inhibitory_datasets.Rdata")

tcf4_ko_interneuron_sub <- readRDS("/datastore_share/Users/neuhaus/Tcf4_KO_cortex_scRNA/results/tcf4_ko_interneuron_sub_seurat.rds")
tcf4_ko_interneuron_sub$Stage <- sample(x = c("e12","e14","e16"), size = ncol(tcf4_ko_interneuron_sub), replace = T)

eRegulon_md_df <- read.table("/data/mayerlab/neuhaus/dorsal_ventral_comp/cfse_network/results_wArchRPeaks/eRegulon_metadata_filtered.tsv", sep = "\t", h=T)

mm10_tfs <- read.table("/datastore_share/Users/neuhaus/for_Ann/mm_mgi_tfs.txt"); mm10_tfs <- mm10_tfs$V1

## UI ##

ui <- page_navbar(
  title = "Mouse Inhibitory Neuron Development",
  bg = "#b2abd2",
  
  ## HOME ##
  nav_panel(
    title = "Home",
    ## HTML static content
    h2("Mouse Inhibitory Neuron Development"),
    h5("Introduction"),
    p("This webserver contains data from: "),
    a(href="https://www.biorxiv.org/content/10.1101/2024.03.18.585524v1", "Kotlyarenko & Bright et al. 2024"),
    p("In order to explore different data modalities we separated this webserver into 3 parts."),
    h5("RNA expression"),
    p("Explore RNA transcription levels of your genes of interest in a merged dataset spanning the whole of inhibitory neuron development. Data shown here contains WT data from e12.5, e14.5 and e16.5; CFSE-data from e12.5+6h, e12.5+96h and e16.5+6h. Additionally we also included WT data from Bandler et al."),
    h5("Tracks"),
    p("Link to UCSC session"),
    h5("Network"),
    p("Create Sub-graphs for direct targets of TF of choice.")
  ),
  
  ## RNA ##
  nav_panel(
    title = "RNA expression",
    sidebarPanel(
      helpText(
        "Create FeaturePlot for gene of interest."
      ),
      textInput(
        inputId = "gene",
        label = "Type your gene of interest:"
      ),
      #selectInput(
      #  inputId = "split_by_stage",
      #  label = "Split plots by stage?",
      #  choices = c("Yes", "No"),
      #  selected = "Yes"
      #),
      selectInput(
        inputId = "dataset",
        label = "Which dataset to use?",
        choices = c("Combined", "WT", "CFSE", "Bandler et al."),
        selected = "Combined"
      )
    ),
    submitButton("Update View", icon("refresh")),
    mainPanel(
      card(
        card_header("Expression FeaturePlot"),
        plotOutput("umap", width = 300, height = 300)
      ),
      card(
        card_header("Expression FeaturePlot by Stage"),
        plotOutput("umap_by_stage", width = 900, height = 300)
      )
    )
  ),
  
  ## Tracks ##
  nav_panel(
    "Tracks", 
    p("First page content.")
  ),
  
  ## Network ##
  nav_panel(
    "Network", 
    sidebarPanel(
      helpText(
        "Create subnetwork for TF of interest."
      ),
      textInput(
        inputId = "tf",
        label = "Type your TF of interest: "
      ),
      selectInput(
        inputId = "only_TF",
        label = "Plot only TFs?",
        choices = c("Yes", "No"),
        selected = "Yes"
      ),
      width = 3
    ),
    submitButton("Update View", icon("refresh")),
    card(plotOutput("network", width = 1500, height = 1500), height = "600px")
  ),
  
  nav_spacer()
)



## server logic ##

server <- function(input, output) {
  ## 1: RNA
  output$umap <- renderPlot({
    feature_plot_fun(seurat_obj = tcf4_ko_interneuron_sub, gene = input$gene) 
  })
  
  output$umap_by_stage <- renderPlot(({
    feature_plot_by_stage_fun(seurat_obj = tcf4_ko_interneuron_sub, gene = input$gene)
  }))
  ## 2: Tracks
  
  ## todo
  
  ## 3: Network
  output$network <- renderPlot({
    only_tfs <- switch(input$only_TF,
                       "Yes" = TRUE,
                       "No" = FALSE)
    network_plot(
      eRegulon_md_df = eRegulon_md_df,
      tf = input$tf,
      mm10_tfs = mm10_tfs,
      only_tfs = only_tfs
    )
  })
}


## run app ##

shinyApp(ui, server)