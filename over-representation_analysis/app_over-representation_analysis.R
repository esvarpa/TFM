library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(shinycssloaders)
library(DT)
library(shiny)
library(clusterProfiler)
library(DOSE)
library(GOplot)
library(enrichplot)
library(pathview)
library(plotly)
library(debrowser)
library(GO.db)
library(DESeq2)
library(seuratOnline)
library(readr)
library(Quandl)
library(wordcloud2)
library(ggupset)
library(ggnewscale)
library(GOSemSim)
library(BiocGenerics)
library(AnnotationDbi)
library(BiocParallel)
library(europepmc)
library(biomaRt)

# BiocInstaller::biocLite(c("org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Dm.eg.db","org.At.tair.db","org.Dr.eg.db","org.Bt.eg.db","org.Ce.eg.db","org.Gg.eg.db","org.Cf.eg.db","org.Ss.eg.db","org.Mmu.eg.db","org.EcK12.eg.db","org.Xl.eg.db","org.Pt.eg.db","org.Ag.eg.db","org.Pf.plasmo.db","org.EcSakai.eg.db"))

jsCode <- "shinyjs.pageCol = function(params){$('body').css('background', params);}"

ui <- tagList(
  dashboardPage(
    skin = "purple",
    dashboardHeader(title = "ClusterProfShiny (ORA)"),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("User Guide", tabName = "introTab", icon = icon("info-circle")),
        menuItem("Input Data", tabName = "datainput", icon = icon("upload")),
        menuItem("enrichGo", tabName = "enrichGoTab", icon = icon("th")),
        menuItem("enrichKegg", tabName = "enrichKeggTab", icon = icon("th")),
        menuItem("Go Plots", tabName = "goplotsTab", icon = icon("bar-chart")),
        menuItem("KEGG Plots", tabName = "keggPlotsTab", icon = icon("bar-chart")),
        menuItem("Pathview Plots", tabName = "pathviewTab", icon = icon("bar-chart")),
        menuItem("Word Clouds", tabName = "wordcloudTab", icon = icon("bar-chart"))
      )
    ),
    dashboardBody(
      shinyjs::useShinyjs(),
      extendShinyjs(text = jsCode, functions = c("pageCol")),
      tags$head(
        tags$style(HTML(
          " .shiny-output-error-validation {color: darkred; }"
        )),
        tags$style(
          type="text/css",
          "#pathview_plot img {max-width: 100%; width: 100%; height: auto}"),
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
      ),
      tabItems(
        source("ui-tab-intro.R", local = TRUE)$value,
        source("ui-tab-inputdata.R", local = TRUE)$value,
        source("ui-tab-enrichGo.R", local = TRUE)$value,
        source("ui-tab-enrichKegg.R", local = TRUE)$value,
        source("ui-tab-goplots.R", local = TRUE)$value,
        source("ui-tab-keggplots.R", local = TRUE)$value,
        source("ui-tab-pathview.R", local = TRUE)$value,
        source("ui-tab-wordcloud.R", local = TRUE)$value
      )
      
    )
    
  ),
  tags$footer(
    wellPanel(
      HTML(
        '
        <p align="center" width="4">Vall dHebron Research Institur</p>
        <p align="center" width="4">Statistics and Bioinformatics Research Group</a> </p>
        <p align="center" width="4">Copyright (C) 2019, code licensed under GPLv3</p>'
        )
      ),
    tags$script(src = "imgModal.js")
    )
  )

options(shiny.maxRequestSize = 60*1024^2)

# Define server 
server <- function(input, output, session) {
  
  source("server-inputData.R",local = TRUE)
  # 
  source("server-enrichGo.R",local = TRUE)
  # 
  source("server-plots.R",local = TRUE)
  
  
  GotoTab <- function(name){
    #updateTabItems(session, "tabs", name)
    
    shinyjs::show(selector = paste0("a[data-value=\"",name,"\"]"))
    
    shinyjs::runjs("window.scrollTo(0, 0)")
  }
}

shinyApp(ui, server)
