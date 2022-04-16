#' genotypic5 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_genotypic5_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Basic diversity statistics per cluster (or population)"),
                          tags$br(),
                          tags$p('Once we have defined the number of clusters
                                 (or populations) in a sample, it is a good idea
                                 to estimate genetic diversity indexes for each
                                 cluster, for example if we want to find out
                                 which cluster is the most or least diverse.
                                 For doing this, we are going to use the program
                                 hierfstat to estimate observed heterozygosity
                                 (Ho) and mean gene diversities (Hs) within
                                 populations. A significant difference among
                                 Ho and Hs means that the population is not in
                                 HWE. Because different populations may not be
                                 in HWE for different reasons, diversity among
                                 populations may be hardly compared on the basis
                                 of Ho. This is because Hs is preferred to compare
                                 the genetic diversity between populations because
                                 with Hs all populations are assumed to be in HWE.
                                 It is also useful to estimate the fixation index
                                 Fis to measure the deviation to the assumption of
                                 HWE within each population.'),
                          tags$br(),
                          shiny::fileInput(ns("filevcf"), "Choose a VCF file",
                                           multiple = F,
                                           accept = ".vcf",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filetxt"), "Choose a TXT file",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading...")
      ),
      mainPanel(tags$h2("Results ML"),
                tags$h4("Summary Hs"),
                tags$p('To summarize Hs values within populations. We can observe
                       that populations are different in their genetic diversity,
                       with the highest values for the domesticated admixed,
                       followed by the Andean wild gene pool AI'),
                tags$h4("Summary Ho"),
                shiny::verbatimTextOutput(ns("summary1")),
                tags$p('To summarize Ho values within populations. We can
                       observe that within populations, Ho is very los, as
                       expected for a selfing species as Lima bean.'),
                tags$h4("Summary Fis"),
                shiny::verbatimTextOutput(ns("summary2")),
                tags$p('To summarize Fis values within populations.
                       Fis is calculated as (Hs-Ho)/Hs. We can see
                       that Fis values are positive for all populations
                       and in some of them close to 1, as expected for
                       a predominantly autogamous species as Lima bean.'),
                shiny::verbatimTextOutput(ns("summary3")),
                tags$p('A useful plot to compare genetic diversity among genepools'),
                shiny::plotOutput(ns('boxplot'))

      )
    )
  )
}

#' genotypic5 Server Functions
#'
#' @import adegenet
#' @import hierfstat
#' @import ape
#' @import poppr
#' @import pegas
#' @import vcfR
#' @import RColorBrewer
#' @import ggplot2
#' @import uwIntroStat
#'
#' @noRd
mod_genotypic5_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    data <- reactive({
      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop("Upload VCF file")
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop("Upload TXT file")
        }
      )
      Genepool <- as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                  pop= Genepool,
                                  NA.char= "NA")
    })

    output$summary1 <- renderPrint({

      LimaBeanData2 = data()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      print(summary(diversity.clusters$Hs))
    })

    output$summary2 <- renderPrint({
      LimaBeanData2 = data()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      print(summary(diversity.clusters$Ho))
    })

    output$summary3 <- renderPrint({
      LimaBeanData2 = data()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      print(summary(diversity.clusters$Fis))
    })

    output$boxplot <- renderPlot({
      LimaBeanData2 = data()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      boxplot(diversity.clusters$Hs, ylab="Hs")
    })

  })
}

## To be copied in the UI
# mod_genotypic5_ui("genotypic5_1")

## To be copied in the server
# mod_genotypic5_server("genotypic5_1")
