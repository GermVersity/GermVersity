#' genotypic UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList

mod_genotypic_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Germplasm bank management"),
                          tags$br(),
                          tags$p('In this demo, we are going to show how
                          to conduct basic genetic diversity
                          analyses of SNP data in a sample of
                          genebank accessions of Lima bean to
                          explore the genetic structure of the
                          sample. In this demo we are going to
                          learn how to calculate distance matrices
                          among individuals and populations,
                          visualize the distance matrices using
                          clustering algorithms (UPGMA and NJ),
                          how to conduct a principal component
                          analyses and a discriminant analysis
                          of principal components to assign
                          individuals into a k number of populations.'),
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
      mainPanel(tags$h2("Results traditional methods"),
                tags$br(),
                tags$h4('Let´s calculate a Euclidean distance matrix
                        between individuals on the basis of the
                        observed allele frequencies within individuals.
                        This distance is not a genetic
                        distance but a geometric distance since it
                        does not assume any evolutionary model.'),
                tags$br(),
                tags$h4("Circular dendrogram"),
                shiny::plotOutput(ns("plot_data"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  )),
                tags$h4("Horizontal dendrogram"),
                plotOutput(ns("plot_data2"), height = "500px",
                           dblclick = "double_click2",
                           brush = brushOpts(
                             id = "brush_plot",
                             resetOnNew = TRUE
                           ))

      )
    )
  )
}

#' genotypic Server Functions
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
mod_genotypic_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    options(shiny.maxRequestSize=30*1024^2)
    ranges <- reactiveValues(x = NULL, y = NULL)

    observeEvent(input$double_click, {
      brush <- input$brush_plot
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })

    output$plot_data = renderPlot({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      Genepool <- as.character(data$Genepool)
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                   pop= Genepool,
                                   NA.char= "NA")

      EuclideanDistance <- dist(LimaBeanData2,
                                method = "euclidean",
                                diag = FALSE,
                                upper = FALSE,
                                p=2)
      NJtree <- ape::nj(EuclideanDistance)

      mycol = c("light blue",
                "gray",
                "green",
                "red",
                "blue",
                "light green",
                "pink")[LimaBeanData2$pop]

      plot(NJtree, tip.color=mycol, type="fan",
           x.lim = ranges$x,
           y.lim = ranges$y)
    })

    ranges1 <- reactiveValues(x = NULL, y = NULL)

    output$plot_data2 = renderPlot({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      Genepool <- as.character(data$Genepool)
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                   pop= Genepool,
                                   NA.char= "NA")

      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS) # convert the vcf file into a genlight object
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      UPGMAtree <- poppr::aboot(LimaBeanData3, tree = "upgma", distance = poppr::bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

      mycol = c("light blue",
                "gray",
                "green",
                "red",
                "blue",
                "light green",
                "pink")[LimaBeanData2$pop]

      plot(UPGMAtree, tip.color=mycol, cex = 1.2,
           x.lim = ranges1$x,
           y.lim = ranges1$y)
    })

    observeEvent(input$double_click2, {
      brush <- input$brush_plot
      if (!is.null(brush)) {
        ranges1$x <- c(brush$xmin, brush$xmax)
        ranges1$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges1$x <- NULL
        ranges1$y <- NULL
      }
    })
  })
}

## To be copied in the UI
# mod_genotypic_ui("genotypic_1")

## To be copied in the server
# mod_genotypic_server("genotypic_1")
