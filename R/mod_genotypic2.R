#' genotypic2 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_genotypic2_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("PCA"),
                          tags$br(),
                          tags$p('PCA (Principal component analysis) is a multivariate
                                  technique that summarizes
                                 the information provided by the genetic markers
                                 (for example, SNPs) into a few set of components.
                                 We can apply PCA to genlight objects with the
                                 function glPca of the package adgenet'),
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
                tags$br(),
                tags$br(),
                tags$h4("PCA plot"),
                shiny::plotOutput(ns("plotPCA"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  ))

      )
    )
  )
}

#' genotypic2 Server Functions
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
mod_genotypic2_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$plotPCA <- renderPlot({
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

      Genepool = as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      LimaBeanPCA <- adegenet::glPca(LimaBeanData3, nf = 3)
      # to carry out a PCA on a genlight object. With the argument nf as NULL,
      #you are asked interactively for the number of principal components to be
      #retained. For this data, three axes were retained.

      # to create a customized PCA plot with the package ggplot2

      Limapcascores <- as.data.frame(LimaBeanPCA$scores)
      Limapcascores$pop <- adegenet::pop(LimaBeanData3)


      set.seed(5)
      colors <- RColorBrewer::brewer.pal(n = adegenet::nPop(LimaBeanData3), name = "Set1")
      p <- ggplot2::ggplot(Limapcascores, ggplot2::aes(x=PC1, y=PC2, colour=pop))
      p <- p + ggplot2::geom_hline(yintercept = 0)
      p <- p + ggplot2::geom_vline(xintercept = 0)
      p <- p + ggplot2::geom_point(size=3)
      p <- p + ggplot2::theme_bw()
      p <- p + ggplot2::scale_color_manual(values=colors ) +
        ggplot2::xlab(sprintf("PC1 %f percent", 100*LimaBeanPCA$eig[1]/sum(LimaBeanPCA$eig))) +
        ggplot2::ylab(sprintf("PC2 %f percent", 100*LimaBeanPCA$eig[2]/sum(LimaBeanPCA$eig)))
      p

    })


  })
}

## To be copied in the UI
# mod_genotypic2_ui("genotypic2_1")

## To be copied in the server
# mod_genotypic2_server("genotypic2_1")
