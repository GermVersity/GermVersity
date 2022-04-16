#' genotypic3 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_genotypic3_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("DPCA"),
                          tags$br(),
                          tags$p('To identify genetic clusters, we can apply
                                 another multivariate approach known as DAPC.
                                 This approach is convenient when we are more
                                 interested in describing the diversity
                                 among groups of individuals than within groups.
                                 This approach is focused on finding discriminant
                                 functions that better describe the differences
                                 among groups, while minimizing the differences
                                 within groups. However, to find the discriminant
                                 functions, DAPC needs the groups to be known a
                                 priori and in many cases we just do not know
                                 how many groups are present in our sample.
                                 To address this issue, the adegenet package
                                 implements the function find.clusters (to find
                                 the number of clusters by running first a
                                 PCA and then using the k-means algorithm)
                                 and the function dapc to establish how are
                                 the relationships among the clusters.'),
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
                tags$h4("DPCA plot"),
                shiny::plotOutput(ns("plotDPCA"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  )),
                tags$h4("Clusters"),
                shiny::dataTableOutput(ns("table1")
                ),
                tags$h4("DPCA custom plot"),
                shiny::plotOutput(ns("plotDPCA2"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  )),
                tags$h4("Membership probabilities"),
                shiny::plotOutput(ns("plotDPCA3"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  ))

      )
    )
  )
}

#' genotypic3 Server Functions
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
mod_genotypic3_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    data <- reactive({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop(safeError(e))
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop(safeError(e))
        }
      )

      Genepool = as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                  pop= Genepool,
                                  NA.char= "NA")

      grp <- adegenet::find.clusters(LimaBeanData2,
                           max.n.clust = 30,
                           choose.n.clust = FALSE,
                           n.clust=5,
                           n.pca=100, choose=FALSE)
      # You are asked interactively for the number of PCs to be
      #retained (we chose 500, more than the maximum for this dataset)
      #and the number of clusters to be retained (we chose K=5 according
      #to the BIC values shown in the graph and previous biological knowledge).
      #Choosing the "right" number of clusters is a complex issue, the BIC
      #graph is just a guide. The researcher should explore different numbers
      #of clusters and choose the number that makes more sense from a biological point of view.

      dapc <- adegenet::dapc(LimaBeanData2,
                   grp$grp,
                   n.pca=100, n.da=4)
    })

    output$plotDPCA <- renderPlot({

      dapc <- data()
      ade4::scatter(dapc) # to obtain a basic scatterplot of the dapc

    })


    output$plotDPCA2 <- renderPlot({

      dapc <- data()

      myCol2 <- c("pink","red","blue","light blue", "green") # to assign colors to each of the five clusters

      ade4::scatter(dapc, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol2, solid=1.0,
                    cex=3,clab=0, leg=TRUE, posi.leg= "bottomleft", scree.pca=TRUE, posi.pca = "topright", ratio.pca=0.3)

    })

    output$table1 <- renderDataTable({
      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop(safeError(e))
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop(safeError(e))
        }
      )

      Genepool = as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                  pop= Genepool,
                                  NA.char= "NA")

      grp <- adegenet::find.clusters(LimaBeanData2,
                           max.n.clust = 30,
                           choose.n.clust = FALSE,
                           n.clust=5,
                           n.pca=100, choose=FALSE)

      table(adegenet::pop(LimaBeanData2), grp$grp)

    })

    output$plotDPCA3 <- renderPlot({
      dapc <- data()

      adegenet::compoplot.dapc(dapc) #to visualize the membership probabilities as a bar plot

    })

  })
}

## To be copied in the UI
# mod_genotypic3_ui("genotypic3_1")

## To be copied in the server
# mod_genotypic3_server("genotypic3_1")
