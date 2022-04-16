#' genotypic1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_genotypic1_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Distance Matrix"),
                          tags$br(),
                          tags$p("Genetic distances measure how genetically
                                 similar are individuals or populations and
                                 are based on a specific evolutionary model.
                                 Genetic distances can be seen as summary
                                 statistics because they take the whole data
                                 set (for example, data from SNP loci) and
                                 summarize the genetic differentiation between
                                 samples (individuals or populations) in a
                                 value. One of the most used genetic distances
                                 is Nei's genetic distance which measures
                                 the genetic distance among populations on
                                 the basis of their genetic identity, namely
                                 the proportion of alleles that are shared
                                 between populations. For our data set, we
                                 have defined gene pools as populations in the
                                 genind object called LimaBeanData2.
                                 To calculate Nei's genetic distances among
                                 populations we need to add the strata to
                                 the genind object (in this case on the basis
                                 of the gene pool) to define the populations.
                                 Finally, we will build with this distance
                                 matrix a UPGMA tree and a neighbor-joining (
                                 NJ) topology, applying 1000 bootstrap permutations
                                 to get statistical support of the groups in the topologies."),
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
      mainPanel(
        tags$br(),
        tags$h4("UPGMA topology"),
        shiny::plotOutput(ns("distance_matrix1"), height = "500px",
                          dblclick = "double_click",
                          brush = brushOpts(
                            id = "brush_plot",
                            resetOnNew = TRUE
                          )),
        tags$h4("NJ topology"),
        plotOutput(ns("distance_matrix2"), height = "500px",
                   dblclick = "double_click2",
                   brush = brushOpts(
                     id = "brush_plot",
                     resetOnNew = TRUE
                   ))

      )
    )
  )
}

#' genotypic1 Server Functions
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
mod_genotypic1_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$distance_matrix1 <- renderPlot({
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
      Genepool <- as.character(data$Genepool)
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char= "NA")

      adegenet::strata(LimaBeanData2) <- data.frame(adegenet::pop(LimaBeanData2)) #to add the strata of a genind object
      adegenet::nameStrata(LimaBeanData2) <- ~Genepool # to assign a name to the strata of a genind object
      #Building of a UPGMA topology with bootstrap support using the function aboot of the poppr package
      set.seed(999)
      poppr::aboot(LimaBeanData2, strata = Genepool, tree= "upgma", distance="nei.dist", sample=1000, cutoff=0, showtree=TRUE,quiet=FALSE)

    })

    output$distance_matrix2 <- renderPlot({
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
      Genepool <- as.character(data$Genepool)
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")

      adegenet::strata(LimaBeanData2) <- data.frame(adegenet::pop(LimaBeanData2)) #to add the strata of a genind object
      adegenet::nameStrata(LimaBeanData2) <- ~Genepool # to assign a name to the strata of a genind object
      #Building of a UPGMA topology with bootstrap support using the function aboot of the poppr package
      set.seed(999)
      #Building of a NJ topology with bootstrap support using the function aboot of the poppr package
      poppr::aboot(LimaBeanData2, strata = Genepool, tree= "nj", distance="nei.dist", sample=1000, cutoff=0, showtree=TRUE,quiet=FALSE)
    })

  })
}

## To be copied in the UI
# mod_genotypic1_ui("genotypic1_1")

## To be copied in the server
# mod_genotypic1_server("genotypic1_1")
