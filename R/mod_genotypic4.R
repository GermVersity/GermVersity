#' genotypic4 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_genotypic4_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Basic diverity statistics per locus"),
                          tags$br(),
                          tags$p('In this section we are going to estimate
                                 basic diversity statistics per locus in the
                                 whole sample. We are going to use the program
                                 adegenet to estimate the number of alleles
                                 per locus (NA), observed heterosygosity per
                                 locus (Hobs) and expected (Hexp) heterozygosity
                                 per locus. For SNP makers, NA might not be
                                 very useful because these loci are expected
                                 to be biallelic in populations (according to
                                 the infinite site mutational model) and we
                                 have also filtered the vcf file to include
                                 only biallelic SNPs. Hobs is the proportion
                                 of heterozygous individuals that were observed
                                 in the locus. Hexp is the heterozygosity we
                                 expect to observe in the locus assuming that
                                 the population (in this case the whole sample)
                                 is in Hardy-Weinberg equilibrium (HWE)in that
                                 locus. A significant difference among Hobs and
                                 Hexp means that the locus is not in HWE'),
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
                tags$h4("Hobs per locus"),
                shiny::plotOutput(ns("plot_d1"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  )),
                tags$h4("Hexp per locus"),
                shiny::plotOutput(ns("plot_d2"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  )),
                tags$h4("Hexp as a function of Hobs per locus"),
                shiny::plotOutput(ns("plot_d3"), height = "500px",
                                  dblclick = "double_click",
                                  brush = brushOpts(
                                    id = "brush_plot",
                                    resetOnNew = TRUE
                                  )),
                tags$h4("Testing the difference among Hexp and Hobs per locus"),
                tags$p('We can apply Bartlett tests to assess whether Hexp is
                       different from Hobs. As the number of loci, and
                       therefore of multiple tests, is high, you may want
                       to correct the p-value accordingly (for example with
                       the Bonferroni correction).'),
                shiny::verbatimTextOutput(ns("barlett"))

      )
    )
  )
}

#' genotypic4 Server Functions
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
mod_genotypic4_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$plot_d1 <- renderPlot({

      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file') #safeError(e)
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )
      Genepool <- as.character(data$Genepool)

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char = "NA")

      diversity1 <- base::summary(LimaBeanData2)#we use the genind object created above
      print(diversity1)#['Hobs'][[1]])
      plot(diversity1['Hobs'][[1]],
           xlab="loci number",
           ylab="Hobs",
           ylim = c(0,0.9))
      #As expected for a selfer as Lima bean, in most loci the number of
      #heterosygous individuals is low

    })

    output$plot_d2 <- renderPlot({

      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file') #safeError(e)
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )
      Genepool <- as.character(data$Genepool)

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char = "NA")
      diversity <- base::summary(LimaBeanData2)#we use the genind object created above
      plot(diversity['Hexp'][[1]],
           xlab="loci number",
           ylab="Hexp",
           ylim = c(0,0.6))
      #A great variation in Hexp is observed. The maximum value for Hexp is 0.5
      #and it is observed when in a locus both alleles have the same frequency of 0.5.

    })

    output$plot_d3 <- renderPlot({

      diversity <- base::summary(LimaBeanData2)#we use the genind object created above

      plot(diversity['Hobs'][[1]], diversity['Hexp'][[1]], xlab="Hobs", ylab="Hexp")
      #This plot shows how Hexp is related to Hobs. For loci in HWE,
      #Hexp = Hobs. In the plot we observe that in most loci, Hexp is higher
      #than Hobs, as expected for selfing species as Lima bean.

    })

    output$barlett <- renderText({
      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file') #safeError(e)
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )
      Genepool <- as.character(data$Genepool)

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char = "NA")
      diversity <- base::summary(LimaBeanData2)#we use the genind object created above

      bartlett.test(list(diversity['Hexp'][[1]], diversity['Hobs'][[1]]))
    })

  })
}

## To be copied in the UI
# mod_genotypic4_ui("genotypic4_1")

## To be copied in the server
# mod_genotypic4_server("genotypic4_1")
