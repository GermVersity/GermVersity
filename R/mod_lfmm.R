#' lfmm UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_lfmm_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Latent factor mixed model analysis"),
                          tags$br(),
                          tags$p('Code to run a latent factor mixed model analysis (LFMM)
                                 using data from 259 accessions of wild Lima bean
                                 (Phaseolus lunatus L.), that belong to four wild gene
                                 pools (two Mesoamerican: MI and MII, and two Andean:AI
                                 and AII), that have been genotyped at 10668 SNP loci.
                                 LFMM applies a regression model to carry association
                                 tests among genetic variants and environmental variables.
                                 Correction for confounding effects, such as population
                                 structure, is done by including (unobserved) latent
                                 factors (set with K) which are estimated in parallel
                                 with the response variables (SNP loci) and the environmental
                                 variable of interest.'),
                          tags$br(),
                          shiny::fileInput(ns("filevcf"), "Choose a VCF file",
                                           multiple = F,
                                           accept = ".vcf",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filelfmm"), "Choose a LFMM file",
                                           multiple = F,
                                           accept = ".lfmm",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filetxt"), "Choose a TXT file (environment)",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filetxt1"), "Choose a TXT file (position)",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          tags$br()
      ),
      mainPanel(
        shiny::tabsetPanel(
          type = 'tabs',
          shiny::tabPanel('Descriptive',
                          tags$h2("Results"),
                          tags$br(),
                          tags$h3('Estimation of K (number of populations)'),
                          tags$br(),
                          plotOutput(ns('plot_data')),
                          tags$br(),
                          tags$p('The screeplot indicates that there are around K=6-8 main components in the data.'),
                          tags$br(),
                          plotOutput(ns('plot_data1'))),
          shiny::tabPanel('False Discovery Rate',
                          tags$br(),
                          tags$h3("Identify LFMM candidates using False Discovery Rate"),
                          tags$br(),
                          tags$p('The ridge_lfmm function returns an object that contains
                                  the following matrices:
                                  Matrix U: matrix of latent variables
                                  Matrix B: matrix of effects of all explanatory variables (loci)
                                  Matrix V: matrix of loadings for all latent variables'),
                          tags$br(),
                          tags$h3('Unadjusted p-values.'),
                          tags$br(),
                          plotOutput(ns('plot_data2')),
                          tags$br(),
                          tags$h3('GIF-adjusted p-values.'),
                          tags$br(),
                          plotOutput(ns('plot_data3')),
                          tags$br(),
                          tags$h3('Unadjusted p-values QQplots.'),
                          tags$br(),
                          plotOutput(ns('plot_data4')),
                          tags$br(),
                          tags$h3('GIF-adjusted p-values QQplots.'),
                          tags$br(),
                          plotOutput(ns('plot_data5'))),
          shiny::tabPanel('Combining results',
                          tags$br(),
                          tags$h3("Combining results with position matrices"),
                          tags$br(),
                          tags$h3('Unadjusted and Calibrated p-values.'),
                          tags$br(),
                          DT::dataTableOutput(ns('table1')),
                          tags$br(),
                          tags$h3('Manhattan Plot'),
                          tags$br(),
                          plotOutput(ns('plot_data6')),
                          tags$br(),
                          tags$h3('Effect size of the environmental variable of interest on each SNP'),
                          tags$br(),
                          plotOutput(ns('plot_data7')))
        )
      )
    )
  )
}

#' lfmm Server Functions
#'
#' @import qvalue
#' @import LEA
#' @import qqman
#' @import DT
#' @import lfmm
#'
#' @noRd
mod_lfmm_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    data <- reactive({

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      pc1 <- prcomp(Y1) # carry out a PCA

    })

    output$plot_data <- renderPlot({

      pc1 <- data()

      plot(pc1$sdev[1:15]^2, ylab = "percentage of variance explained", xlab = 'Axes', ) #plot results

    })

    output$plot_data1 <- renderPlot({
      pc1 <- data()

      screeplot(pc1, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")
    })

    output$plot_data2 <- renderPlot({
      tryCatch(
        {
          data1 <- LEA::vcf2lfmm(input$filevcf$datapath, force = TRUE)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      tryCatch(
        {
          bio10 <- read.table(input$filetxt$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                   X = bio10,
                                   K = 8)

      pv <- lfmm::lfmm_test(Y = Y1,
                            X = bio10,
                            lfmm = mod.lfmm,
                            calibrate = "gif")
      hist(pv$pvalue[,1], main="Unadjusted p-values", xlab = 'P value')
    })

    output$plot_data3 <- renderPlot({
      tryCatch(
        {
          data1 <- LEA::vcf2lfmm(input$filevcf$datapath, force = TRUE)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      tryCatch(
        {
          bio10 <- read.table(input$filetxt$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                   X = bio10,
                                   K = 8)

      pv <- lfmm::lfmm_test(Y = Y1,
                            X = bio10,
                            lfmm = mod.lfmm,
                            calibrate = "gif")

      hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values", xlab = 'P value calibrated')
    })

    output$plot_data4 <- renderPlot({
      tryCatch(
        {
          data1 <- LEA::vcf2lfmm(input$filevcf$datapath, force = TRUE)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      tryCatch(
        {
          bio10 <- read.table(input$filetxt$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )
      mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                   X = bio10,
                                   K = 8)

      pv <- lfmm::lfmm_test(Y = Y1,
                            X = bio10,
                            lfmm = mod.lfmm,
                            calibrate = "gif")

      pvaluesUncal <- pv$pvalue
      qqplot(rexp(length(pvaluesUncal), rate = log(10)),
             -log10(pvaluesUncal),
             xlab = "Expected quantile",
             pch = 19, cex = .4)
      abline(0,1)
    })

    output$plot_data5 <- renderPlot({
      tryCatch(
        {
          data1 <- LEA::vcf2lfmm(input$filevcf$datapath, force = TRUE)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      tryCatch(
        {
          bio10 <- read.table(input$filetxt$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )
      mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                   X = bio10,
                                   K = 8)
      pv <- lfmm::lfmm_test(Y = Y1,
                            X = bio10,
                            lfmm = mod.lfmm,
                            calibrate = "gif")

      pvalues <- pv$calibrated.pvalue
      qqplot(rexp(length(pvalues), rate = log(10)),
             -log10(pvalues), xlab = "Expected quantile",
             pch = 19, cex = .4)
      abline(0,1)
    })

    output$table1 <- DT::renderDataTable({
      tryCatch(
        {
          bio10 <- read.table(input$filetxt$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      tryCatch(
        {
          chr_pos_10668snps <- read.delim(input$filetxt1$datapath)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                   X = bio10,
                                   K = 8)

      pv <- lfmm::lfmm_test(Y = Y1,
                            X = bio10,
                            lfmm = mod.lfmm,
                            calibrate = "gif")
      qv <- qvalue::qvalue(pv$calibrated.pvalue)$qvalues

      lima_lfmm <- cbind(chr_pos_10668snps, pv$calibrated.pvalue, qv)
      colnames(lima_lfmm)[4:5] <- c("calibrated.pvalues", "pvalues")

      DT::datatable(lima_lfmm,
                    filter = 'top', extensions = c('Buttons', 'Scroller'),
                    options = list(scrollY = 650,
                                   scrollX = 500,
                                   deferRender = TRUE,
                                   scroller = TRUE,
                                   # paging = TRUE,
                                   # pageLength = 25,
                                   buttons = list('excel',
                                                  list(extend = 'colvis', targets = 0, visible = FALSE)),
                                   dom = 'lBfrtip',
                                   fixedColumns = TRUE),
                    rownames = FALSE)

    })

    output$plot_data6 <- renderPlot({
      tryCatch(
        {
          bio10 <- read.table(input$filetxt$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      tryCatch(
        {
          chr_pos_10668snps <- read.delim(input$filetxt1$datapath)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                   X = bio10,
                                   K = 8)

      pv <- lfmm::lfmm_test(Y = Y1,
                            X = bio10,
                            lfmm = mod.lfmm,
                            calibrate = "gif")
      qv <- qvalue::qvalue(pv$calibrated.pvalue)$qvalues

      lima_lfmm <- cbind(chr_pos_10668snps, pv$calibrated.pvalue, qv)
      colnames(lima_lfmm)[4:5] <- c("calibrated.pvalues", "pvalues")

      qqman::manhattan(lima_lfmm, chr="chrom", bp="bp", snp="snp",
                       p="calibrated.pvalues",suggestiveline = FALSE,
                       genomewideline = -log10(8e-05)) #genomewideline corresponds to qvalues < 0.1
    })

    output$plot_data7 <- renderPlot({
      tryCatch(
        {
          bio10 <- read.table(input$filetxt$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          Y1 <- read.table(input$filelfmm$datapath, quote="\"", comment.char="", na.strings="9")
        },
        error = function(e){
          stop('Upload LFMM file')
        }
      )

      tryCatch(
        {
          chr_pos_10668snps <- read.delim(input$filetxt1$datapath)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                   X = bio10,
                                   K = 8)

      b.values <- lfmm::effect_size(Y1, bio10, mod.lfmm)
      hist(b.values, main = '')
    })

  })
}

## To be copied in the UI
# mod_lfmm_ui("lfmm_1")

## To be copied in the server
# mod_lfmm_server("lfmm_1")
