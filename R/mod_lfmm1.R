#' lfmm1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_lfmm1_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Identify LFMM candidates using False Discovery Rate"),
                          tags$br(),
                          tags$p('The ridge_lfmm function returns an object that contains
                                  the following matrices:
                                  Matrix U: matrix of latent variables
                                  Matrix B: matrix of effects of all explanatory variables (loci)
                                  Matrix V: matrix of loadings for all latent variables'),
                          tags$br(),
                          shiny::fileInput(ns("filevcf"), "Choose a VCF file",
                                           multiple = F,
                                           accept = ".vcf",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filelfmm"), "Choose a LFMM file",
                                           multiple = F,
                                           accept = ".lfmm",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filetxt"), "Choose a TXT file",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          tags$br()
      ),
      mainPanel(tags$h2("Results"),
                tags$br(),
                tags$h3('Unadjusted p-values.'),
                tags$br(),
                plotOutput(ns('plot_data')),
                tags$br(),
                tags$h3('GIF-adjusted p-values.'),
                tags$br(),
                plotOutput(ns('plot_data1')),
                tags$br(),
                tags$h3('Unadjusted p-values QQplots.'),
                tags$br(),
                plotOutput(ns('plot_data2')),
                tags$br(),
                tags$h3('GIF-adjusted p-values QQplots.'),
                tags$br(),
                plotOutput(ns('plot_data3'))

      )
    )
  )
}

#' lfmm1 Server Functions
#'
#' @import LEA
#' @import lfmm
#'
#' @noRd
mod_lfmm1_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$plot_data <- renderPlot({
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

      output$plot_data1 <- renderPlot({
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

        pvaluesUncal <- pv$pvalue
        qqplot(rexp(length(pvaluesUncal), rate = log(10)),
               -log10(pvaluesUncal),
               xlab = "Expected quantile",
               pch = 19, cex = .4)
        abline(0,1)
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

        pvalues <- pv$calibrated.pvalue
        qqplot(rexp(length(pvalues), rate = log(10)),
               -log10(pvalues), xlab = "Expected quantile",
               pch = 19, cex = .4)
        abline(0,1)
      })

  })
}

## To be copied in the UI
# mod_lfmm1_ui("lfmm1_1")

## To be copied in the server
# mod_lfmm1_server("lfmm1_1")
