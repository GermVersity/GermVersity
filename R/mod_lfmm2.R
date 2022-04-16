#' lfmm2 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_lfmm2_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Combining results with position matrices"),
                          tags$br(),
                          tags$p(''),
                          tags$br(),
                          shiny::fileInput(ns("filetxt"), "Choose a TXT file",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filelfmm"), "Choose a LFMM file",
                                           multiple = F,
                                           accept = ".lfmm",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filetxt1"), "Choose a TXT file",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          tags$br()
      ),
      mainPanel(tags$h2("Results"),
                tags$br(),
                tags$h3('Unadjusted and Calibrated p-values.'),
                tags$br(),
                DT::dataTableOutput(ns('table1')),
                tags$br(),
                tags$h3('Manhattan Plot'),
                tags$br(),
                plotOutput(ns('plot_data1')),
                tags$br(),
                tags$h3('Effect size of the environmental variable of interest on each SNP'),
                tags$br(),
                plotOutput(ns('plot_data2'))

      )
    )
  )
}

#' lfmm2 Server Functions
#'
#' @import LEA
#' @import lfmm
#' @import DT
#'
#' @noRd
mod_lfmm2_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

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
      qv <- qvalue(pv$calibrated.pvalue)$qvalues

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

    output$plot_data1 <- renderPlot({
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
      qv <- qvalue(pv$calibrated.pvalue)$qvalues

      lima_lfmm <- cbind(chr_pos_10668snps, pv$calibrated.pvalue, qv)
      colnames(lima_lfmm)[4:5] <- c("calibrated.pvalues", "pvalues")

      qqman::manhattan(lima_lfmm, chr="chrom", bp="bp", snp="snp",
                p="calibrated.pvalues",suggestiveline = FALSE,
                genomewideline = -log10(8e-05)) #genomewideline corresponds to qvalues < 0.1
    })

    output$plot_data2 <- renderPlot({
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
# mod_lfmm2_ui("lfmm2_1")

## To be copied in the server
# mod_lfmm2_server("lfmm2_1")
