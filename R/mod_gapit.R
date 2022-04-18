#' gapit UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gapit_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("GAPIT"),
                          tags$br(),
                          tags$p('GAPIT is a Genome Association and Prediction
                                 Integrated Tool freely available for Public
                                 since 2011 (Lipka et al. 2012). It has been
                                 updated frequently to incorporate the state
                                 of art methods for Genome Wide Association
                                 Study (GWAS). Currently, Dr. Jiabo Wang is
                                 leading the new development, GAPIT version
                                 3. In addition to the MLM methods, the new
                                 version implemented two new GWAS methods as
                                 FarmCPU (Liu et al. 2016) and BLINK (Huang
                                 et al. 2019). López-Hernández & Cortés (2019)
                                 explored this tool to obtain the Genome-Environment
                                 Associations (GEA) in crops (Common bean)
                                 using Genotyping by Sequencing (GBS) and
                                 Environmental Indices of Heat Stress.'),
                          tags$br(),
                          shiny::fileInput(ns("filexlsx"), "Choose a Excel (xlsx) file",
                                           multiple = F,
                                           accept = ".xlsx",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filetxt"), "Choose a TXT file",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          tags$br(),
                          shiny::actionButton(ns('gob'), 'Click to continue..')
      ),
      mainPanel(tags$h2("Results Downloading..."),
                tags$br(),
                tags$h3('See the files with the result in the next folder: ~/GAPIT/'),
                tags$br(),
                shiny::verbatimTextOutput(ns('ready')),
                tags$br(),
                tags$h3('Manhattan Circular plot for MLM'),
                shiny::actionButton(ns('Manhattan'), 'See plot..'),
                tags$br(),
                shiny::plotOutput(ns('plot_data1')),
                tags$h3('QQ plot for MLM'),
                shiny::actionButton(ns('QQ'), 'See plot..'),
                tags$br(),
                shiny::plotOutput(ns('plot_data2'))

      )
    )
  )
}

#' gapit Server Functions
#'
#' @import GAPIT3
#' @import bigmemory
#' @import genetics
#' @import gplots
#' @import LDheatmap
#' @import moments
#' @import readxl
#' @import readr
#'
#' @noRd
mod_gapit_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    values <- reactiveValues(ploting = NULL, plott = NULL)

    observeEvent(input$gob,{
      tryCatch(
        {
          Lima = readxl::read_excel(input$filexlsx$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          data = input$filetxt$datapath
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      myY = as.data.frame(Lima)
      myG <- as.data.frame(readr::read_delim(data,
                                      delim = "\t", escape_double = FALSE,
                                      col_names = FALSE, trim_ws = TRUE))

      PCA_number <- 8
      kinship <- "Vanraden"
      kinship_cluster <- "ward"
      MAF <- 0.05
      defaultW <- getOption("warn")
      dirname <-  '~/GAPIT/MLM'
      if (!dir.exists(dirname))dir.create(dirname,recursive=TRUE)

      # set working folder path for MLM
      setwd("~/GAPIT/MLM")
      # Run the MLM model
      myGAPIT_MLM <- GAPIT3::GAPIT(
        Y = myY,
        G = myG,
        model = "MLM",
        PCA.total = PCA_number,
        file.output = TRUE,
        Geno.View.output = FALSE,
        kinship.cluster = kinship_cluster,
        SNP.MAF = MAF
      )
      options(warn = defaultW)

      dirname <-  '~/GAPIT/Farm-CPU'
      if (!dir.exists(dirname))dir.create(dirname,recursive=TRUE)

      # set working folder path for Farm-CPU
      setwd('~/GAPIT/Farm-CPU')
      # Run the Farm-CPU model
      myGAPIT_MLM <- GAPIT3::GAPIT(
        Y=myY,
        G=myG,
        model="Farm-CPU",
        PCA.total=PCA_number,
        file.output=T,Geno.View.output=F,
        # kinship.algorithm=kinship,
        kinship.cluster=kinship_cluster,
        SNP.MAF=MAF
      )
      options(warn = defaultW)

      dirname <-  '~/GAPIT/Blink'
      if (!dir.exists(dirname))dir.create(dirname,recursive=TRUE)

      # set working folder path for Blink
      setwd('~/GAPIT/Blink')
      # Run the Blink model
      myGAPIT_MLM <- GAPIT3::GAPIT(
        Y=myY,
        G=myG,
        model="Blink",
        PCA.total=PCA_number,
        file.output=T,Geno.View.output=F,
        # kinship.algorithm=kinship,
        kinship.cluster=kinship_cluster,
        SNP.MAF=MAF
      )
      options(warn = defaultW)
      values$ploting <- 1

    })

    output$ready <- renderPrint({
      if (values$ploting == 1){
        print('Please click in See plot!')
      }else if(values$ploting == NULL){
        print('Please wait!')
      }
    })

    observeEvent(input$Manhattan, {
      if(values$ploting == 1){
        values$plott = 1
      }
    })

    output$plot_data1 <- renderPlot({
      if(values$plott == 1){
        image <- magick::image_read_pdf("~/GAPIT/MLM/GAPIT.Manhattan.Multiple.Plot.circular.pdf")
        file_i <- magick::image_scale(image, geometry = "700")
        plot(file_i)
      }
      else{
        file_i <- NULL
      }
    })

    output$plot_data2 <- renderPlot({
      if(values$plott == 1){
      image = magick::image_read_pdf("~/GAPIT/MLM/GAPIT.Multiple.QQ.plot.symphysic .pdf")
      file_i <- magick::image_scale(image, geometry = "700")
      plot(file_i)
    }else{file_i <- NULL}
    })


  })
}

## To be copied in the UI
# mod_gapit_ui("gapit_1")

## To be copied in the server
# mod_gapit_server("gapit_1")
