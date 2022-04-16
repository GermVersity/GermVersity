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
      shiny::sidebarPanel(tags$h3("atent factor mixed model analysis"),
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

                          shiny::fileInput(ns("filelfmm"), "Choose a LFMM file",
                                           multiple = F,
                                           accept = ".lfmm",
                                           buttonLabel = "Uploading..."),

                          tags$br(),
                          shiny::actionButton(ns('gob'), 'Click to continue..')
      ),
      mainPanel(tags$h2("Results"),
                tags$br(),
                tags$h3('Estimation of K (number of populations)'),
                tags$br(),
                plotOutput(ns('plot_data')),
                tags$br(),
                tags$p('The screeplot indicates that there are around K=6-8 main components in the data.'),
                tags$br(),
                plotOutput(ns('plot_data1'))

      )
    )
  )
}

#' lfmm Server Functions
#'
#' @import qvalue
#' @import LEA
#' @import qqman
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

  })
}

## To be copied in the UI
# mod_lfmm_ui("lfmm_1")

## To be copied in the server
# mod_lfmm_server("lfmm_1")
