#' SDM UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_SDM_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Spatial Distribution Modeling (SDM)"),
                          tags$br(),
                          tags$p("Exploration of the spatial distribution modeling
                                 using the Maxent Algorithm (Phillips et al. 2006,
                                 Phillips et al. 2008, Phillips et al. 2017)
                                 by means of the R-package ENMeval (Kass et al.
                                 2021, https://doi.org/10.1111/2041-210X.13628)"),
                          tags$br(),
                          shiny::fileInput(ns("filexlsx"), "Choose a Excel file",
                                           multiple = F,
                                           accept = ".xlsx",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filebio"), "Choose the Layers files",
                                           multiple = T,
                                           accept = c(".bil",".hdr"),
                                           buttonLabel = "Uploading..."),
                          tags$hr(),
                          tags$p('Use WorldClim data (2022)'),
                          shiny::actionButton('WC', 'Click me!')
      ),
      mainPanel(tags$h2("Results"),
                tags$br(),
                tags$h4("Responses of the variables (predictors)"),
                shiny::plotOutput(ns("plot_data"), height = "500px")
                #tags$h4("Plotting the prediction modeling using the ENMevaluate results"),
                #plotOutput(ns("plot_data2"), height = "500px"),
                #tags$h4("Plotting the prediction modeling using the GLM results"),
                #plotOutput(ns("plot_data2"), height = "500px")

      )
    )
  )
}

#' SDM Server Functions
#' @import magrittr
#' @import maptools
#' @import sp
#' @import raster
#' Change library
#' @import geodata
#' @import sdm
#' @import sf
#' @import rworldxtra
#' Maxent Algorithm in R
#' @import maxnet
#' @import ENMeval
#' Biological databases
#' @import BIEN
#' @import rgbif
#' Data Manipulation
#' @import readxl
#' @import tidyverse
#' Machine Learning
#' @import caret
#' Plotting
#' @import ggplot2
#' @import ggspatial
#' @import viridis
#' @import ecospat
#' @noRd
mod_SDM_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$plot_data <- renderPlot({
      tryCatch(
        {
          Specie = readxl::read_xlsx(input$filexlsx$datapath)
        },
        error = function(e){
          stop('Upload Excel file')
        }
      )

      # Deleting duplicates
      recordsSpecie <- unique(Specie)
      # Deleting occurrence points with NA
      recordsSpecie <- na.omit(Specie)
      # Deleting ID column
      recordsSpecie <- recordsSpecie[,-1]

      # We transform into sf format
      Lunatus <- recordsSpecie |> sf::st_as_sf(coords = c(1, 2),
                                            crs = "+proj=longlat +ellps=WGS84
                                      +datum=WGS84 +no_defs +towgs84=0,0,0")
      # Generate a buffer
      Hull <- Lunatus |>
        sf::st_union() |>
        sf::st_convex_hull()
      Buffer <- Hull |>
        sf::st_buffer(dist = 1) |>
        sf::st_as_sf()

      Bioclimatic <- raster::getData("worldclim", res = 2.5, var = "bio", path = tempdir())

      # Crop layers using the buffer
      Bioclimatic <- Bioclimatic |>
        raster::crop(Buffer) |>
        raster::trim()

      # Selección del número de background points
      Number_background_points = 5000
      # Run ENMevaluate
      Results <- ENMeval::ENMevaluate(occs =  recordsSpecie, envs = Bioclimatic,
                             n.bg = Number_background_points,
                             algorithm = 'maxnet', partitions = 'block',
                             tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:2))
      # Modeling results
      Results@results


      ## Best Model Prediction
      Models <- Results@results
      Models$ID <- 1:nrow(Models)
      Models <- Models |>
        dplyr::arrange(AICc)
      BestModels <- Results@models[[Models$ID[1]]]
      Prediction <- raster::predict(Bioclimatic, BestModels, type = "cloglog")

      plot(BestModels, c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10",
                         "bio11","bio12"), type = "cloglog")

    })

  })
}

## To be copied in the UI
# mod_SDM_ui("SDM_1")

## To be copied in the server
# mod_SDM_server("SDM_1")
