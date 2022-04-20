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
                shiny::plotOutput(ns("plot_data"), height = "500px"),
                tags$h4("Plotting the prediction modeling using the ENMevaluate results"),
                shiny::plotOutput(ns("plot_data2"), height = "500px"),
                tags$h4("Plotting the prediction modeling using the GLM results"),
                plotOutput(ns("plot_data3"), height = "500px")

      )
    )
  )
}


#' SDM Server Functions
#'
#' @import magrittr
#' @import maptools
#' @import sp
#' @import raster
#' @import geodata
#' @import sdm
#' @import sf
#' @import rworldxtra
#' @import maxnet
#' @import ENMeval
#' @import BIEN
#' @import rgbif
#' @import readxl
#' @import tidyverse
#' @import caret
#' @import ggplot2
#' @import ggspatial
#' @import viridis
#' @import ecospat
#'
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
      Number_background_points = 1000
      # Run ENMevaluate
      Results <- ENMeval::ENMevaluate(occs =  recordsSpecie, envs = Bioclimatic,
                                      n.bg = Number_background_points,
                                      algorithm = 'maxnet', partitions = 'block',
                                      tune.args = list(fc = c("L","LQ"), rm = 1:2)) #,"LQH","H"

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

    output$plot_data2 <- renderPlot({
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
      Number_background_points = 1000
      # Run ENMevaluate
      Results <- ENMeval::ENMevaluate(occs =  recordsSpecie, envs = Bioclimatic,
                                      n.bg = Number_background_points,
                                      algorithm = 'maxnet', partitions = 'block',
                                      tune.args = list(fc = c("L","LQ"), rm = 1:2)) #,"LQH","H"
      # Modeling results
      Results@results

      ## Best Model Prediction
      Models <- Results@results
      Models$ID <- 1:nrow(Models)
      Models <- Models |>
        dplyr::arrange(AICc)
      BestModels <- Results@models[[Models$ID[1]]]
      Prediction <- raster::predict(Bioclimatic, BestModels, type = "cloglog")

      gdd_ggplot2 <- Prediction$layer |> as("SpatialPixelsDataFrame") |>
        as.data.frame()

      # Plotting using ggplot2
      plot_1 <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = gdd_ggplot2,
                           ggplot2::aes(x = x, y = y, fill = layer)) +
        viridis::scale_fill_viridis(option="turbo") +
        ggplot2::xlab("Longitude") +
        ggplot2::ylab("Latitude")+
        ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.5,
                                                                linetype = 'solid',
                                                                colour = "gray14"),
              panel.grid.minor = ggplot2::element_line(size = 0.25,
                                                       linetype = 'solid',
                                                       colour = "gray14"),
              panel.background = ggplot2::element_rect(fill = "black",
                                                       colour = "black",
                                                       size = 0.5,
                                                       linetype = "solid"))+
        ggplot2::labs(fill = "Probability") +
        ggplot2::labs(title= "Spatial Distribution Modeling") +
        ggplot2::labs(subtitle = "MaxEnt By ENMeval (maxnet)",
                      caption = " Maxent Algorithm (Phillips et al. 2006, Phillips et al. 2008, Phillips et al. 2017)\n using the R-package ENMeval (Kass et al. 2021, https://doi.org/10.1111/2041-210X.13628)")+
        ggplot2::theme(plot.caption.position = "plot",
                       plot.caption = ggplot2::element_text(hjust = 0)) +
        ggplot2::labs(fill = "Probability")+
        ggspatial::annotation_north_arrow(location = "tr",
                                          which_north = "true",
                                          pad_x = ggplot2::unit(0.4, "in"),
                                          pad_y = ggplot2::unit(0.4, "in"),
                                          style = ggspatial::north_arrow_nautical(fill = c("grey40", "gray49"),
                                                                                  line_col = "gray30"))
      plot_1
    })

    output$plot_data3 <- renderPlot({
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

      # Option 1: Extracting Environmental Data from Occurrences Points of Species
      Occurrences<-raster::extract(Bioclimatic, recordsSpecie)
      # Removing Missing Data
      Occurrences<-na.omit(Occurrences)
      # Removing Repeating Data from Points
      Occurrences<-unique(Occurrences)

      # Selecting pseudo-absence Points
      background_data <- raster::sampleRandom(Bioclimatic,length(recordsSpecie$Longitude),cells=T)
      # Removing Repeating Data from Layer Cells
      background_data <- unique(background_data)
      # Removing the Number of Cells Column
      background_data <- background_data[,-1]
      # Merging Data
      Environmental_values <- data.frame(rbind(Occurrences,background_data))
      # Labeling the presences as '1' and background as '0'
      y <- c(rep(1,nrow(Occurrences)), rep(0,nrow(background_data)))
      Data_01 <- cbind(y, Environmental_values)
      # Table as Data.Frame
      Data_SDM <- as.data.frame(Data_01)
      Data_SDM <- na.omit(Data_SDM)

      # Set Seed
      set.seed(2021)
      training <- caret::createDataPartition(y = Data_SDM$y, p = 0.8, list = FALSE, times = 1)

      # Training
      Data_training <- Data_SDM[training, ]
      Data_training_modeling <- Data_training
      Data_training_modeling$y <- as.factor(Data_training_modeling$y)
      levels(Data_training_modeling$y) <- c("back", "pres")
      # Training Data Ratio
      prop.table(table(Data_training$y))
      #Testing
      Data_test  <- Data_SDM[-training, ]
      Data_test_modeling <- Data_test
      Data_test_modeling$y <- as.factor(Data_test_modeling$y)
      levels(Data_test_modeling$y) <- c("back", "pres")
      # Test Data Ratio
      prop.table(table(Data_test$y))

      # Definition of partitions and repetitions
      partitions  <- 5
      repetitions <- 3

      # Definition of Training
      control_tr <- caret::trainControl(method = "repeatedcv", number = partitions,
                                 repeats = repetitions,
                                 returnResamp = "final", verboseIter = FALSE,
                                 classProbs = TRUE)
      # Model fit
      SDM_glm_model <- caret::train(y ~ ., data = Data_training_modeling,
                             method = "glm",
                             preProcess = c('center', 'scale'),
                             trControl = control_tr,
                             metric = 'Accuracy')
      # Model output
      #SDM_glm_model

      # Testing
      Prediction_raw <- raster::predict(SDM_glm_model, newdata = Data_test_modeling)
      # Confusio nMatrix
      confusionMatrix <- caret::confusionMatrix(data = Prediction_raw, reference = Data_test_modeling$y,positive = "pres")
      #confusionMatrix

      # Error Test
      error_test <- mean(Prediction_raw != Data_test_modeling$y)

      # Extrapolating Model to others Zones
      mapGLM <- raster::predict(Bioclimatic,SDM_glm_model, progress="text",type = "prob",args=c("extrapolate=T", "doclamp=TRUE"))
      mapGLM <- 1-mapGLM

      # Projections to Data.Frame Format
      mapGLM.ggplot2 <- mapGLM |>
        as("SpatialPixelsDataFrame") |>
        as.data.frame()

      # Plotting using ggplot2
      plot_2 <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = mapGLM.ggplot2,
                           ggplot2::aes(x = x, y = y, fill = layer)) +
        viridis::scale_fill_viridis(option="turbo") +
        ggplot2::xlab("Longitude") +
        ggplot2::ylab("Latitude")+
        # annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.5,
                                                                linetype = 'solid',
                                                                colour = "gray14"),
              panel.grid.minor = ggplot2::element_line(size = 0.25,
                                                       linetype = 'solid',
                                                       colour = "gray14"),
              panel.background = ggplot2::element_rect(fill = "black",
                                                       colour = "black",
                                                       size = 0.5,
                                                       linetype = "solid"))+
        ggplot2::labs(fill = "Probability") +
        ggplot2::labs(title= "Spatial Distribution Modeling") +
        ggplot2::labs(subtitle = "Generalized Linear Models By glm() function",
                      caption = " glm() function from the basic R-package using caret (https://www.cienciadedatos.net/machine-learning-r.h)")+
        ggplot2::theme(plot.caption.position = "plot",
                       plot.caption = ggplot2::element_text(hjust = 0)) +
        ggplot2::labs(fill = "Probability") +
        ggspatial::annotation_north_arrow(location = "tr",
                                          which_north = "true",
                                          pad_x = ggplot2::unit(0.4, "in"),
                                          pad_y = ggplot2::unit(0.4, "in"),
                                          style = ggspatial::north_arrow_nautical(fill = c("grey40", "gray49"),
                                                                                  line_col = "gray30"))
      plot_2
    })

  })
}

## To be copied in the UI
# mod_SDM_ui("SDM_1")

## To be copied in the server
# mod_SDM_server("SDM_1")
