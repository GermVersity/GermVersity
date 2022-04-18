#' descriptor UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_descriptor_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(tags$h3("Oescriptor optimization"),
                   tags$p('Distribution of minimum depth and number of trees for
                          each of the descriptors evaluated. It is observed in
                          order of importance of the descriptor, the number of
                          nodes associated with each descriptor and the interaction
                          between them.'),
                   shiny::fileInput(ns("filecsv"), "Choose a CVS file",
                                    multiple = F,
                                    accept = c("text/csv",
                                               "text/comma-separated-values, text/plain",
                                               ".csv"),
                                    buttonLabel = "Uploading...")
      ),
      mainPanel(
        tags$h3("ML Results"),
        plotOutput(ns('min_depth')),
        tags$hr(),
        plotOutput(ns('multi')),
        tags$hr(),
        plotOutput(ns('inter'))
      )
    )
  )
}

#' descriptor Server Functions
#'
#' @import readr
#' @import readxl
#' @import knitr
#' @import randomForest
#' @import randomForestExplainer
#' @import dplyr
#' @import ggplot2
#' @import DT
#' @import caret
#'
#' @noRd
mod_descriptor_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    model <- reactive({
      tryCatch(
        {
          # Import data from module import_data
          racesR <- read.csv(input$filecsv$datapath)
        },
        error = function(e){
          stop('Upload CSV file')
        }
      )

      ## Name data ####

      data <- racesR

      ### Convert double data to factor ###
      indx <- sapply(data, is.double)
      data[indx] <- lapply(data[indx], function(data) as.factor(as.double(data)))

      ### Convert factor data to double ###

      data$P100S <- as.numeric(data$P100S)
      data$DAF <- as.numeric(data$DAF)
      data$DAM <- as.numeric(data$DAM)
      data$LONG_SEM <- as.numeric(data$LONG_SEM)
      data$ANCHO_SEM <- as.numeric(data$ANCHO_SEM)
      data$GRUESO_SEM <- as.numeric(data$GRUESO_SEM)
      data$Races <- as.factor(data$Races)

      ######### ML ######################

      Datos <- as.data.frame(data[,-1])
      missingcols <- sapply(Datos, function(x) { any(is.na(x)) })
      tcontrol <- caret::trainControl(method="repeatedcv", number=50, repeats=3)
      set.seed(123)

      # replace data by keeping only those variables that don't have missing data
      data <- Datos[ , !missingcols]

      data %>% dplyr::group_by(Races) %>% summarise(n = n())

      # create training and test sets
      inTrain <- caret::createDataPartition(y = data$Races, p = 0.6, list = FALSE)

      # subset
      training <- data[inTrain, ]
      testing <- data[-inTrain, ]

      training$Races <- as.factor(training$Races)
      testing$Races <- as.factor(testing$Races)

      metric <- "Accuracy"
      mtry <- sqrt(ncol(training))
      modFitN <- caret::train(Races~., method = "rf",
                              data = training[], trControl = tcontrol,
                              metric=metric, tuneLength=15)

      predictions.rfN<- stats::predict(modFitN, newdata = testing[])

      Datos$Races<- as.factor(Datos$Races)

      modFit.rf <- randomForest::randomForest(Races ~ ., data = Datos, localImp = TRUE, mtry=modFitN$bestTune$mtry)

      VI <- varImp(modFit.rf, scale=T)
      VI$Variable <- row.names(VI)
      return(modFit.rf)
    })

    output$min_depth <- renderPlot({
      modFit.rf <- model()
      min_depth_frame <- randomForestExplainer::min_depth_distribution(modFit.rf)
      randomForestExplainer::plot_min_depth_distribution(min_depth_frame, mean_sample = "top_trees")
    })

    output$multi <- renderPlot({
      modFit.rf <- model()
      randomForestExplainer::plot_multi_way_importance(modFit.rf, size_measure = "no_of_nodes")
    })

    output$inter <- renderPlot({
      modFit.rf <- model()
      importance_frame <- randomForestExplainer::measure_importance(modFit.rf,  measures = NULL)
      vars <- randomForestExplainer::important_variables(importance_frame, k = 10, measures = c("mean_min_depth", "no_of_trees"))
      interactions_frame <- randomForestExplainer::min_depth_interactions(modFit.rf, vars)
      randomForestExplainer::plot_min_depth_interactions(interactions_frame)
    })

  })
}

## To be copied in the UI
# mod_descriptor_ui("descriptor_1")

## To be copied in the server
# mod_descriptor_server("descriptor_1")
