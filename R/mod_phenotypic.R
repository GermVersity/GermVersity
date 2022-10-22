#' phenotypic UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phenotypic_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(tags$h3("Phenotypic analysis in genebanks"),
                   tags$p('This work is in the framework of the project "
                          To use AI tools for multidimensional data analysis
                          to understand and promote biodiversity in plant
                          germplasm banks (Phaseolus case)" led by AGROSAVIA,
                          UNAL and Alianza Bioversity-CIAT from a cooperation
                          agreement. the project objective is to "To use AI
                          tools for multidimensional data analysis to understand
                          and promote biodiversity in plant germplasm banks
                          (Phaseolus case).'),
                   shiny::fileInput(ns("filecsv"), "Choose a CVS file",
                                    multiple = F,
                                    accept = c("text/csv",
                                               "text/comma-separated-values, text/plain",
                                               ".csv"),
                                    buttonLabel = "Uploading...")
      ),
      mainPanel(
        tags$h3("Model"),
        verbatimTextOutput(ns('modelF')),
        plotOutput(ns("modFit")),
        tags$h3('Confusion matrix'),
        tags$p('Table of confusion of the training data of
                       the common bean races classification.'),
        plotOutput(ns('con_matrix'))
      )
    )
  )
}

#' phenotypic Server Functions
#'
#' @import readr
#' @import dplyr
#' @import randomForest
#' @import randomForestExplainer
#' @import ggplot2
#' @import DT
#' @import caret
#' @import magrittr
#'
#' @noRd
mod_phenotypic_server <- function(id){
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

      data %>% dplyr::group_by(Races) %>% dplyr::summarise(n = dplyr::n())

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
    })

    output$modelF <- renderPrint({
      modFitN <- model()

      print(modFitN)

    })

    output$modFit <- renderPlot({

      modFitN <- model()

      plot(modFitN)

    })

    output$con_matrix <- renderPlot({

      modFitN <- model()

      modFit.rfN <- randomForest::randomForest(Races ~., data = training[],
                                               mtry= modFitN$bestTune$mtry)

      ### Confusion matrix ###

      TablaE <- as.table(modFit.rfN$confusion)
      TablaE <- as.data.frame(TablaE)

      TablaE <- TablaE[1:36, ]

      ### Training confusion ###

      ggplot2::ggplot(TablaE, aes(x=Var1, y=Var2, fill=Freq)) +
        ggplot2::geom_tile(color="black") +
        ggplot2::theme_bw() +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_distiller(palette="Greys", direction=1) +
        ggplot2::guides(fill=F) + # removing legend for `fill`
        ggplot2::labs(title = "Training Confusion Matrix") + # using a title instead
        ggplot2::geom_text(aes(label=Freq), color="black") +
        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank())

    })

  })
}

## To be copied in the UI
# mod_phenotypic_ui("phenotypic_1")

## To be copied in the server
# mod_phenotypic_server("phenotypic_1")
