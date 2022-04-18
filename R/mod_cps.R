#' cps UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_cps_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Correcting for population stratification"),
                          tags$br(),
                          tags$p('We are going to use the approach of Zhao et al.
                                 (2012) which consists in first, regressing
                                 (at each SNP locus and by applying a generalized
                                 linear model, GLM), the environmental variable
                                 and the genotype against population information.
                                 The residuals of the models are then used to
                                 adjust the values of the environmental variable
                                 and the genotype before performing the Random
                                 Forest analysis.'),
                          tags$br(),

                          shiny::fileInput(ns("filelfmm"), "Choose a LFMM file",
                                           multiple = F,
                                           accept = ".lfmm",
                                           buttonLabel = "Uploading..."),
                          tags$br(),

                          shiny::fileInput(ns("filetxt1"), "Choose a TXT file (genepool)",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          tags$br(),

                          shiny::fileInput(ns("filetxt2"), "Choose a TXT file (environment)",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),
                          tags$br(),

                          shiny::fileInput(ns("filetxt3"), "Choose a TXT file (position)",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading..."),

                          tags$br()
      ),
      mainPanel(
        shiny::tabsetPanel(
          type = 'tabs',
          shiny::tabPanel('Distribution',
                          tags$br(),
                          tags$h3('Explore the overall distribution of the phenotype'),
                          plotOutput(ns('plot_data')),
                          tags$br(),
                          tags$h3('Should all be zero (0)'),
                          verbatimTextOutput(ns('corrected'))
                          ),
          shiny::tabPanel('Random Forest',
                          tags$br(),
                          tags$p('RF algorithms should be optimized in order to
                                 maximize the percentage of variance explained
                                 (PVE, for continuos variables) and to ensure results
                                 are repeatable. We are going to adjuts two parameters:
                                 ntree (number of trees grown per forest) and mtry
                                 (number of predictors to randomly sample at each node).
                                 We should establish the values of ntree and mtry where
                                 PVE reaches a plateau. We will try values of ntree from
                                 100 to 10000 (in increments of 100), and values of mtry
                                 of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p (
                                 p=number of loci). '),
                          tags$br(),
                          plotOutput(ns('plot_data2'))
                          ),
          shiny::tabPanel('Correlation',
                          tags$br(),
                          plotOutput(ns('hist1'))),
          shiny::tabPanel('Results SNP',
                          tags$br(),
                          DT::dataTableOutput(ns('data_table')),
                          tags$br(),
                          plotOutput(ns('manhattan')))
        )

      )
    )
  )
}

#' cps Server Functions
#'
#' @import randomForest
#' @import Pomona
#' @import ranger
#' @import parallel
#'
#' @noRd
mod_cps_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$plot_data <- renderPlot({
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
          bio10 <- read.table(input$filetxt2$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          PopData <- read.table(input$filetxt1$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      pop <- PopData$pop

      #Create dataframe with genetic data, environmental data and SNP matrix
      data<-as.data.frame(cbind(pop,bio10,Y1))

      colnames(data)[1:2]<-c("pop","bio10")

      # First, explore the overall distribution of the phenotype
      hist(data$bio10, xlab = 'Phenotype', main = '')

    })

    output$corrected <- renderPrint({

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
          bio10 <- read.table(input$filetxt2$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          PopData <- read.table(input$filetxt1$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      pop <- PopData$pop

      #Create dataframe with genetic data, environmental data and SNP matrix
      data<-as.data.frame(cbind(pop,bio10,Y1))

      data_corrected <- data

      data_corrected[,2:ncol(data_corrected)] <- NA

      for (i in 2:ncol(data)){
        LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
        data_corrected[,i] <- LM_SNP_i$residuals
        colnames(data_corrected)[i]<-colnames(data)[i]
        #if(i%%50==0) print(i)
      }

      data_corrected[,ncol(data_corrected)]-LM_SNP_i$residuals
    })


    output$plot_data2 <- renderPlot({

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
          bio10 <- read.table(input$filetxt2$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          PopData <- read.table(input$filetxt1$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      pop <- PopData$pop

      #Create dataframe with genetic data, environmental data and SNP matrix
      data <- as.data.frame(cbind(pop,bio10,Y1))

      colnames(data)[1:2] <- c("pop","bio10")

      data_corrected <- data

      data_corrected[,2:ncol(data_corrected)] <- NA

      for (i in 2:ncol(data)){
        LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
        data_corrected[,i] <- LM_SNP_i$residuals
        colnames(data_corrected)[i]<-colnames(data)[i]
        #if(i%%50==0) print(i)
      }

      p <- ncol(data_corrected)-2

      results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)
      for (i in seq(from = 100, to = 1000 , by = 100)){
        print(i)
        for (j in c(sqrt(p), 2*sqrt(p), 0.1*p, 0.2*p, p/3, p)){
          rf_ij <- ranger::ranger(x = data_corrected[,3:ncol(data_corrected)],
                                  y = data_corrected$bio10,
                                  num.trees=i,
                                  mtry=j)
          results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$r.squared,1)))
        }
      }

      # Clean up the file format
      results_optimization<-as.data.frame(results_optimization)
      colnames(results_optimization)<-c("ntree", "mtry","PVE")

      #write.csv(results_optimization,file="lima_bean_results.optimization.csv", row.names = FALSE, quote=FALSE)

      # Now plot results to see if there's a plateau


      plot(results_optimization$ntree[results_optimization$mtry == sqrt(p)],results_optimization$PVE[results_optimization$mtry == sqrt(p)], type="l", col="black", xlab="ntree",ylab="PVE", ylim=range(results_optimization$PVE))
      lines(results_optimization$ntree[results_optimization$mtry == 2*sqrt(p)],results_optimization$PVE[results_optimization$mtry == 2*sqrt(p)], col="blue")
      lines(results_optimization$ntree[results_optimization$mtry == 0.1*p],results_optimization$PVE[results_optimization$mtry == 0.1*p], col="green")
      lines(results_optimization$ntree[results_optimization$mtry == 0.2*p],results_optimization$PVE[results_optimization$mtry == 0.2*p], col="purple")
      lines(results_optimization$ntree[results_optimization$mtry == p/3],results_optimization$PVE[results_optimization$mtry == p/3], col="orange")
      lines(results_optimization$ntree[results_optimization$mtry == p],results_optimization$PVE[results_optimization$mtry == p], col="red")


    })

    output$hist1 <- renderPlot({
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
          bio10 <- read.table(input$filetxt2$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          PopData <- read.table(input$filetxt1$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          pos <- read.table(input$filetxt3$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      pop <- PopData$pop

      #Create dataframe with genetic data, environmental data and SNP matrix
      data <- as.data.frame(cbind(pop,bio10,Y1))

      colnames(data)[1:2] <- c("pop","bio10")

      data_corrected <- data

      data_corrected[,2:ncol(data_corrected)] <- NA

      for (i in 2:ncol(data)){
        LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
        data_corrected[,i] <- LM_SNP_i$residuals
        colnames(data_corrected)[i]<-colnames(data)[i]
        #if(i%%50==0) print(i)
      }

      #Run 10 random forests with 1000 trees each and the optimized mtry value (0.1)
      r2vim <- Pomona::var.sel.r2vim(x = data_corrected[,3:ncol(data_corrected)],
                                     y=data_corrected$bio10,
                                     no.runs = 10,
                                     factor = 1,
                                     ntree = 1000,
                                     mtry.prop = 0.1,
                                     method = "ranger",
                                     type = "regression")

      #calculate correlation between runs
      runs <- as.data.frame(cbind(r2vim$info$vim.run.1,
                                  r2vim$info$vim.run.2,
                                  r2vim$info$vim.run.3,
                                  r2vim$info$vim.run.4,
                                  r2vim$info$vim.run.5,
                                  r2vim$info$vim.run.6,
                                  r2vim$info$vim.run.7,
                                  r2vim$info$vim.run.8,
                                  r2vim$info$vim.run.9,
                                  r2vim$info$vim.run.10))

      hist(cor(runs), col="blue", xlab = 'Correlation', main= '')

    })

    output$data_table <- DT::renderDataTable({
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
          bio10 <- read.table(input$filetxt2$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          PopData <- read.table(input$filetxt1$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          pos <- read.table(input$filetxt3$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      pop <- PopData$pop

      #Create dataframe with genetic data, environmental data and SNP matrix
      data <- as.data.frame(cbind(pop,bio10,Y1))

      colnames(data)[1:2] <- c("pop","bio10")

      data_corrected <- data

      data_corrected[,2:ncol(data_corrected)] <- NA

      for (i in 2:ncol(data)){
        LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
        data_corrected[,i] <- LM_SNP_i$residuals
        colnames(data_corrected)[i]<-colnames(data)[i]
        #if(i%%50==0) print(i)
      }

      #Run 10 random forests with 1000 trees each and the optimized mtry value (0.1)
      r2vim <- Pomona::var.sel.r2vim(x = data_corrected[,3:ncol(data_corrected)],
                                     y=data_corrected$bio10,
                                     no.runs = 10,
                                     factor = 1,
                                     ntree = 1000,
                                     mtry.prop = 0.1,
                                     method = "ranger",
                                     type = "regression")
      # Build a data frame with results
      df<-data.frame(chrom = pos$chrom,
                     bp = pos$bp,
                     min.rel.vim = r2vim$info$rel.vim.min,
                     snp_name = pos$snp)
      DT::datatable(df,
                    filter = 'top',
                    extensions = c('Buttons', 'Scroller'),
                    options = list(scrollY = 650,
                                   scrollX = 500,
                                   deferRender = TRUE,
                                   scroller = TRUE,
                                   # paging = TRUE,
                                   # pageLength = 25,
                                   buttons = list('excel',
                                                  list(extend = 'colvis',
                                                       targets = 0,
                                                       visible = FALSE)),
                                   dom = 'lBfrtip',
                                   fixedColumns = TRUE),
                    rownames = FALSE)
    })

    output$manhattan <- renderPlot({
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
          bio10 <- read.table(input$filetxt2$datapath, quote="\"", comment.char="")
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          PopData <- read.table(input$filetxt1$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      tryCatch(
        {
          pos <- read.table(input$filetxt3$datapath, header = T)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )

      pop <- PopData$pop

      #Create dataframe with genetic data, environmental data and SNP matrix
      data <- as.data.frame(cbind(pop,bio10,Y1))

      colnames(data)[1:2] <- c("pop","bio10")

      data_corrected <- data

      data_corrected[,2:ncol(data_corrected)] <- NA

      for (i in 2:ncol(data)){
        LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
        data_corrected[,i] <- LM_SNP_i$residuals
        colnames(data_corrected)[i]<-colnames(data)[i]
        #if(i%%50==0) print(i)
      }

      #Run 10 random forests with 1000 trees each and the optimized mtry value (0.1)
      r2vim <- Pomona::var.sel.r2vim(x = data_corrected[,3:ncol(data_corrected)],
                                     y=data_corrected$bio10,
                                     no.runs = 10,
                                     factor = 1,
                                     ntree = 1000,
                                     mtry.prop = 0.1,
                                     method = "ranger",
                                     type = "regression")
      # Build a data frame with results
      df<-data.frame(chrom = pos$chrom,
                     bp = pos$bp,
                     min.rel.vim = r2vim$info$rel.vim.min,
                     snp_name = pos$snp)

      #plot Manhattan plot with variable importance for each SNP
      plot(x=1:nrow(df),y=r2vim$info$rel.vim.min, xlab = 'SNP', ylab = 'Importance')
      abline(h=1, col="red")
    })

  })
}

## To be copied in the UI
# mod_cps_ui("cps_1")

## To be copied in the server
# mod_cps_server("cps_1")
