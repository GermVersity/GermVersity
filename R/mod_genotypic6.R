#' genotypic6 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_genotypic6_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Genetic divergence among populations"),
                          tags$br(),
                          tags$p('We may also be interested in finding out how genetically
                                 divergent are the populations. For this, we can calculate
                                 the fixation index Fst, which measures the difference
                                 between the mean gene diversity within populations (Hs)
                                 and the total genetic diversity expected in the whole
                                 population (Ht). The higher the difference between Ht
                                 and Hs, the higher the fixation index Fst and therefore
                                 the higher the populations we are comparing. We will
                                 use the program hierfstat to calculate Fst indexes.'),
                          tags$br(),
                          shiny::fileInput(ns("filevcf"), "Choose a VCF file",
                                           multiple = F,
                                           accept = ".vcf",
                                           buttonLabel = "Uploading..."),
                          shiny::fileInput(ns("filetxt"), "Choose a TXT file",
                                           multiple = F,
                                           accept = ".txt",
                                           buttonLabel = "Uploading...")
      ),
      mainPanel(tags$h2("Results ML"),
                tags$p('A useful plot to explore the overall difference between
                       Hs and Ht. We can see that Hs is lower than Ht, therefore
                       we can envision that there will be certain degree of
                       differentiation among populations (in this case we are
                       comparing genepools in Lima bean) as measured by Fst
                       (Fst=(Ht-Hs)/Ht).'),
                shiny::plotOutput(ns('boxplot')),
                tags$p('to calculate Fis and Fst on one level hierarchy
                       (the populations -or genepools- grouped within the
                       total population) by the method of Weir and Cockerham (1984)'),
                shiny::verbatimTextOutput(ns("summary1")),
                tags$p('To obtain pariwise Fst values among genepools'),
                shiny::verbatimTextOutput(ns("summary2")),
                tags$p('Cluster Validity by NbCLust (and factoextra)
                       using 30 indices from the scientific literature'),
                shiny::plotOutput(ns('plot1')),
                tags$p('Cluster Validity by OptCluster (an improvement of ClValid)'),
                shiny::plotOutput(ns('plot2')),
                tags$p('Vizualization of OptCluster Output by means of UPGMA dendogram'),
                shiny::plotOutput(ns('plot3')),
      )
    )
  )
}

#' genotypic6 Server Functions
#'
#' @import adegenet
#' @import hierfstat
#' @import ape
#' @import poppr
#' @import NbClust
#' @import pegas
#' @import vcfR
#' @import RColorBrewer
#' @import ggplot2
#' @import uwIntroStat
#' @import factoextra
#' @import optCluster
#'
#' @noRd
mod_genotypic6_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    data <- reactive({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )
      Genepool <- as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                  pop= Genepool,
                                  NA.char= "NA")
    })

    data2 <- reactive({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )
      Genepool <- as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool
      LimaBeanPCA <- adegenet::glPca(LimaBeanData3, nf = 3)
      Limapcascores <- as.data.frame(LimaBeanPCA$scores)
      Limapcascores$pop <- adegenet::pop(LimaBeanData3)
    })

    output$boxplot <- renderPlot({

      LimaBeanData2 <- data()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      boxplot(diversity.clusters$perloc[,2:3])
    })

    output$summary1 <- renderPrint({

      LimaBeanData2 <- data()

      global.Fst.weir_cock <- hierfstat::wc(LimaBeanData2)
      #to calculate Fis and Fst on one level hierarchy (the populations -
      #or genepools- grouped within the total population) by the method of Weir and Cockerham (1984).

      # We can see that genetic differentiation among
      #genepools is high (Fst=0.59), and that deviation from HWE is also high,
      #with deficit of heterozygous (Fis=0.69).

      print(global.Fst.weir_cock)
    })

    output$summary2 <- renderPrint({

      LimaBeanData2 <- data()

      pairwise.fst.genepools <- hierfstat::genet.dist(LimaBeanData2, method="WC84")
      # To obtain pariwise Fst values among genepools
      #We can see that the highest Fst values (highest genetic
      #differentiation) were observed among genepools Dom_MI vs Dom_AI,
      #Dom_MI vs Wild_AII and Dom_AI vs wild_AII. The lowest genetic
      #differentiation was observed among genepools Dom_ADM vs wild MII,
      #Wild_AI vs Dom_AI (Andean wild ancestor and their domesticated counterpart)
      #and Wild_MI and Dom_MI (Mesoamerican wild ancestor and their domesticated counterpart).

      print(pairwise.fst.genepools)
    })

    output$plot1 <- renderPlot({

      Limapcascores <- data2()

      distance = "euclidean" #"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
      algorithm = "kmeans" #"ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans"
      index = c("all")
      res.nbclust <- NbClust::NbClust(as.data.frame(Limapcascores[,-4]), distance = distance,
                             min.nc = 4, max.nc = 10,
                             method = algorithm, index =index)
      factoextra::fviz_nbclust(res.nbclust) + theme_minimal()
      res.nbclust$Best.nc[1,]
    })

    output$plot2 <- renderPlot({

      Limapcascores <- data2()

      algorithms =  c("pam","kmeans","clara","diana","agnes") #agglomerative =agnes &hclust "kmeans","clara","diana","agnes","som","sota"
      # Select the k range
      k = 4:10
      # Select the valitadion method: "CE", "GA"
      validation_method="CE"
      # run
      PVCA <-as.data.frame(Limapcascores[,1:3])
      set.seed(2022)
      norm1 <- optCluster::optCluster(as.matrix.data.frame(Limapcascores[,-4]), k, clMethods = algorithms,
                          rankMethod = validation_method, hierMethod = "ward")
      summary(norm1)

      putput <- optCluster::optAssign(norm1)
      CLUSTER <- as.data.frame(putput$cluster)
      colnames(CLUSTER) <- "CLUSTER"

      final_cluster_data <- cbind(PVCA$PC1, PVCA$PC2, CLUSTER)
      rownames(final_cluster_data) <- rownames(Limapcascores)
      final_cluster_data <- as.data.frame(final_cluster_data)
      final_cluster_data$CLUSTER <- as.factor(final_cluster_data$CLUSTER)
      colnames(final_cluster_data) <- c("V1","V2","CLUSTER")

      my_pal <- c("darkgreen","darkblue", "orangered","darkred","lightslateblue","orange","purple4","darkred","green","red","pink","yellow","black","deeppink4","darkturquoise","khaki3")
      my_fill <- c("darkgreen","darkblue", "orangered","darkred","lightslateblue","orange","purple4","darkred","green","red","pink","yellow","black","deeppink4","darkturquoise","khaki3")
      p1 <- ggplot2::ggplot(final_cluster_data, aes(x = V1, y = V2,color = CLUSTER))
      p1 <- p1 + ggplot2::geom_point(size = 3, aes(fill = CLUSTER),alpha =0.5)
      p1 <- p1 + ggplot2::geom_hline(yintercept = 0)
      p1 <- p1 + ggplot2::geom_vline(xintercept = 0)
      p1 <- p1 + ggplot2::geom_point(size=3)
      p1 <- p1 + ggplot2::theme_bw()
      p1 <- p1 + ggplot2::scale_color_manual(values=c(my_pal))
      p1 <- p1 + ggplot2::scale_fill_manual(values=c(paste(my_fill)))+
        ggplot2::xlab(sprintf("PC1 %f percent", 100*LimaBeanPCA$eig[1]/sum(LimaBeanPCA$eig))) +
        ggplot2::ylab(sprintf("PC2 %f percent", 100*LimaBeanPCA$eig[2]/sum(LimaBeanPCA$eig)))
      p1
    })

    output$plot3 <- renderPlot({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop('Upload TXT file')
        }
      )
      Genepool <- as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool
      Limapcascores <- data2()
      algorithms =  c("pam","kmeans","clara","diana","agnes") #agglomerative =agnes &hclust "kmeans","clara","diana","agnes","som","sota"
      # Select the k range
      k = 4:10
      # Select the valitadion method: "CE", "GA"
      validation_method="CE"
      # run
      PVCA <-as.data.frame(Limapcascores[,1:3])
      set.seed(2022)
      norm1 <- optCluster::optCluster(as.matrix.data.frame(Limapcascores[,-4]), k, clMethods = algorithms,
                                      rankMethod = validation_method, hierMethod = "ward")
      summary(norm1)

      putput <- optCluster::optAssign(norm1)
      CLUSTER <- as.data.frame(putput$cluster)
      colnames(CLUSTER) <- "CLUSTER"

      final_cluster_data <- cbind(PVCA$PC1, PVCA$PC2, CLUSTER)
      rownames(final_cluster_data) <- rownames(Limapcascores)
      final_cluster_data <- as.data.frame(final_cluster_data)
      final_cluster_data$CLUSTER <- as.factor(final_cluster_data$CLUSTER)
      colnames(final_cluster_data) <- c("V1","V2","CLUSTER")
      # metodo distancia
      adegenet::pop(LimaBeanData3) <- final_cluster_data$CLUSTER
      tree <- poppr::aboot(LimaBeanData3, tree = "upgma", distance = nei.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
      # cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")
      ape::plot.phylo(tree, cex = 0.3, font = 2, adj = 0, tip.color =  my_pal[adegenet::pop(LimaBeanData3)])
      nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.3,font = 3, xpd = TRUE)
      axis(side = 1)
      title(xlab = "Genetic distance (proportion of loci that are different)")
    })

  })
}

## To be copied in the UI
# mod_genotypic6_ui("genotypic6_1")

## To be copied in the server
# mod_genotypic6_server("genotypic6_1")
