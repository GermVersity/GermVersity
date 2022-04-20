#' genotypic UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList

mod_genotypic_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(tags$h3("Germplasm bank management"),
                          tags$br(),
                          tags$p('In this demo, we are going to show how
                          to conduct basic genetic diversity
                          analyses of SNP data in a sample of
                          genebank accessions of Lima bean to
                          explore the genetic structure of the
                          sample. In this demo we are going to
                          learn how to calculate distance matrices
                          among individuals and populations,
                          visualize the distance matrices using
                          clustering algorithms (UPGMA and NJ),
                          how to conduct a principal component
                          analyses and a discriminant analysis
                          of principal components to assign
                          individuals into a k number of populations.'),
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
      mainPanel(
        tabsetPanel(
          type = 'tabs',
          shiny::tabPanel('Traditional Methods',
                          tags$h2("Results traditional methods"),
                          tags$br(),
                          tags$h4('Let´s calculate a Euclidean distance matrix
                                  between individuals on the basis of the
                                  observed allele frequencies within individuals.
                                  This distance is not a genetic
                                  distance but a geometric distance since it
                                  does not assume any evolutionary model.'),
                          tags$br(),
                          tags$h4("Circular dendrogram"),
                          shiny::plotOutput(ns("plot_data"), height = "500px",
                                            dblclick = "double_click",
                                            brush = brushOpts(
                                              id = "brush_plot",
                                              resetOnNew = TRUE)),
                          tags$h4("Horizontal dendrogram"),
                          shiny::plotOutput(ns("plot_data2"), height = "500px",
                                            dblclick = "double_click2",
                                            brush = brushOpts(
                                              id = "brush_plot",
                                              resetOnNew = TRUE))
          ),
          shiny::tabPanel('Distance Matrix',
                          tags$br(),
                          tags$p("Genetic distances measure how genetically
                                 similar are individuals or populations and
                                 are based on a specific evolutionary model.
                                 Genetic distances can be seen as summary
                                 statistics because they take the whole data
                                 set (for example, data from SNP loci) and
                                 summarize the genetic differentiation between
                                 samples (individuals or populations) in a
                                 value. One of the most used genetic distances
                                 is Nei's genetic distance which measures
                                 the genetic distance among populations on
                                 the basis of their genetic identity, namely
                                 the proportion of alleles that are shared
                                 between populations. For our data set, we
                                 have defined gene pools as populations in the
                                 genind object called LimaBeanData2.
                                 To calculate Nei's genetic distances among
                                 populations we need to add the strata to
                                 the genind object (in this case on the basis
                                 of the gene pool) to define the populations.
                                 Finally, we will build with this distance
                                 matrix a UPGMA tree and a neighbor-joining (
                                 NJ) topology, applying 1000 bootstrap permutations
                                 to get statistical support of the groups in the topologies."),
                          tags$br(),
                          tags$h4("UPGMA topology"),
                          shiny::plotOutput(ns("distance_matrix1"), height = "500px"),
                          tags$h4("NJ topology"),
                          shiny::plotOutput(ns("distance_matrix2"), height = "500px")
          ),
          shiny::tabPanel('PCA',
                          tags$br(),
                          tags$p('PCA (Principal component analysis) is a multivariate
                                  technique that summarizes
                                 the information provided by the genetic markers
                                 (for example, SNPs) into a few set of components.
                                 We can apply PCA to genlight objects with the
                                 function glPca of the package adgenet'),
                          tags$br(),
                          tags$h4("PCA plot"),
                          shiny::plotOutput(ns("plotPCA"), height = "500px")
          ),
          shiny::tabPanel('DPCA',
                          tags$br(),
                          tags$p('To identify genetic clusters, we can apply
                                 another multivariate approach known as DAPC.
                                 This approach is convenient when we are more
                                 interested in describing the diversity
                                 among groups of individuals than within groups.
                                 This approach is focused on finding discriminant
                                 functions that better describe the differences
                                 among groups, while minimizing the differences
                                 within groups. However, to find the discriminant
                                 functions, DAPC needs the groups to be known a
                                 priori and in many cases we just do not know
                                 how many groups are present in our sample.
                                 To address this issue, the adegenet package
                                 implements the function find.clusters (to find
                                 the number of clusters by running first a
                                 PCA and then using the k-means algorithm)
                                 and the function dapc to establish how are
                                 the relationships among the clusters.'),
                          tags$br(),
                          tags$h4("DPCA plot"),
                          shiny::plotOutput(ns("plotDPCA"), height = "500px"),
                          tags$h4("Clusters"),
                          shiny::dataTableOutput(ns("table1")),
                          tags$h4("DPCA custom plot"),
                          shiny::plotOutput(ns("plotDPCA2"), height = "500px"),
                          tags$h4("Membership probabilities"),
                          shiny::plotOutput(ns("plotDPCA3"), height = "500px")
          ),
          shiny::tabPanel('BDSpl',
                          tags$h3("Basic diverity statistics per locus"),
                          tags$br(),
                          tags$p('In this section we are going to estimate
                                 basic diversity statistics per locus in the
                                 whole sample. We are going to use the program
                                 adegenet to estimate the number of alleles
                                 per locus (NA), observed heterosygosity per
                                 locus (Hobs) and expected (Hexp) heterozygosity
                                 per locus. For SNP makers, NA might not be
                                 very useful because these loci are expected
                                 to be biallelic in populations (according to
                                 the infinite site mutational model) and we
                                 have also filtered the vcf file to include
                                 only biallelic SNPs. Hobs is the proportion
                                 of heterozygous individuals that were observed
                                 in the locus. Hexp is the heterozygosity we
                                 expect to observe in the locus assuming that
                                 the population (in this case the whole sample)
                                 is in Hardy-Weinberg equilibrium (HWE)in that
                                 locus. A significant difference among Hobs and
                                 Hexp means that the locus is not in HWE'),
                          tags$br(),
                          tags$h4("Hobs per locus"),
                          shiny::plotOutput(ns("plot_d1"), height = "500px"),
                          tags$h4("Hexp per locus"),
                          shiny::plotOutput(ns("plot_d2"), height = "500px"),
                          tags$h4("Hexp as a function of Hobs per locus"),
                          shiny::plotOutput(ns("plot_d3"), height = "500px"),
                          tags$h4("Testing the difference among Hexp and Hobs per locus"),
                          tags$p('We can apply Bartlett tests to assess whether Hexp is
                                 different from Hobs. As the number of loci, and
                                 therefore of multiple tests, is high, you may want
                                 to correct the p-value accordingly (for example with
                                 the Bonferroni correction).'),
                          shiny::verbatimTextOutput(ns("barlett"))

          ),
          shiny::tabPanel('BDSpc',
                          tags$h3("Basic diversity statistics per cluster (or population)"),
                          tags$br(),
                          tags$p('Once we have defined the number of clusters
                                 (or populations) in a sample, it is a good idea
                                 to estimate genetic diversity indexes for each
                                 cluster, for example if we want to find out
                                 which cluster is the most or least diverse.
                                 For doing this, we are going to use the program
                                 hierfstat to estimate observed heterozygosity
                                 (Ho) and mean gene diversities (Hs) within
                                 populations. A significant difference among
                                 Ho and Hs means that the population is not in
                                 HWE. Because different populations may not be
                                 in HWE for different reasons, diversity among
                                 populations may be hardly compared on the basis
                                 of Ho. This is because Hs is preferred to compare
                                 the genetic diversity between populations because
                                 with Hs all populations are assumed to be in HWE.
                                 It is also useful to estimate the fixation index
                                 Fis to measure the deviation to the assumption of
                                 HWE within each population.'),
                          tags$br(),
                          tags$h4("Summary Hs"),
                          tags$p('To summarize Hs values within populations. We can observe
                                 that populations are different in their genetic diversity,
                                 with the highest values for the domesticated admixed,
                                 followed by the Andean wild gene pool AI'),
                          tags$h4("Summary Ho"),
                          shiny::verbatimTextOutput(ns("summary1")),
                          tags$p('To summarize Ho values within populations. We can
                                 observe that within populations, Ho is very los, as
                                 expected for a selfing species as Lima bean.'),
                          tags$h4("Summary Fis"),
                          shiny::verbatimTextOutput(ns("summary2")),
                          tags$p('To summarize Fis values within populations.
                                 Fis is calculated as (Hs-Ho)/Hs. We can see
                                 that Fis values are positive for all populations
                                 and in some of them close to 1, as expected for
                                 a predominantly autogamous species as Lima bean.'),
                          shiny::verbatimTextOutput(ns("summary3")),
                          tags$p('A useful plot to compare genetic diversity among genepools'),
                          shiny::plotOutput(ns('boxplot'))
          ),
          shiny::tabPanel('Results ML',
                          tags$h3("Genetic divergence among populations"),
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
                          tags$p('A useful plot to explore the overall difference between
                                 Hs and Ht. We can see that Hs is lower than Ht, therefore
                                 we can envision that there will be certain degree of
                                 differentiation among populations (in this case we are
                                 comparing genepools in Lima bean) as measured by Fst
                                 (Fst=(Ht-Hs)/Ht).'),
                          shiny::plotOutput(ns('boxplot1')),
                          tags$p('to calculate Fis and Fst on one level hierarchy
                                 (the populations -or genepools- grouped within the
                                 total population) by the method of Weir and Cockerham (1984)'),
                          shiny::verbatimTextOutput(ns("summary4")),
                          tags$p('To obtain pariwise Fst values among genepools'),
                          shiny::verbatimTextOutput(ns("summary5")),
                          tags$p('Cluster Validity by NbCLust (and factoextra)
                                using 30 indices from the scientific literature'),
                          shiny::plotOutput(ns('plot1')),
                          tags$p('Cluster Validity by OptCluster (an improvement of ClValid)'),
                          shiny::plotOutput(ns('plot2')),
                          tags$p('Vizualization of OptCluster Output by means of UPGMA dendogram'),
                          shiny::plotOutput(ns('plot3'))
          )
        )
      )
    )
  )
}

#' genotypic Server Functions
#'
#' @import adegenet
#' @import hierfstat
#' @import ape
#' @import poppr
#' @import pegas
#' @import vcfR
#' @import RColorBrewer
#' @import ggplot2
#'
#' @noRd
mod_genotypic_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    options(shiny.maxRequestSize=30*1024^2)
    ranges <- reactiveValues(x = NULL, y = NULL)

    observeEvent(input$double_click, {
      brush <- input$brush_plot
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })

    output$plot_data = renderPlot({
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
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                   pop= Genepool,
                                   NA.char= "NA")

      EuclideanDistance <- dist(LimaBeanData2,
                                method = "euclidean",
                                diag = FALSE,
                                upper = FALSE,
                                p=2)
      NJtree <- ape::nj(EuclideanDistance)

      mycol = c("light blue",
                "gray",
                "green",
                "red",
                "blue",
                "light green",
                "pink")[LimaBeanData2$pop]

      plot(NJtree, tip.color=mycol, type="fan",
           x.lim = ranges$x,
           y.lim = ranges$y)
    })

    ranges1 <- reactiveValues(x = NULL, y = NULL)

    output$plot_data2 = renderPlot({
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
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                   pop= Genepool,
                                   NA.char= "NA")

      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS) # convert the vcf file into a genlight object
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      UPGMAtree <- poppr::aboot(LimaBeanData3, tree = "upgma", distance = poppr::bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

      mycol = c("light blue",
                "gray",
                "green",
                "red",
                "blue",
                "light green",
                "pink")[LimaBeanData2$pop]

      plot(UPGMAtree, tip.color=mycol, cex = 1.2,
           x.lim = ranges1$x,
           y.lim = ranges1$y)
    })

    observeEvent(input$double_click2, {
      brush <- input$brush_plot
      if (!is.null(brush)) {
        ranges1$x <- c(brush$xmin, brush$xmax)
        ranges1$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges1$x <- NULL
        ranges1$y <- NULL
      }
    })

    output$distance_matrix1 <- renderPlot({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop(safeError(e))
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop(safeError(e))
        }
      )
      Genepool <- as.character(data$Genepool)
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char= "NA")

      adegenet::strata(LimaBeanData2) <- data.frame(adegenet::pop(LimaBeanData2)) #to add the strata of a genind object
      adegenet::nameStrata(LimaBeanData2) <- ~Genepool # to assign a name to the strata of a genind object
      #Building of a UPGMA topology with bootstrap support using the function aboot of the poppr package
      set.seed(999)
      poppr::aboot(LimaBeanData2, strata = Genepool, tree= "upgma", distance="nei.dist", sample=1000, cutoff=0, showtree=TRUE,quiet=FALSE)

    })

    output$distance_matrix2 <- renderPlot({
      tryCatch(
        {
          LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop(safeError(e))
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop(safeError(e))
        }
      )
      Genepool <- as.character(data$Genepool)
      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")

      adegenet::strata(LimaBeanData2) <- data.frame(adegenet::pop(LimaBeanData2)) #to add the strata of a genind object
      adegenet::nameStrata(LimaBeanData2) <- ~Genepool # to assign a name to the strata of a genind object
      #Building of a UPGMA topology with bootstrap support using the function aboot of the poppr package
      set.seed(999)
      #Building of a NJ topology with bootstrap support using the function aboot of the poppr package
      poppr::aboot(LimaBeanData2, strata = Genepool, tree= "nj", distance="nei.dist", sample=1000, cutoff=0, showtree=TRUE,quiet=FALSE)
    })

    output$plotPCA <- renderPlot({
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

      Genepool = as.character(data$Genepool)
      LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
      adegenet::ploidy(LimaBeanData3) <- 2
      adegenet::pop(LimaBeanData3) <- Genepool

      LimaBeanPCA <- adegenet::glPca(LimaBeanData3, nf = 3)
      # to carry out a PCA on a genlight object. With the argument nf as NULL,
      #you are asked interactively for the number of principal components to be
      #retained. For this data, three axes were retained.

      # to create a customized PCA plot with the package ggplot2

      Limapcascores <- as.data.frame(LimaBeanPCA$scores)
      Limapcascores$pop <- adegenet::pop(LimaBeanData3)


      set.seed(5)
      colors <- RColorBrewer::brewer.pal(n = adegenet::nPop(LimaBeanData3), name = "Set1")
      p <- ggplot2::ggplot(Limapcascores, ggplot2::aes(x=PC1, y=PC2, colour=pop))
      p <- p + ggplot2::geom_hline(yintercept = 0)
      p <- p + ggplot2::geom_vline(xintercept = 0)
      p <- p + ggplot2::geom_point(size=3)
      p <- p + ggplot2::theme_bw()
      p <- p + ggplot2::scale_color_manual(values=colors ) +
        ggplot2::xlab(sprintf("PC1 %f percent", 100*LimaBeanPCA$eig[1]/sum(LimaBeanPCA$eig))) +
        ggplot2::ylab(sprintf("PC2 %f percent", 100*LimaBeanPCA$eig[2]/sum(LimaBeanPCA$eig)))
      p

    })

    output$plot_d1 <- renderPlot({

      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file') #safeError(e)
        }
      )

      tryCatch(
        {
          PopData = read.table(input$filetxt$datapath, header = T)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )
      Genepool <- as.character(PopData$Genepool)

      require("adegenet")

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")

      diversity <- summary(LimaBeanData2)#we use the genind object created above
      print(class(diversity))
      plot(diversity['Hobs'][[1]],
           xlab="loci number",
           ylab="Hobs")
      #As expected for a selfer as Lima bean, in most loci the number of
      #heterosygous individuals is low

    })

    output$plot_d2 <- renderPlot({

      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file') #safeError(e)
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )
      Genepool <- as.character(data$Genepool)

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char = "NA")
      diversity <- base::summary(LimaBeanData2)#we use the genind object created above
      plot(diversity['Hexp'][[1]],
           xlab="loci number",
           ylab="Hexp",
           ylim = c(0,0.6))
      #A great variation in Hexp is observed. The maximum value for Hexp is 0.5
      #and it is observed when in a locus both alleles have the same frequency of 0.5.

    })

    output$plot_d3 <- renderPlot({
      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file') #safeError(e)
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )
      Genepool <- as.character(data$Genepool)

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char = "NA")
      diversity <- base::summary(LimaBeanData2)#we use the genind object created above

      plot(diversity$Hobs, diversity$Hexp, xlab="Hobs", ylab="Hexp")
      #This plot shows how Hexp is related to Hobs. For loci in HWE,
      #Hexp = Hobs. In the plot we observe that in most loci, Hexp is higher
      #than Hobs, as expected for selfing species as Lima bean.

    })

    output$barlett <- renderText({
      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop('Upload VCF file') #safeError(e)
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath)
        },
        error = function(e){
          stop('Upload VCF file')
        }
      )
      Genepool <- as.character(data$Genepool)

      LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                         pop = Genepool,
                                         NA.char = "NA")
      diversity <- base::summary(LimaBeanData2)#we use the genind object created above

      bartlett.test(list(diversity$Hobs, diversity$Hexp))
    })

    data1 <- reactive({
      tryCatch(
        {
          LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
        },
        error = function(e){
          stop("Upload VCF file")
        }
      )

      tryCatch(
        {
          data = read.table(input$filetxt$datapath, header = TRUE)
        },
        error = function(e){
          stop("Upload TXT file")
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

    output$summary1 <- renderPrint({

      LimaBeanData2 = data1()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      print(summary(diversity.clusters$Hs))
    })

    output$summary2 <- renderPrint({
      LimaBeanData2 = data1()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      print(summary(diversity.clusters$Ho))
    })

    output$summary3 <- renderPrint({
      LimaBeanData2 = data1()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      print(summary(diversity.clusters$Fis))
    })

    output$boxplot <- renderPlot({
      LimaBeanData2 = data1()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      boxplot(diversity.clusters$Hs, ylab="Hs")
    })

    data3 <- reactive({
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

    output$boxplot1 <- renderPlot({

      LimaBeanData2 <- data3()

      diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
      boxplot(diversity.clusters$perloc[,2:3])
    })

    output$summary4 <- renderPrint({

      LimaBeanData2 <- data3()

      global.Fst.weir_cock <- hierfstat::wc(LimaBeanData2)
      #to calculate Fis and Fst on one level hierarchy (the populations -
      #or genepools- grouped within the total population) by the method of Weir and Cockerham (1984).

      # We can see that genetic differentiation among
      #genepools is high (Fst=0.59), and that deviation from HWE is also high,
      #with deficit of heterozygous (Fis=0.69).

      print(global.Fst.weir_cock)
    })

    output$summary5 <- renderPrint({

      LimaBeanData2 <- data3()

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
# mod_genotypic_ui("genotypic_1")

## To be copied in the server
# mod_genotypic_server("genotypic_1")
