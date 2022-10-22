#' home_page UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_home_page_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        tags$h3(class = "title",
                "Application in development..."),
        tags$h3(),
        tags$p(class = "texto_inicial",
               "This project transforms the way we address environmental
               challenges by more systematically and consistently studying
               the biodiversity of crop gene pools. Platforms for the
               discovery of agrobiodiversity, such as the one proposed here
               using AI, will allow the adaptability of crops not to depend
               exclusively on their intrinsic diversity, probably without
               sufficient previous adaptations to unusual climates.This project
               will harness modern AI technology to assist in the identification,
               characterization and utilization of new adaptive sources among
               thousands of landraces and related wild species, all stored in genebanks.."),
        tags$h4(class = "version",
                "Version 1.8")
      ),
      mainPanel(
        tags$div(class = "images",
                 tags$img(src = "https://alimentro.agrosavia.co/Content/imagenes/logo-agrosavia.png", width = "230px", height = "50px"),
                 tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0a/Logotipo_de_la_Universidad_Nacional_de_Colombia.svg/1200px-Logotipo_de_la_Universidad_Nacional_de_Colombia.svg.png", width = "150px", height = "150px"),
                 tags$img(src = "https://upload.wikimedia.org/wikipedia/en/thumb/f/f4/Kansas_Jayhawks_logo.svg/1200px-Kansas_Jayhawks_logo.svg.png", width = "150px", height = "150px"),
                 tags$img(src = "https://climaloca.org/wp-content/uploads/alliance_logo_standard_cropped.png", width = "150px",  height = "70px")),
        tags$h1(class = "title",
                "Welcome to GermVersity"),
        tags$h1(),
        tags$h4("Artificial intelligence as a tool for
                diversity analysis of plant genetic
                resources in genebanks"),
        tags$p(class = "texto_inicial",
               "Plant genetic resources conserved in genebanks have
               the potential to harbor agronomically important traits
               that can be used to overcome current and future
               environmental challenges. However, despite their
               importance, these resources remain largely underutilized
               around the world. During the last years, the development
               of genome sequencing technologies has opened new avenues
               in the discovery of relevant genes for various agronomic
               traits in multiple crop species. Genebanks and the
               scientific community have generated massive amounts of
               phenotypic, genotypic, climatic, and geospatial data
               related to germplasm accessions. Currently, the challenge
               is to efficiently use these high-dimensional datasets to
               discover the true potential of plant accessions for breeding
               purposes. Although artificial intelligence (AI) tools may
               outperform competing methods in extracting hidden insights
               from high-dimensional data, their use in studying genebank
               accessions is still in its infancy. This situation could
               be due to the fact that these tools are relatively new
                and applications in the field of plant genetic resources
               are still incipient. To address this issue, we propose to
               use artificial intelligence tools to comprehensively
               analyze the phenotypic, genomic, and environmental data
               associated with plant germplasm accessions held in genebanks
               to gain a deep understanding of their genetic diversity,
               discover useful germplasm in climate change scenarios and
               promote their use.")
      )
    )
  )
}

#' home_page Server Functions
#'
#' @noRd
mod_home_page_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_home_page_ui("home_page_1")

## To be copied in the server
# mod_home_page_server("home_page_1")
