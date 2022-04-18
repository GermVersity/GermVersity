#' team UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_team_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarLayout(position = "right",
                  sidebarPanel(
                    tags$h5(class="title",
                            "Versión 2.0.1"),
                    tags$p(class = "footer",
                           "© Copyright GermVersity. All rights reserved")
                  ),
                  mainPanel(
                    tags$div(
                      class = "container1",
                      tags$h1(class = "title_i",
                              "Members GermVersity"),
                      tags$div(class = "container",
                               tags$div(class = "Juan",
                                        tags$h3("Juan Camilo Henao-Rojas"),
                                        tags$img(src = "https://alimentro.agrosavia.co/Content/imagenes/logo-agrosavia.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Joaquin",
                                        tags$h3("Joaquin Guillermo-Ramirez"),
                                        tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0a/Logotipo_de_la_Universidad_Nacional_de_Colombia.svg/1200px-Logotipo_de_la_Universidad_Nacional_de_Colombia.svg.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Andres",
                                        tags$h3("Andres Cortes-Vera"),
                                        tags$img(src = "https://alimentro.agrosavia.co/Content/imagenes/logo-agrosavia.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Maria",
                                        tags$h3("Maria Isabel Chacon-Sanchez"),
                                        tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0a/Logotipo_de_la_Universidad_Nacional_de_Colombia.svg/1200px-Logotipo_de_la_Universidad_Nacional_de_Colombia.svg.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Diego",
                                        tags$h3("Diego Felipe Conejo"),
                                        tags$img(src = "https://climaloca.org/wp-content/uploads/alliance_logo_standard_cropped.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Marlon",
                                        tags$h3("Marlon E. Cobos"),
                                        tags$img(src = "https://upload.wikimedia.org/wikipedia/en/thumb/f/f4/Kansas_Jayhawks_logo.svg/1200px-Kansas_Jayhawks_logo.svg.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Luis",
                                        tags$h3("Luis Felipe Lopez"),
                                        tags$img(src = "https://alimentro.agrosavia.co/Content/imagenes/logo-agrosavia.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Paula",
                                        tags$h3("Paula Helena Reyes Herrera"),
                                        tags$img(src = "https://alimentro.agrosavia.co/Content/imagenes/logo-agrosavia.png", width = "200px")),
                               tags$h1(),
                               tags$div(class = "Kevin",
                                        tags$h3("Kevin Steven Quiroga Benavides"),
                                        tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0a/Logotipo_de_la_Universidad_Nacional_de_Colombia.svg/1200px-Logotipo_de_la_Universidad_Nacional_de_Colombia.svg.png", width = "200px"))

                      )
                    )

                  )
    )
  )
}

#' team Server Functions
#'
#' @noRd
mod_team_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_team_ui("team_1")

## To be copied in the server
# mod_team_server("team_1")
