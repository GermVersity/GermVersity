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
    fluidPage(
      tags$h1(class = "title", 'Miembros del equipo'),
      tags$br(),
      tags$div(class = "cardsteam",
               tags$div(class = "card",
                        tags$img(src = "www/Joaquin.jpg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Joaquin Guillermo Ramirez Gil')),
                                 tags$img(src = 'www/unalblack.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0002-0162-3598")
                        )
               ),
               tags$div(class = "card",
                        tags$img(src = "www/Camilo.jpg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Juan Camilo Henao-Rojas')),
                                 tags$img(src = 'www/agrosavia.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0003-0007-6809")
                        )
               ),
               tags$div(class = "card",
                        tags$img(src = "www/Andres.jpg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Andrés Cortés-Vera')),
                                 tags$img(src = 'www/agrosavia.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0003-4178-0675")
                        )
               )

      ),
      tags$br(),
      tags$br(),
      tags$div(class = "cardsteam",
               tags$div(class = "card",
                        tags$img(src = "www/Isabel.jpeg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Maria Isabel Chacon-Sanchez')),
                                 tags$img(src = 'www/unalblack.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0001-7781-6129")
                        )
               ),
               tags$div(class = "card",
                        tags$img(src = "www/Diego.jpg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Diego Felipe Conejo')),
                                 tags$img(src = 'www/ciat.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0001-7129-4016")
                        )
               ),
               tags$div(class = "card",
                        tags$img(src = "www/Marlon.jpeg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Marlon E. Cobos')),
                                 tags$img(src = 'www/kansas.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0002-2611-1767")
                        )
               )
      ),
      tags$br(),
      tags$br(),
      tags$div(class = "cardsteam",
               tags$div(class = "card",
                        tags$img(src = "www/Luis.jpeg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Luis Felipe Lopez')),
                                 tags$img(src = 'www/agrosavia.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0002-4967-6955")
                        )
               ),
               tags$div(class = "card",
                        tags$img(src = "www/Paula.jpg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Paula Helena Reyes Herrera')),
                                 tags$img(src = 'www/agrosavia.png', width = '100px'),
                                 tags$h6("ORCID:")
                        )
               ),
               tags$div(class = "card",
                        tags$img(src = "www/Kevin.jpg", width = "100px"),
                        tags$div(class = "containerteam",
                                 tags$h4(tags$b('Kevin Steven Quiroga Benavides')),
                                 tags$img(src = 'www/unalblack.png', width = '100px'),
                                 tags$h6("ORCID: https://orcid.org/0000-0002-2748-942X")
                        )
               )
      ),
      tags$br(),
      tags$div(class = "footer",
               tags$div(class = 'footer_1',
                        tags$h5(class="title", "Versión 3.0.1"),
                        tags$p(class = "footer_2", "© Copyright GermVersity. All rights reserved")))
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
