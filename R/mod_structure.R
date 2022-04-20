#' structure UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_structure_ui <- function(id){
  ns <- NS(id)
  tagList(
    shinydashboard::dashboardPage(
      skin = 'purple',
      shinydashboard::dashboardHeader(
        title = 'GermVersity',
        titleWidth = 300,
        shinydashboard::dropdownMenu(type = 'messages',
                     icon = icon('github'),
                     shinydashboard::messageItem(
                       from = 'Github',
                       message = '',
                       icon = icon('github'),
                       href = 'https://github.com/Viinky-Kevs/microsoftAI'
                     ))
      ),
      shinydashboard::dashboardSidebar(
        width = 300,
        shinydashboard::sidebarMenu(
          id = 'sidebar',
          style = "position: relative; overflow: visible;",
          tags$div(
            tags$img(height = "280px",
                     alt="Logo",
                     style = 'display: flex; justify-content: center',
                     src="www/Logo.png")
          ),
          shinydashboard::menuItem("Home",
                                   tabName = 'home',
                                   icon = icon('home')),
          shinydashboard::menuItem('Module 1',
                                   tabName = 'module1',
                                   icon = icon('bezier-curve'),
                                   #shinydashboard::menuSubItem('Upload files',
                                    #           tabName = 'upload',
                                    #           icon = icon('upload')),
                                   shinydashboard::menuSubItem('Genotypic module',
                                            tabName = 'genotypic1',
                                            icon = icon('dna')),
                                   shinydashboard::menuSubItem('Phenotypic data',
                                               tabName = 'phenotypic',
                                               icon = icon('leaf')),
                                   shinydashboard::menuSubItem('Descriptor optimization',
                                               tabName = 'optimization',
                                               icon = icon('chart-bar'))
          ),
          shinydashboard::menuItem('Module 2',
                                   tabName = 'module2',
                                   icon = icon('bezier-curve'),
                                   shinydashboard::menuSubItem('Spatial DM',
                                                               tabName = 'spatialdm',
                                                               icon = icon('globe'))
                                   #,icon = icon('globe'))
          ),
          shinydashboard::menuItem('Module 3',
                                   tabName = 'module3',
                                   icon = icon('bezier-curve'),
                                   shinydashboard::menuSubItem('GAPIT',
                                                               tabName = 'gapit',
                                                               icon = icon('archway')),
                                   shinydashboard::menuSubItem('CPS',
                                                               tabName = 'cps',
                                                               icon = icon('bell')),
                                   shinydashboard::menuSubItem('LFMM',
                                                              tabName = 'lfmm',
                                                              icon = icon('brain'))
                                   ),
          shinydashboard::menuItem('Team',
                                   tabName = 'team',
                                   icon = icon('users')),
          HTML(paste0(
            "<br><br><br><br><br><br><br><br><br>",
            "<br><br><br><br><br><br><br><br><br>",
            "<br><br><br><br><br><br><br><br><br>",
            "<script>",
            "var today = new Date();",
            "var yyyy = today.getFullYear();",
            "</script>",
            "<p style = 'text-align: center;'><small>&copy; - <a href='#' target='_blank'>GermVersity </a> - <script>document.write(yyyy);</script></small></p>")
          )
        )

      ),
      shinydashboard::dashboardBody(
        shinydashboard::tabItems(
          shinydashboard::tabItem(
            tabName = 'home',
            fluidPage(
              mod_home_page_ui("home_page_1")
            )
          ),
          shinydashboard::tabItem(
            tabName = 'genotypic1',
            fluidPage(
              mod_genotypic_ui("genotypic_1")
            )
          ),
          shinydashboard::tabItem(
            tabName = 'phenotypic',
            fluidPage(
              mod_phenotypic_ui("phenotypic_1")
            )
          ),
          shinydashboard::tabItem(
            tabName = 'optimization',
            fluidPage(
              mod_descriptor_ui("descriptor_1")
            )
          ),
          shinydashboard::tabItem(
            tabName = 'team',
            fluidPage(
              mod_team_ui("team_1")
            )
          ),
          shinydashboard::tabItem(
            tabName = 'gapit',
            fluidPage(
              mod_gapit_ui("gapit_1")
            )
          ),
          shinydashboard::tabItem(
            tabName = 'spatialdm',
            fluidPage(
            mod_SDM_ui("SDM_1")
           )
          ),
          shinydashboard::tabItem(
            tabName = 'cps',
            fluidPage(
              mod_cps_ui("cps_1")
            )
          ),
          shinydashboard::tabItem(
            tabName = 'lfmm',
            fluidPage(
              mod_lfmm_ui("lfmm_1")
            )
          )
        )
      )
    )
  )
}

#' structure Server Functions
#'
#' @noRd
mod_structure_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
  })
}

## To be copied in the UI
# mod_structure_ui("structure_1")

## To be copied in the server
# mod_structure_server("structure_1")
