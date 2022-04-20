#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options

run_GermVersity <- function(...) {
  with_golem_options(app = shinyApp(ui = app_ui,
                                    server = app_server,
                                    options = list(launch.browser = TRUE)),
                     golem_opts = list(...))
}
