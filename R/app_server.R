#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  ## Call modules
  mod_genotypic_server("genotypic_1")
  mod_phenotypic_server("phenotypic_1")
  mod_SDM_server("SDM_1")
  mod_gapit_server("gapit_1")
  mod_cps_server("cps_1")
  mod_lfmm_server("lfmm_1")

  ## Verify function
  print('Backend succesfully work!')
}
