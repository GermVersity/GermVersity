#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  ## Call modules
  mod_genotypic_server("genotypic_1")
  mod_genotypic1_server("genotypic1_1")
  mod_genotypic2_server("genotypic2_1")
  mod_genotypic3_server("genotypic3_1")
  mod_genotypic4_server("genotypic4_1")
  mod_genotypic5_server("genotypic5_1")
  mod_genotypic6_server("genotypic6_1")
  mod_descriptor_server("descriptor_1")
  mod_phenotypic_server("phenotypic_1")
  mod_SDM_server("SDM_1")
  mod_gapit_server("gapit_1")
  mod_cps_server("cps_1")
  mod_lfmm_server("lfmm_1")
  mod_lfmm1_server("lfmm1_1")
  mod_lfmm2_server("lfmm2_1")

  ## Verify function
  print('Backend succesfully work!')
}
