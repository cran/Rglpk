.onLoad <- function(lib, pkg){
  ## initialize print routines
  .C("Rglpk_initialize", PACKAGE = pkg)
  out <- .C("Rglpk_get_engine_version", GLPK_version = character(1L), PACKAGE = pkg)
  writeLines(paste("Using the GLPK callable library version", out$GLPK_version))
}
