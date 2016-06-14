.check_control_parameters <- function( control ){
    stopifnot(is.list(control))
    ## Default paramenters (currently only verbosity supported)
    out <- list(verbose = FALSE,
                presolve = FALSE,
                tm_limit = 0L,
                canonicalize_status = TRUE)
    out[names(control)] <- control
    if (!is.null(out$verbose)) {
        out$verbose <- as.integer(out$verbose)
        if (!out$verbose %in% c(0L, 1L)) {
            warning("Improper value for 'verbose' parameter. Using default.")
            out$verbose <- 0L
        }
    }
    if (!is.null(out$presolve)) {
        out$presolve <- as.integer(out$presolve)
        if (!out$presolve %in% c(0L, 1L)) {
            warning("Improper value for 'presolve' parameter. Using default.")
            out$verbose <- 0L
        }
    }
    if (!is.null(out$canonicalize_status)) {
        out$canonicalize_status <- as.logical(out$canonicalize_status)
    }
    if( !is.null(control$tm_limit) )
        out$tm_limit <- as.integer(out$tm_limit)
    out
}
