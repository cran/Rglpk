## Experimental thus not exported yet (submitted by Florian Schwendinger)
## Simple linear program.
## maximize:   2 x_1 + 4 x_2 + 3 x_3
## subject to: 3 x_1 + 4 x_2 + 2 x_3 <= 60
##             2 x_1 +   x_2 + 2 x_3 <= 40
##               x_1 + 3 x_2 + 2 x_3 <= 80
##               x_1, x_2, x_3 are non-negative real numbers

## obj <- c(2, 4, 3)
## mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
## dir <- c("<=", "<=", "<=")
## rhs <- c(60, 40, 80)
## max <- TRUE

## Rglpk_write_file("lp1_mps_fixed.mps", "MPS_fixed", obj, mat, dir, rhs, max = max)
## Rglpk_write_file("lp1_mps_free.mps", "MPS_free", obj, mat, dir, rhs, max = max)
## Rglpk_write_file("lp1_CPLEX_LP.mps", "CPLEX_LP", obj, mat, dir, rhs, max = max)
## Rglpk_write_file("lp1_MathProg.mps", "MathProg", obj, mat, dir, rhs, max = max)

Rglpk_write_file <- function(file, type = c("MPS_fixed", "MPS_free", "CPLEX_LP", "MathProg"), obj, mat, dir, rhs, bounds = NULL, types = NULL, max = FALSE){
    file_type <- which(match.arg(type) == c("MPS_fixed", "MPS_free", "CPLEX_LP", "MathProg"))

    out <- Rglpk_call( obj = obj, mat = mat, dir = dir, rhs = rhs, bounds = bounds, types = types, max = max,
                       canonicalize_status = logical(1L), presolve = logical(1L), time_limit = integer(1L), verb = logical(1L), ## default values should be ignored
                       file = file, file_type = file_type )
    invisible( out$status )
}


