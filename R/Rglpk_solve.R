## the R-ported GNU Linear Programming kit
## solve function --- C Interface

Rglpk_solve_LP <-
function(obj, mat, dir, rhs, types = NULL, max = FALSE,
         bounds = NULL, verbose = FALSE)
{
  ## validate direction of optimization
  if(!identical(max, TRUE) && !identical(max, FALSE))
    stop("'Argument 'max' must be either TRUE or FALSE.")
  direction_of_optimization <- as.integer(max)

  ## validate verbosity flag
  if(!identical(verbose, TRUE) && !identical(verbose, FALSE))
    stop("'Argument 'verbose' must be either TRUE or FALSE.")
  verb <- as.integer(verbose)
 
  ## match direction of constraints
  n_of_constraints <- length(dir)
  ## match relational operators to requested input
  direction_of_constraints <- match(dir, c("<", "<=", ">", ">=", "=="))
  if(any(is.na(direction_of_constraints)))
    stop("Argument 'dir' must be either '<', '<=', '>', '>=' or '=='.")
  
  n_of_objective_vars <- length(obj)

  constraint_matrix <- as.simple_triplet_matrix(mat)

  ## types of objective coefficients
  ## Default: "C"
  if(is.null(types))
    types <- "C"
  ## check if valid types
  if(any(is.na(match(types, c("I", "B", "C"), nomatch = NA))))
    stop("'types' must be either 'B', 'C' or 'I'.")
  ## replicate types to fit number of columns
  types <- rep(types, length.out = n_of_objective_vars)
  ## need a TRUE/FALSE integer/binary representation
  integers <- types == "I"
  binaries <- types == "B"
  
  ## do we have a mixed integer linear program?
  is_integer <- any(binaries | integers)

  ## bounds of objective coefficients
  bounds <- as.glp_bounds(as.list(bounds), n_of_objective_vars)

  ## call the C interface - this actually runs the solver
  x <- glp_call_interface(obj, n_of_objective_vars, constraint_matrix$i,
                          constraint_matrix$j, constraint_matrix$v,
                          length(constraint_matrix$v),
                          rhs, direction_of_constraints, n_of_constraints,
                          is_integer,
                          integers, binaries,
                          direction_of_optimization, bounds[,1L],
                          bounds[,2L], bounds[,3L], verb)
  out <- list(optimum=NA, solution=NA, status=NA)
  out$optimum <- x$lp_optimum
  out$solution <- x$lp_objective_vars_values
  ## match status of solution
  ## 0 -> optimal solution (5 in GLPK) else 1
  out$status <- as.integer(x$lp_status != 5L)
  out
}

## this function calls the C interface
glp_call_interface <-
function(lp_objective_coefficients, lp_n_of_objective_vars,
         lp_constraint_matrix_i, lp_constraint_matrix_j, lp_constraint_matrix_v,
         lp_n_of_values_in_constraint_matrix, lp_right_hand_side,
         lp_direction_of_constraints, lp_n_of_constraints, lp_is_integer,
         lp_objective_var_is_integer, lp_objective_var_is_binary,
         lp_direction_of_optimization,
         lp_bounds_type, lp_bounds_lower, lp_bounds_upper,
         verbose)
{
  out <- .C("R_glp_solve",
            lp_direction_of_optimization= as.integer(lp_direction_of_optimization),
            lp_n_of_constraints         = as.integer(lp_n_of_constraints),
            lp_direction_of_constraints = as.integer(lp_direction_of_constraints),
            lp_right_hand_side          = as.double(lp_right_hand_side),
            lp_n_of_objective_vars      = as.integer(lp_n_of_objective_vars),
            lp_objective_coefficients   = as.double(lp_objective_coefficients),
            lp_objective_var_is_integer = as.integer(lp_objective_var_is_integer),
            lp_objective_var_is_binary  = as.integer(lp_objective_var_is_binary),
            lp_is_integer               = as.integer(lp_is_integer),
            lp_n_of_values_in_constraint_matrix = as.integer(lp_n_of_values_in_constraint_matrix),
            lp_constraint_matrix_i      = as.integer(lp_constraint_matrix_i),
            lp_constraint_matrix_j      = as.integer(lp_constraint_matrix_j),
            lp_constraint_matrix_values = as.double(lp_constraint_matrix_v),
            lp_bounds_type             = as.integer(lp_bounds_type),
            lp_bounds_lower             = as.double(lp_bounds_lower),
            ## lp_n_of_bounds_l            = as.integer(length(lp_lower_bounds_i)),
            lp_bounds_upper             = as.double(lp_bounds_upper), 
            ## lp_n_of_bounds_u            = as.integer(length(lp_upper_bounds_i)),
            lp_optimum                  = double(1),
            lp_objective_vars_values    = double(lp_n_of_objective_vars),
            lp_verbosity                = as.integer(verbose),
            lp_status                   = integer(1),
            NAOK = TRUE, PACKAGE = "Rglpk")
  out
}
