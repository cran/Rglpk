## Reads linear programs from MPS files
## uses GLPK's facilities for reading these files
## Interface to GLPK's MPS file reader

## Input: Path to a file specifying a mathematical program (MP),
##        the model specification language


## Output: an object of class 'MP_data_from_file' describing the MP 

## Mathematical Programming (MP) Data Object
##$file                     ... absolute path to data file 
##$type                     ... file type (currently 'MPS-fixed', 'MPS-free', 'CPLEX LP'
##$objective_coefficients   ... a vector of the objective coefficients
##$constraint_matrix        ... specifies the constraint matrixin simple triplet form
##$direction_of_constraints ... contains the direction of constraints
##$right_hand_side          ... vector of right hand side values
##$objective_var_is_integer ... a vector of logicals specifying which objective variable is of type 'integer'
##$objective_var_is_binary  ... a vector of logicals specifying which objective variable is of type 'binary'
##$direction_of_optimization... can be either 'min' or 'max'
##$bounds                   ... upper and lower bounds of objective variables
##$n_objective_vars         ... number of objective variables    
##$n_integer_vars           ... number of variables which are of type 'integer'
##$n_binary_vars            ... number of variables which are of type 'binary'
##$n_constraints            ... number of constraints    
##$n_values_in_constraint_matrix ... number of values in constraint matrix
##$problem_name             ... name of the problem     
##


Rglpk_read_file <- function(file, type = c("MPS_fixed", "MPS_free", "CPLEX_LP"), ignore_first_row = FALSE, verbose = FALSE){
  if(!file.exists(file))
    stop(paste("There is no file called", file, "!"))
  ## which file type to read from
  type <- match.arg(type)
  type_db <- c("MPS_fixed" = 1L,
               "MPS_free"  = 2L,
               "CPLEX_LP"  = 3L
               )
  type <- type_db[type]
  obj <- list(file = tools::file_path_as_absolute(file),
              type = type)
  ## read files in a two step approach 
  ## first, retrieve meta data like number of objective variables etc.
  ## we need this to reserve memory in R accordingly
  meta_data <- glp_get_meta_data_from_file(obj, verbose)
  ## second, read all remaining data
  milp_data <- glp_retrieve_MP_from_file(meta_data, ignore_first_row, verbose)
  ## merge everything together
  MP_data <- glp_merge_MP_data(meta_data, milp_data)
  ## Post processing
  MP_data$type <- names(type_db[type_db == MP_data$type])
  ## unbounded (FREE, 1L) variable gets an upper bound of Inf
  ## and a direction '<='
  dir_db <- c("<=" = 1L, ">=" = 2L, "<=" = 3L, "DB" = 4L, "==" = 5L)
  MP_data$direction_of_constraints <- names(dir_db[MP_data$direction_of_constraints])
  ## default is to have only continuous variables
  ## if any is binary or integer set the value accordingly
  types <- rep("C", length.out = MP_data$n_objective_vars)
  if(any(MP_data$objective_var_is_integer)){
    types[MP_data$objective_var_is_integer] <- "I"
  }
  if(any(MP_data$objective_var_is_binary)){
    types[MP_data$objective_var_is_binary] <- "B"
  }
  ## build object we want to return
  ## First add MILP to the object
  out <- MILP(objective = MP_data$objective_coefficients,
              constraints = list(MP_data$constraint_matrix,
                                 MP_data$direction_of_constraints,
                                 MP_data$right_hand_side),
              bounds = MP_data$bounds,
              types = types,
              maximum = MP_data$maximize
              )
  out$n_objective_vars <- MP_data$n_objective_vars
  out$n_integer_vars <- MP_data$n_integer_vars
  out$n_binary_vars <- MP_data$n_binary_vars
  out$n_constraints <- MP_data$n_constraints
  ##out$n_values_in_constraint_matrix <- MP_data$n_values_in_constraint_matrix
  out$file_type <- MP_data$type
  out$file_name <- MP_data$file
  class(out) <- c("MP_data_from_file", class(out))
  out
}

## First parse file to get some meta data of the LP/MILP
## (number of constraints/objective variables, direction of optimization, ...)
glp_get_meta_data_from_file <- function(x, verbose){   
  res <- .C("Rglpk_read_file",
            file                          = as.character(x$file),
            type                          = as.integer(x$type),
            direction_of_optimization     = integer(1L),
            n_constraints                 = integer(1L),         
            n_objective_vars              = integer(1L),
            n_values_in_constraint_matrix = integer(1L),
            n_integer_vars                = integer(1L),
            n_binary_vars                 = integer(1L),
            verbosity                     = as.integer(verbose),
            PACKAGE = "Rglpk")
  res
}

## Retrieve all missing elements of the LP/MILP
glp_retrieve_MP_from_file <- function(x, ignore_first_row, verbose = FALSE){
  res <- .C("Rglpk_retrieve_MP_from_file",
            file                     = as.character(x$file),
            type                     = as.integer(x$type),
            n_constraints            = x$n_constraints,         
            n_objective_vars         = x$n_objective_vars,
            ##n_values_in_constraint_matrix = x$n_values_in_constraint_matrix,
            ##n_integer_vars           = x$n_integer_vars,
            ##n_binary_vars            = x$n_binary_vars,
            objective_coefficients   = double(x$n_objective_vars),
            constraint_matrix_i      = integer(x$n_values_in_constraint_matrix),
            constraint_matrix_j      = integer(x$n_values_in_constraint_matrix),
            constraint_matrix_values = double(x$n_values_in_constraint_matrix),
            direction_of_constraints = integer(x$n_constraints),
            right_hand_side          = double(x$n_constraints),           
            objective_var_is_integer = integer(x$n_objective_vars),
            objective_var_is_binary  = integer(x$n_objective_vars),
            bounds_type              = integer(x$n_objective_vars),
            bounds_lower             = double(x$n_objective_vars),
            bounds_upper             = double(x$n_objective_vars),
            lp_ignore_first_row      = as.integer(ignore_first_row),
            verbosity                = as.integer(verbose),
            PACKAGE = "Rglpk")
  ## lp_is_integer               = as.integer(lp_is_integer),

  ## replace infinity values
  res$bounds_lower <- replace(res$bounds_lower, res$bounds_lower == -.Machine$double.xmax, -Inf)
  res$bounds_upper <- replace(res$bounds_upper, res$bounds_upper == .Machine$double.xmax, Inf)
  ## in MPS definition first row is sometimes problematic. E.g., in MIPLIB2003
  ## it has to be removed!
  if(ignore_first_row){
    res$n_constraints <- res$n_constraints - 1
    ## zeros values in the constraint matrix have to be removed, these
    ## are the values from the first row
    to_remove <- which(res$constraint_matrix_values == 0)
    res$constraint_matrix_i <- res$constraint_matrix_i[-to_remove] - 1
    res$constraint_matrix_j <- res$constraint_matrix_j[-to_remove]
    res$constraint_matrix_values <- res$constraint_matrix_values[-to_remove]
    res$right_hand_side <- res$right_hand_side[-1]
    #res$direction_of_constraints <- res$direction_of_constraints[-length(res$right_hand_side)]
  }
  res
}
                        
glp_merge_MP_data <- function(x, y){
  out <- list(objective_coefficients        = y$objective_coefficients,
              constraint_matrix             = simple_triplet_matrix(
                                                   y$constraint_matrix_i,
                                                   y$constraint_matrix_j,
                                                   y$constraint_matrix_values,
                                                   y$n_constraints,
                                                   y$n_objective_vars),
              direction_of_constraints      = y$direction_of_constraints,
              right_hand_side               = y$right_hand_side,
              objective_var_is_integer      = as.logical(y$objective_var_is_integer),
              objective_var_is_binary       = as.logical(y$objective_var_is_binary),
              ## minimization if GLP_MIN (1L) or max if GLP_MAX (2L)
              maximize                      = x$direction_of_optimization == 2L,
              bounds                        = list(lower = list(ind = 1L:x$n_objective_vars,
                                                                val = y$bounds_lower),
                                                   upper = list(ind = 1L:x$n_objective_vars,
                                                                val = y$bounds_upper)),
              n_objective_vars              = x$n_objective_vars,
              n_integer_vars                = x$n_integer_vars,
              n_binary_vars                 = x$n_binary_vars,
              ## here from y because it might have changed -> ignore_first_row_parameter
              n_constraints                 = y$n_constraints,
              n_values_in_constraint_matrix = x$n_values_in_constraint_matrix,
              ## problem_name                  = x$problem_name,
              file                          = x$file,
              type                          = x$type
              )
  out
}
