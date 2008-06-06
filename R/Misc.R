## Miscellaneous functions used in this package

## Get the right representation of integer vector.
## This function takes an index vector and returns a vector of logicals
## with length 'n'.

#glp_integers <-
#function(x, n)
#{
#  if(!all(x <= n))
#stop("Indices must not exceed the number of objective coefficients.")
#  out <- logical(n)
#  out[x] <- TRUE
#  out
#}

print.MP_data_from_file <- function(x, ...){
  if(!inherits(x, "MP_data_from_file"))
     stop("'x' must be of class 'MP_data_from_file'")
  if(x$n_integer_vars > 0){
    writeLines(paste("A mixed integer linear program with",
                     x$n_objective_vars, "objective variables,"))
    writeLines(paste(x$n_integer_vars, "are integer and",
                     x$n_binary_vars, "of which are binary variables."))
    writeLines(paste("This problem has", x$n_constraints,
                     "constraints with", length(x$constraints[[1]]$v),
                     "non-zero values in the constraint matrix."))
  } else{
    writeLines(paste("A linear program with", x$n_objective_vars, "objective variables."))
    writeLines(paste("This problem has", x$n_constraints,
                     "constraints with", dim(x$constraints)[1] ,
                     "non-zero values in the constraint matrix."))
  } 
}
