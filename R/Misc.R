## Miscellaneous functions used in this package

## Get the right representation of integer vector.
## This function takes an index vector and returns a vector of logicals
## with length 'n'.

glp_integers <-
function(x, n)
{
  if(!all(x <= n))
    stop("Indices must not exceed the number of objective coefficients.")
  out <- logical(n)
  out[x] <- TRUE
  out
}
