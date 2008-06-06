\name{Rglpk_read_file}
\alias{Rglpk_read_file}
\alias{print.MP_data_from_file}
\title{Interface to GLPK's file reader}
\description{
  High level R interface to the MPS and CPLEX\_LP reader of the GNU Linear
  Programming Kit (GLPK).
}
\usage{
## MPS file reader
Rglpk_read_file(file, type = c("MPS_fixed", "MPS_free", "CPLEX_LP"),
ignore_first_row = FALSE, verbose = FALSE) 

## print method
\method{print}{MP_data_from_file}(x, \ldots)
}
\arguments{
  \item{file}{a character specifying the relative or absolute path to
    the model file.}
  \item{type}{a character specifying the file format. This can be either
    \code{"MPS_fixed"}, \code{"MPS_free"} and \code{"CPLEX_LP"}.}
  \item{ignore_first_row}{a logical indicating whether the first row of
    the model file should be ignored or not.
    Default: \code{FALSE}.}
  \item{verbose}{a logical for turning on/off additional solver output.
    Default: \code{FALSE}.}
  \item{x}{an object of class \code{"MP_data_from_file"}.}
  \item{\ldots}{further arguments passed on to the print method.}
}
\details{\code{Rglpk_read_file()} takes the path to a file as an
  argument and calls GLPK's file reader. The description of the linear or
  mixed integer linear program is returned as an object of class
  \code{"MP_data_from_file"}.
}
\value{
  \code{Rglpk_read_file()} returns the specification of a mixed integer
  linear program defined in \code{file} as an object of class
  \code{"MP_data_from_file"}.
}
%%\references{}
\author{Stefan Theussl}
%%\seealso{}
\examples{
## not run:
## x <- Rglpk_read_file(some_file.MPS)
## x
## Rglpk_solve_LP(x$objective, x$constraints[[1]], x$constraints[[2]],
##                x$constraints[[3]], x$types, x$maximum, x$bounds)
}
\keyword{IO}
\keyword{utilities}
