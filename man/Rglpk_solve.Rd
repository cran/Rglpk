\name{Rglpk_solve_LP}
\alias{Rglpk_solve_LP}
\title{Linear and Mixed Integer Programming Solver Using GLPK}
\description{
  High level R interface to the GNU Linear Programming Kit for solving
  linear as well as mixed integer linear programming problems (MILPs).
}
\usage{
Rglpk_solve_LP(obj, mat, dir, rhs, types = NULL, max = FALSE,
               bounds = NULL, verbose = FALSE)
}
\arguments{
  \item{obj}{a vector with the objective coefficients}
  \item{mat}{a vector or a matrix of the constraint coefficients}
  \item{dir}{a character vector with the directions of the constraints.
    Each element must be one of \code{"<"}, \code{"<="}, \code{">"},
    \code{">="}, or \code{"=="}.}
  \item{rhs}{the right hand side of the constraints}
  \item{types}{a vector indicating the types of the objective
    variables. \code{types} can be either \code{"B"} for binary,
    \code{"C"} for continuous or \code{"I"} for integer. By default all
    variables are of type \code{"C"}.}
  \item{max}{a logical giving the direction of the optimization.
    \code{TRUE} means that the objective is to maximize the objective
    function, \code{FALSE} (default) means to minimize it.}
  \item{bounds}{\code{NULL} (default) or a list with elements
    \code{upper} and \code{lower} containing the indices and
    corresponding bounds of the objective variables.  The default for
    each variable is a bound between 0 and \code{Inf}.}
  \item{verbose}{a logical for turning on/off additional solver output:
    Default: \code{FALSE}.}
}
\details{
  The GNU Linear Programming Kit is open source.  The current version
  can be found at
  \url{http://www.gnu.org/software/glpk/glpk.html}.  Package \pkg{Rglpk}
  provides a high level solver function using the low level C interface
  of the GLPK solver.  There also exists an R interface done by Lopaka Lee
  which ports all low level C Interface routines of the GLPK API to R (R
  package \pkg{glpk}). 
}
\value{
  A list containing the optimal solution, with the following components.
  \item{solution}{the vector of optimal coefficients}
  \item{objval}{the value of the objective function at the optimum}
  \item{status}{an integer with status information about the solution
    returned: 0 if the optimal solution was found, a non-zero value
    otherwise.}
}
\references{
  GNU Linear Programming Kit
  (\url{http://www.gnu.org/software/glpk/glpk.html}).
  
  GLPK Interface to R 
  (\url{http://cran.R-project.org/package=Rglpk}).
}
\author{Stefan Theussl and Kurt Hornik}
\seealso{
  \code{\link[lpSolve]{lp}} in package \pkg{lpSolve};
  \code{\link[Rsymphony]{Rsymphony_solve_LP}} in package
  \pkg{Rsymphony}.
}
\examples{
## Simple linear program.
## maximize:   2 x_1 + 4 x_2 + 3 x_3
## subject to: 3 x_1 + 4 x_2 + 2 x_3 <= 60
##             2 x_1 +   x_2 +   x_3 <= 40
##               x_1 + 3 x_2 + 2 x_3 <= 80
##               x_1, x_2, x_3 are non-negative real numbers

obj <- c(2, 4, 3)
mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
dir <- c("<=", "<=", "<=")
rhs <- c(60, 40, 80)
max <- TRUE

Rglpk_solve_LP(obj, mat, dir, rhs, max = max)

## Simple mixed integer linear program.
## maximize:    3 x_1 + 1 x_2 + 3 x_3
## subject to: -1 x_1 + 2 x_2 +   x_3 <= 4
##                      4 x_2 - 3 x_3 <= 2
##                x_1 - 3 x_2 + 2 x_3 <= 3
##                x_1, x_3 are non-negative integers
##                x_2 is a non-negative real number

obj <- c(3, 1, 3)
mat <- matrix(c(-1, 0, 1, 2, 4, -3, 1, -3, 2), nrow = 3)
dir <- c("<=", "<=", "<=")
rhs <- c(4, 2, 3)
types <- c("I", "C", "I")
max <- TRUE

Rglpk_solve_LP(obj, mat, dir, rhs, types, max)

## Same as before but with bounds replaced by
## -Inf <  x_1 <= 4
##    0 <= x_2 <= 100
##    2 <= x_3 <  Inf

bounds <- list(lower = list(ind = c(1L, 3L), val = c(-Inf, 2)),
               upper = list(ind = c(1L, 2L), val = c(4, 100)))
Rglpk_solve_LP(obj, mat, dir, rhs, types, max, bounds)

## Examples from the GLPK manual

## 1.3.1
## maximize:   10 x_1 + 6 x_2 + 4 x_3
## subject to:    x_1 +   x_2 +   x_3 <= 100
##             10 x_1 + 4 x_2 + 5 x_3 <= 600
##              2 x_1 + 2 x_2 + 6 x_3 <= 300
##                x_1,    x_2,    x_3 are non-negative real numbers

obj <- c(10, 6, 4)
mat <- matrix(c(1, 10, 2, 1, 4, 2, 1, 5, 6), nrow = 3)
dir <- c("<=", "<=", "<=")
rhs <- c(100, 600, 300)
max <- TRUE

Rglpk_solve_LP(obj, mat, dir, rhs, max = max)


}
\keyword{optimize}
