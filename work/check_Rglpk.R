library("Rglpk")
verbose = TRUE

## check if Rglpk can solve simple LPs
obj <- c(2, 4, 3)
mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
dir <- c("<=", "<=", "<=")
rhs <- c(60, 40, 80)
max <- TRUE

Rglpk_solve_LP(obj, mat, dir, rhs, max = max)


## check if Rglpk can solve simple MILPs
obj <- c(3, 1, 3)
mat <- matrix(c(-1, 0, 1, 2, 4, -3, 1, -3, 2), nrow = 3)
dir <- c("<=", "<=", "<=")
rhs <- c(4, 2, 3)
types <- c("I", "C", "I")
max <- TRUE

Rglpk_solve_LP(obj, mat, dir, rhs, types, max)


bounds <- list(lower = list(ind = c(1L, 3L), val = c(-Inf, 2)),
               upper = list(ind = c(1L, 2L), val = c(4, 100)))
Rglpk_solve_LP(obj, mat, dir, rhs, types, max, bounds)


## check if Rglpk still reads file in MPS format
MPS <- Rglpk_read_file(file = "./misc07.mps", type = "MPS_free",
                       ignore_first_row = TRUE, verbose = verbose)
MPS

## solve problem
Rglpk:::.Rglpk_solve(MPS, control = list(verbose = verbose))


## check if Rglpk still reads file in CPLEX/LP format
CPLEX <- Rglpk_read_file(file = "./u5.cplx", type = "CPLEX_LP", verbose = verbose)
CPLEX

## solve problem
Rglpk:::.Rglpk_solve(CPLEX, control = list(verbose = verbose))
