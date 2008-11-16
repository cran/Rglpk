require(glpk)

## the R-ported GNU Linear Programming kit
## solve function
glp_solve <- function(x, ...){
  if (!inherits(x, "lp"))
        stop("'lp' must inherit from class \"lp\"")
  ## create problem object
  p.obj<-lpx_create_prob()
  ## add the optimization direction to the problem object
  if(x$direction=="max") lpx_set_obj_dir(p.obj, LPX_MAX)
  else lpx_set_obj_dir(p.obj, LPX_MIN)
  ## is it an integer problem?
  if(x$is_integer)
    lpx_set_class(p.obj, LPX_MIP)
  ## add rows to the problem object
  lpx_add_rows(p.obj, x$n_of_constraints)
  for(i in 1:x$n_of_constraints)
    switch(x$direction_of_constraints[i],
           "default"=stop("'direction_of_constraints' must be either '<', '<=', ..."),
           "<"  = lpx_set_row_bnds(p.obj, i, LPX_UP, 0.0, x$right_hand_side[i]),
           "<=" = lpx_set_row_bnds(p.obj, i, LPX_UP, 0.0, x$right_hand_side[i]),
           ">"  = lpx_set_row_bnds(p.obj, i, LPX_LO, x$right_hand_side[i], 0.0),
           ">=" = lpx_set_row_bnds(p.obj, i, LPX_LO, x$right_hand_side[i], 0.0),
           "==" = lpx_set_row_bnds(p.obj, i, LPX_FX, x$right_hand_side[i],
                                                     x$right_hand_side[i])
           )
  ## add columns to the problem object
  lpx_add_cols(p.obj, x$n_of_objective_vars)
  for(i in 1:x$n_of_objective_vars){
    lpx_set_col_bnds(p.obj, i, LPX_LO, 0.0, 0.0)
    lpx_set_obj_coef(p.obj, i, x$objective_coefficients[i])
    if(x$objective_var_is_integer[i])
      lpx_set_col_kind(p.obj, i, LPX_IV)
  }
  ## load the matrix
  lpx_load_matrix(p.obj, x$n_of_values_in_constraint_matrix, x$constraint_matrix_i,
                  x$constraint_matrix_j, x$constraint_matrix_v)

  ## run simplex algorithm
  lpx_simplex(p.obj)
  opt<-lpx_get_obj_val(p.obj)
  obj<-NULL
  ## rc<-NULL
  ## sp<-NULL
  for(i in 1:x$n_of_objective_vars){
    obj<-c(obj, lpx_get_col_prim(p.obj,i))
    ## rc<-c(rc,lpx_get_col_dual(lp,i))
  }
  ##for(i in 1:length(b))
  ##  sp<-c(sp,lpx_get_row_dual(lp,i))
  if(x$is_integer){
    lpx_integer(p.obj)
    opt<-lpx_mip_obj_val(p.obj)
    obj<-NULL
    ##rc<-NA
    ##sp<-NA
    for(i in 1:x$n_of_objective_vars){
    obj<-c(obj,lpx_mip_col_val(p.obj,i))
    }
  }
  lpx_delete_prob(p.obj)
  x$optimum <- opt
  x$solution <- obj
  x
}
