/* This is the GLPK C Interface
 */

#include "Rglpk.h"
#include <stdio.h>
#include <setjmp.h>

// this is the solve function called from R
void R_glp_solve (int *lp_direction, int *lp_number_of_constraints,
                  int *lp_direction_of_constraints, double *lp_right_hand_side,
                  int *lp_number_of_objective_vars,
                  double *lp_objective_coefficients,
                  int *lp_objective_var_is_integer, 
                  int *lp_objective_var_is_binary,
                  int *lp_is_integer,                     //should be boolean
                  int *lp_number_of_values_in_constraint_matrix,
                  int *lp_constraint_matrix_i, int *lp_constraint_matrix_j,
                  double *lp_constraint_matrix_values,
                  int *lp_bounds_type, double *lp_bounds_lower,
                  double *lp_bounds_upper,
                  double *lp_optimum,
                  int *lp_col_stat,
                  double *lp_objective_vars_values,
                  double *lp_objective_dual_values,
                  int *lp_row_stat,
                  double *lp_row_prim_aux,
                  double *lp_row_dual_aux,
                  int *lp_verbosity,
                  int *lp_presolve,
                  int *lp_time_limit,
                  int *lp_status,
                  int *write_fmt,
                  char **fname,
                  int *write_sensitivity_report,
                  char **fname_sensitivity_report) {

  // GLPK problem object
  glp_prob *lp;
  // GLPK simplex control object
  glp_smcp control_sm;
  // GLPK mixed integer control object
  glp_iocp control_io;
  int i, kl, ku;
  jmp_buf env;

  // Patch provided by Xypron: A far jump is used to return if an
  // error occurs. Prior to that R crashed.
  if (setjmp(env)) {
    error("An error occured inside the GLPK library.");
  } else {
    glp_error_hook(Rglpk_error_hook, &env);

    // create problem object 
    lp = glp_create_prob();

    // Turn on/off Terminal Output
    if(*lp_verbosity==1)
      glp_term_out(GLP_ON);
    else
      glp_term_out(GLP_OFF);
    
    // direction of optimization
    if(*lp_direction==1)
      glp_set_obj_dir(lp, GLP_MAX);
    else
      glp_set_obj_dir(lp, GLP_MIN);
    
    // is it a mixed integer problem? -- seems to be an R glpk function
    //if(lp_integer)
    //lpx_set_class(lp, LPX_MIP);
    // add rows to the problem object
    if( *lp_number_of_constraints > 0 ){
      glp_add_rows(lp, *lp_number_of_constraints);
      for(i = 0; i < *lp_number_of_constraints; i++)
        switch(lp_direction_of_constraints[i]){
        case 1: 
          glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, lp_right_hand_side[i]);
          break;
        case 2: 
          glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, lp_right_hand_side[i]);
          break;
        case 3: 
          glp_set_row_bnds(lp, i+1, GLP_LO, lp_right_hand_side[i], 0.0);
          break;
        case 4: 
          glp_set_row_bnds(lp, i+1, GLP_LO, lp_right_hand_side[i], 0.0);
          break;
        case 5: 
          glp_set_row_bnds(lp, i+1, GLP_FX, lp_right_hand_side[i],
                           lp_right_hand_side[i]);
          break;
        }
    }
    
    // add columns to the problem object
    glp_add_cols(lp, *lp_number_of_objective_vars);
    kl = ku = 0;
    for(i = 0; i < *lp_number_of_objective_vars; i++) {
      glp_set_col_bnds(lp, i+1, lp_bounds_type[i], lp_bounds_lower[i], lp_bounds_upper[i]);
      // set objective coefficients and integer if necessary
      glp_set_obj_coef(lp, i+1, lp_objective_coefficients[i]);
      if (lp_objective_var_is_integer[i])
        glp_set_col_kind(lp, i+1, GLP_IV);
      if (lp_objective_var_is_binary[i])
        glp_set_col_kind(lp, i+1, GLP_BV);
    }
    // load the matrix
    // IMPORTANT: as glp_load_matrix requires triplets as vectors of the
    // form: ia[1] ... ia[n], we have to pass the pointer to the adress
    // [-1] of the corresponding vector 
    if( *lp_number_of_constraints > 0 ){
      glp_load_matrix(lp, *lp_number_of_values_in_constraint_matrix,
                      &lp_constraint_matrix_i[-1],
                      &lp_constraint_matrix_j[-1], &lp_constraint_matrix_values[-1]);
    }

    // write lp to file
    // mps_fixed := 1L, mps_free := 2L
    if ( *write_fmt > 0 ) {
        const char *out_name = fname[0];
        if ( *write_fmt < 3 ) {
            *lp_status = glp_write_mps(lp, *write_fmt, NULL, out_name);
        } else if ( *write_fmt == 3 ) {
            *lp_status = glp_write_lp(lp, NULL, out_name);
        } else {
            int future_flag = 0;
            *lp_status = glp_write_prob(lp, future_flag, out_name);
        }
        glp_delete_prob(lp);
        return;
    }
    
    // set optimizer control parameters
    glp_init_smcp(&control_sm);
    if (*lp_time_limit > 0) {
      control_sm.tm_lim = *lp_time_limit;
    }
    if (*lp_presolve == 1) {
      control_sm.presolve = GLP_ON;
    }
    
    // run simplex method to solve linear problem
    glp_simplex(lp, &control_sm);
    
    // retrieve status of optimization
    *lp_status = glp_get_status(lp);
    // retrieve optimum
    *lp_optimum = glp_get_obj_val(lp);
    // retrieve values of objective vars
    for(i = 0; i < *lp_number_of_objective_vars; i++) {
      lp_col_stat[i] = glp_get_col_stat(lp, i+1);
      lp_objective_vars_values[i] = glp_get_col_prim(lp, i+1);
      lp_objective_dual_values[i] = glp_get_col_dual(lp, i+1);
    }
    // retrieve primal/dual multipliers
    for(i = 0; i < *lp_number_of_constraints; i++) {
      lp_row_stat[i] = glp_get_row_stat(lp, i+1);
      lp_row_prim_aux[i] = glp_get_row_prim(lp, i+1);
      lp_row_dual_aux[i] = glp_get_row_dual(lp, i+1);
    }
    if(*lp_is_integer) {
      // set optimizer control parameters
      glp_init_iocp(&control_io);
      if (*lp_time_limit > 0) {
        control_io.tm_lim = *lp_time_limit;
      }
      if (*lp_presolve == 1) {
        control_io.presolve = GLP_ON;
      }
      // optimize
      glp_intopt(lp, &control_io);
      // retrieve status of optimization
      *lp_status = glp_mip_status(lp);
      
      // retrieve MIP optimum
      *lp_optimum = glp_mip_obj_val(lp);
      // retrieve MIP values of objective vars
      for(i = 0; i < *lp_number_of_objective_vars; i++){
        lp_objective_vars_values[i] = glp_mip_col_val(lp, i+1);
      }
      // retrieve MIP auxiliary variable values
      for(i = 0; i < *lp_number_of_constraints; i++) {
        lp_row_prim_aux[i] = glp_mip_row_val(lp, i+1);
      }
    }

    // write sensitivity analysis report
    if (*write_sensitivity_report == 1) {
      const char *out_name = fname_sensitivity_report[0];
      glp_print_ranges(lp, 0, NULL, 0, out_name);
    }

    // delete problem object
    glp_delete_prob(lp);
  }
}


