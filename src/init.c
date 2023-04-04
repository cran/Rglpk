#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

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
                  char **fname_sensitivity_report);
void Rglpk_initialize(void);
void Rglpk_get_engine_version(char **GLPK_version);
void R_glp_read_file (char **file, int *type, 
                      int *lp_direction_of_optimization,
                      int *lp_n_constraints, int *lp_n_objective_vars,
                      int *lp_n_values_in_constraint_matrix,
                      int *lp_n_integer_vars, int *lp_n_binary_vars, 
                      char **lp_prob_name,
                      char **lp_obj_name,
                      int *lp_verbosity);
void Rglpk_delete_prob(void);
void Rglpk_retrieve_MP_from_file (char **file, int *type,
                                  int *lp_n_constraints,
                                  int *lp_n_objective_vars,
                                  double *lp_objective_coefficients,
                                  int *lp_constraint_matrix_i,
                                  int *lp_constraint_matrix_j,
                                  double *lp_constraint_matrix_values,
                                  int *lp_direction_of_constraints,
                                  double *lp_right_hand_side,
                                  double *lp_left_hand_side,
                                  int *lp_objective_var_is_integer,
                                  int *lp_objective_var_is_binary,
                                  int *lp_bounds_type,
                                  double *lp_bounds_lower,
                                  double *lp_bounds_upper,
                                  int *lp_ignore_first_row,
                                  int *lp_verbosity,
                                  char **lp_constraint_names,
                                  char **lp_objective_vars_names
                                  );

static const R_CMethodDef CEntries[] = {
    {"R_glp_solve", (DL_FUNC) &R_glp_solve, 31},
    {"Rglpk_initialize", (DL_FUNC) &Rglpk_initialize, 0},
    {"Rglpk_get_engine_version", (DL_FUNC) &Rglpk_get_engine_version, 1},
    {"R_glp_read_file", (DL_FUNC) &R_glp_read_file, 11},
    {"Rglpk_delete_prob", (DL_FUNC) &Rglpk_delete_prob, 0},
    {"Rglpk_retrieve_MP_from_file", (DL_FUNC) &Rglpk_retrieve_MP_from_file, 20},
    {NULL, NULL, 0}
};

void R_init_Rglpk(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
