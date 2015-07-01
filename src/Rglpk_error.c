#include "Rglpk.h"
#include "setjmp.h"

/*
 * Path provided by Xypron
 * This hook function will be called if an error occured when
 * calling the glpk library
 */
void Rglpk_error_hook(void *in) {
  /* free glpk memory */
  glp_free_env();
  /* set print hook for terminal */
  Rglpk_initialize();
  /* safely return */
  longjmp(*((jmp_buf*)in), 1);
}
