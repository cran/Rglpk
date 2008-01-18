/* glplib05.c (terminal output and error handling) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 01, 02, 03, 04, 05, 06, 07 Andrew Makhorin,
*  Department for Applied Informatics, Moscow Aviation Institute,
*  Moscow, Russia. All rights reserved. E-mail: <mao@mai2.rcnet.ru>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "glplib.h"

/***********************************************************************
*  NAME
*
*  xputs - write character string to the terminal
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  void xputs(const char *s);
*
*  DESCRIPTION
*
*  The routine xputs writes the character string s to the terminal.
*  It does not write the terminating '\0' character and does not append
*  the '\n' character.
*
*  The routine xputs is a base low-level routine used by all other glpk
*  routines to perform terminal output. */

void xputs(const char *s)
{     LIBENV *env = lib_link_env();
      /* pass the string to the user-defined routine */
      if (env->term_hook != NULL)
         if (env->term_hook(env->term_info, s) != 0) goto skip;
      /* write the string to the terminal */
      if (env->term_out) fputs(s, stdout);
      /* write the string to the hardcopy file */
      if (env->log_file != NULL) fputs(s, env->log_file);
skip: return;
}

/***********************************************************************
*  NAME
*
*  xprintf - write formatted output to the terminal
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  void xprintf(const char *fmt, ...);
*
*  DESCRIPTION
*
*  The routine xprintf uses the format control string fmt and optional
*  parameter list to perform fomatting in the same way as the standard
*  function printf does and then writes the formatted output to the
*  terminal using the routine xputs. The number of characters written at
*  once must not exceed 4,000. */

static void xvprintf(const char *fmt, va_list arg)
{     char buf[4000+1];
#ifdef HAVE_VSNPRINTF
      vsnprintf(buf, sizeof(buf), fmt, arg);
#else
      vsprintf(buf, fmt, arg);
      xassert(strlen(buf) < sizeof(buf));
#endif
      xputs(buf);
      return;
}

void xprintf(const char *fmt, ...)
{     va_list arg;
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      return;
}

void xprint1(const char *fmt, ...)
{     /* (obsolete) */
      va_list arg;
      va_start(arg, fmt);
      xvprintf(fmt, arg), xputs("\n");
      va_end(arg);
      return;
}

/***********************************************************************
*  NAME
*
*  lib_term_hook - install hook to intercept terminal output
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  void lib_term_hook(int (*func)(void *info, const char *s),
*     void *info);
*
*  DESCRIPTION
*
*  The routine lib_term_hook installs the user-defined hook routine to
*  intercept all terminal output performed by glpk routines.
*
*  This feature can be used to redirect the terminal output to other
*  destination, for example to a file or a text window.
*
*  The parameter func specifies the user-defined hook routine. It is
*  called from an internal printing routine, which passes to it two
*  parameters: info and s. The parameter info is a transit pointer,
*  specified in the corresponding call to the routine lib_term_hook;
*  it may be used to pass some information to the hook routine. The
*  parameter s is a pointer to the null terminated character string,
*  which is intended to be written to the terminal. If the hook routine
*  returns zero, the printing routine writes the string s to the
*  terminal in a usual way; otherwise, if the hook routine returns
*  non-zero, no terminal output is performed.
*
*  To uninstall the hook routine the parameters func and info should be
*  specified as NULL. */

void lib_term_hook(int (*func)(void *info, const char *s), void *info)
{     LIBENV *env = lib_link_env();
      if (func == NULL)
      {  env->term_hook = NULL;
         env->term_info = NULL;
      }
      else
      {  env->term_hook = func;
         env->term_info = info;
      }
      return;
}

void lib_print_hook(int (*func)(void *info, char *buf), void *info)
{     /* (obsolete) */
#if 0 /* iso c complains */
      int (*hook)(void *, const char *) = (void *)func;
#else
      int (*hook)(void *, const char *) = (int(*)(void *, const char *))
         (func);
#endif
      lib_term_hook(hook, info);
      return;
}

/***********************************************************************
*  NAME
*
*  xfault - display error message and terminate execution
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  void xfault(const char *fmt, ...);
*
*  DESCRIPTION
*
*  The routine xfault uses the format control string fmt and optional
*  parameter list to perform formatting in the same way as the standard
*  function printf does, writes the formatted error message to the
*  terminal using the routine xputs, and then abnormally terminates the
*  program. The message length must not exceed 4,000 characters. */

void xfault(const char *fmt, ...)
{     va_list arg;
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      abort();
      /* no return */
}

void xfault1(const char *fmt, ...)
{     /* (obsolete) */
      va_list arg;
      va_start(arg, fmt);
      xvprintf(fmt, arg), xputs("\n");
      va_end(arg);
      abort();
      /* no return */
}

void lib_fault_hook(int (*func)(void *info, char *buf), void *info)
{     /* (obsolete) */
      xassert(func == func);
      xassert(info == info);
      return;
}

/***********************************************************************
*  NAME
*
*  xassert - check for logical condition
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  void xassert(int expr);
*
*  DESCRIPTION
*
*  The routine xassert (implemented as a macro) checks for a logical
*  condition specified by the parameter expr. If the condition is false
*  (i.e. the value of expr is zero), the routine displays an error
*  message and abnormally terminates the program. */

void _xassert(const char *expr, const char *file, int line)
{     xfault("GLPK internal error: %s; file %s, line %d\n",
         expr, file, line);
      /* no return */
}

/* eof */
