/* glplpx20.c */

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

#include "glpapi.h"
#include "glpgmp.h"
#include "glpmpl.h"
#define print xprint1
#define fault old_fault

/*--------------------------------------------------------------------*/
/* This program is a stand-alone LP/MIP solver. For pure LP problems  */
/* either the simplex method or the primal-dual interior point method */
/* can be used. For MIP problems the branch-and-bound procedure based */
/* on the simplex method is used.                                     */
/*--------------------------------------------------------------------*/

static int format = 4;
/* type of input text file:
   0 - fixed MPS
   1 - CPLEX LP
   2 - GNU MathProg
   3 - GNU LP
   4 - free MPS */

static const char *in_file = NULL;
/* name of input text file */

static const char *in_data = NULL;
/* name of optional input text file, which contains data section; NULL
   means no separate data section is provided */

static const char *display = NULL;
/* name of optional output text file, to which display output is sent;
   NULL means the output is sent to stdout */

static const char *in_bas = NULL;
/* name of input text file, which contains initial LP basis in MPS
   format; NULL means no initial LP basis is provided */

static const char *out_bas = NULL;
/* name of output text file, to which final LP basis should be written
   in MPS format; NULL means no output */

static const char *in_res = NULL;
/* name of input text file, which contains solution in raw format;
   NULL means no solution is provided */

static const char *out_res = NULL;
/* name of output text file, to which final solution should be written
   in raw format; NULL means no output */

static int dir = 0;
/* optimization direction flag:
   0       - not specified
   LPX_MIN - minimization
   LPX_MAX - maximization */

static int scale = 1;
/* if this flag is set, automatic scaling before solving the problem is
   performed; otherwise scaling is not used */

static int method = 0;
/* which method should be used for solving the problem:
   0 - simplex
   1 - interior point
   2 - branch-and-bound */

static int exact = 0;
/* if this flag is set, use lpx_exact rather than lpx_simplex */

#if 0 /* 12/XI-2006 */
static int primal = 0;
/* if this flag is set, use lpx_primal rather than lpx_simplex */
#endif

static int intopt = 0;
/* if this flag is set, use lpx_intopt rather than lpx_integer */

static const char *out_sol = NULL;
/* name of output text file, to which the final solution should be sent
   in plain text format; NULL means no solution output */

static const char *out_bnds = NULL;
/* name of output text file, to which sensitivity bounds should be sent
   in plain text format; NULL means no sensitivity output */

static int tmlim = -1;
/* solution time limit, in seconds */

static int memlim = -1;
/* available memory limit, in megabytes */

static int check = 0;
/* if this flag is set, only input data checking is required */

static int orig = 0;
/* if this flag is set, try to use original names of rows and columns;
   otherwise use plain names */

static const char *out_mps = NULL;
/* name of output text file, to which the problem should be written in
   fixed MPS format; NULL means no MPS output */

static const char *out_freemps = NULL;
/* name of output text file, to which the problem should be written in
   free MPS format; NULL means no MPS output */

static const char *out_cpxlp = NULL;
/* name of output text file, to which the problem should be written in
   CPLEX LP format; NULL means no CPLEX LP output */

static const char *out_txt = NULL;
/* name of output text file, to which the problem should be written in
   plain text format; NULL means no plain text output */

static const char *out_glp = NULL;
/* name of output text file, to which the problem should be written in
   GNU LP format; NULL means no GNU LP output */

static const char *out_pb = NULL;
/* name of output text file, to which the problem should be written in
    OPB format; NULL means no OPB output */

static const char *out_npb = NULL;
/* name of output text file, to which the problem should be written in
    nomralized OPB format; NULL means no normalized OPB output */

static const char *log_file = NULL;
/* name of output text file, to which a hardcopy of all screen output
   should be written; NULL means no hardcopying */

static const char *newname = NULL;
/* new name which has to be assigned to the problem */

static bf_type = 1;
/* which LP basis factorization should be used:
   1 - LP basis factorization option:
   2 - LU + Schur complement + Bartels-Golub update
   3 - LU + Schur complement + Givens rotation update */

static int basis = 1;
/* which initial basis should be used:
   0 - standard initial basis
   1 - advanced initial basis
   2 - Bixby's initial basis */

static int price = 1;
/* which pricing technique should be used:
   0 - textbook pricing
   1 - steepest edge pricing */

static int relax = 1;
/* if this flag is set, the solver uses two-pass ratio test (for both
   primal and dual simplex) proposed by P.Harris; otherwise the solver
   uses the standard "textbook" ratio test */

static int presol = 1;
/* if this flag is set, the solver uses the LP presolver; otherwise the
   LP presolver is not used */

static int xcheck = 0;
/* if this flag is set, the solver checks the final basis using exact
   (bignum) arithmetic */

static int nomip = 0;
/* if this flag is set, the solver considers all integer variables as
   continuous (this allows solving MIP problem as pure LP) */

static int branch = 2;
/* which branching technique should be used:
   0 - branch on first variable
   1 - branch on last variable
   2 - branch using heuristic by Driebeck and Tomlin
   3 - branch on most fractional variable */

static int btrack = 3;
/* which backtracking technique should be used:
   0 - select most recent node (depth first search)
   1 - select earliest node (breadth first search)
   2 - select node using the best projection heuristic
   3 - select node with best local bound */

static double mip_gap = 0.0;
/* relative MIP gap tolerance */

static int binarize = 0;
/* if this flag is set, the solver replaces general integer variables
   by binary ones */

static int use_cuts = 0;
/* if this flag is set, the solver tries to generate cutting planes */

/*----------------------------------------------------------------------
-- display_help - display help.
--
-- This routine displays help information about the program as required
-- by the GNU Coding Standards. */

static void display_help(const char *my_name)
{     print("Usage: %s [options...] filename", my_name);
      print("");
      print("General options:");
#if 0
      print("   --glp             read LP/MIP model in GNU LP format");
#endif
      print("   --mps             read LP/MIP problem in Fixed MPS form"
         "at");
      print("   --freemps         read LP/MIP problem in Free MPS forma"
         "t (default)");
      print("   --cpxlp           read LP/MIP problem in CPLEX LP forma"
         "t");
      print("   --math            read LP/MIP model written in GNU Math"
         "Prog modeling");
      print("                     language");
      print("   -m filename, --model filename");
      print("                     read model section and optional data "
         "section from");
      print("                     filename (the same as --math)");
      print("   -d filename, --data filename");
      print("                     read data section from filename (for "
         "--math only);");
      print("                     if model file also has data section, "
         "that section");
      print("                     is ignored");
      print("   -y filename, --display filename");
      print("                     send display output to filename (for "
         "--math only);");
      print("                     by default the output is sent to term"
         "inal");
      print("   -r filename, --read filename");
      print("                     read solution from filename rather to"
         " find it with");
      print("                     the solver");
      print("   --min             minimization");
      print("   --max             maximization");
      print("   --scale           scale problem (default)");
      print("   --noscale         do not scale problem");
      print("   --simplex         use simplex method (default)");
      print("   --interior        use interior point method (for pure L"
         "P only)");
      print("   -o filename, --output filename");
      print("                     write solution to filename in printab"
         "le format");
      print("   -w filename, --write filename");
      print("                     write solution to filename in plain t"
         "ext format");
      print("   --bounds filename");
      print("                     write sensitivity bounds to filename "
         "in printable");
      print("                     format (LP only)");
      print("   --tmlim nnn       limit solution time to nnn seconds");
      print("   --memlim nnn      limit available memory to nnn megabyt"
         "es");
      print("   --check           do not solve problem, check input dat"
         "a only");
      print("   --name probname   change problem name to probname");
      print("   --plain           use plain names of rows and columns ("
         "default)");
      print("   --orig            try using original names of rows and "
         "columns");
      print("                     (default for --mps)");
#if 0
      print("   --wglp filename   write problem to filename in GNU LP f"
         "ormat");
#endif
      print("   --wmps filename   write problem to filename in Fixed MP"
         "S format");
      print("   --wfreemps filename");
      print("                     write problem to filename in Free MPS"
         " format");
      print("   --wcpxlp filename write problem to filename in CPLEX LP"
         " format");
      print("   --wtxt filename   write problem to filename in printabl"
         "e format");
      print("   --wpb filename    write problem to filename in OPB form"
         "at");
      print("   --wnpb filename   write problem to filename in normaliz"
         "ed OPB format");
      print("   --log filename    write copy of terminal output to file"
         "name");
      print("   -h, --help        display this help information and exi"
         "t");
      print("   -v, --version     display program version and exit");
      print("");
      print("LP basis factorization option:");
      print("   --luf             LU + Forrest-Tomlin update");
      print("                     (faster, less stable; default)");
      print("   --cbg             LU + Schur complement + Bartels-Golub"
         " update");
      print("                     (slower, more stable)");
      print("   --cgr             LU + Schur complement + Givens rotati"
         "on update");
      print("                     (slower, more stable)");
      print("");
      print("Options specific to simplex method:");
      print("   --std             use standard initial basis of all sla"
         "cks");
      print("   --adv             use advanced initial basis (default)")
         ;
      print("   --bib             use Bixby's initial basis");
      print("   --bas filename    read initial basis from filename in M"
         "PS format");
      print("   --steep           use steepest edge technique (default)"
         );
      print("   --nosteep         use standard \"textbook\" pricing");
      print("   --relax           use Harris' two-pass ratio test (defa"
         "ult)");
      print("   --norelax         use standard \"textbook\" ratio test")
         ;
      print("   --presol          use presolver (default; assumes --sca"
         "le and --adv)");
      print("   --nopresol        do not use presolver");
      print("   --exact           use simplex method based on exact ari"
         "thmetic");
      print("   --xcheck          check final basis using exact arithme"
         "tic");
      print("   --wbas filename   write final basis to filename in MPS "
         "format");
      print("");
      print("Options specific to MIP:");
      print("   --nomip           consider all integer variables as con"
         "tinuous");
      print("                     (allows solving MIP as pure LP)");
      print("   --first           branch on first integer variable");
      print("   --last            branch on last integer variable");
      print("   --drtom           branch using heuristic by Driebeck an"
         "d Tomlin");
      print("                     (default)");
      print("   --mostf           branch on most fractional varaible");
      print("   --dfs             backtrack using depth first search");
      print("   --bfs             backtrack using breadth first search")
         ;
      print("   --bestp           backtrack using the best projection h"
         "euristic");
      print("   --bestb           backtrack using node with best local "
         "bound");
      print("                     (default)");
      print("   --mipgap tol      set relative gap tolerance to tol");
      print("   --intopt          use advanced MIP solver");
      print("   --binarize        replace general integer variables by "
         "binary ones");
      print("                     (assumes --intopt)");
      print("   --cover           generate mixed cover cuts");
      print("   --clique          generate clique cuts");
      print("   --gomory          generate Gomory's mixed integer cuts")
         ;
      print("   --mir             generate MIR (mixed integer rounding)"
         " cuts");
      print("   --cuts            generate all cuts above (assumes --in"
         "topt)");
#ifdef _GLP_USE_MIPOPT
      print("   --mipopt          use external MIP solver");
#endif
      print("");
      print("For description of the MPS and CPLEX LP formats see Refere"
         "nce Manual.");
      print("For description of the modeling language see \"GLPK: Model"
         "ing Language");
      print("GNU MathProg\". Both documents are included in the GLPK di"
         "stribution.");
      print("");
      print("See GLPK web page at <http://www.gnu.org/software/glpk/glp"
         "k.html>.");
      print("");
      print("Please report bugs to <bug-glpk@gnu.org>.");
      return;
}

/*----------------------------------------------------------------------
-- display_version - display version.
--
-- This routine displays version information for the program as required
-- by the GNU Coding Standards. */

static void display_version(void)
{     print("GLPSOL: GLPK LP/MIP Solver, Version %s", lib_version());
      print("");
      print("Copyright (C) 2000, 01, 02, 03, 04, 05, 06, 07 Andrew Makh"
         "orin,");
      print("Department for Applied Informatics, Moscow Aviation Instit"
         "ute,");
      print("Moscow, Russia. All rights reserved. E-mail: <mao@mai2.rcn"
         "et.ru>.");
      print("");
      print("This program is free software; you may redistribute it und"
         "er the terms of");
      print("the GNU General Public License. This program has absolutel"
         "y no warranty.");
      return;
}

/*----------------------------------------------------------------------
-- parse_cmdline - parse command-line parameters.
--
-- This routine parses parameters specified in the command line. */

#define p(str) (strcmp(argv[k], str) == 0)

static int parse_cmdline(int argc, const char *argv[])
{     int k;
      for (k = 1; k < argc; k++)
      {  if (p("--mps"))
            format = 0;
         else if (p("--cpxlp") || p("--lpt"))
            format = 1;
         else if (p("--math") || p("-m") || p("--model"))
            format = 2;
         else if (p("--glp"))
            format = 3;
         else if (p("--freemps"))
            format = 4;
         else if (p("-d") || p("--data"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No data input file specifed");
               return 1;
            }
            if (in_data != NULL)
            {  print("Only one data input file allowed");
               return 1;
            }
            in_data = argv[k];
         }
         else if (p("-y") || p("--display"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No display output file specifed");
               return 1;
            }
            if (display != NULL)
            {  print("Only one display output file allowed");
               return 1;
            }
            display = argv[k];
         }
         else if (p("-r") || p("--read"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No input solution file specifed");
               return 1;
            }
            if (display != NULL)
            {  print("Only one input solution file allowed");
               return 1;
            }
            in_res = argv[k];
         }
         else if (p("--min"))
            dir = LPX_MIN;
         else if (p("--max"))
            dir = LPX_MAX;
         else if (p("--scale"))
            scale = 1;
         else if (p("--noscale"))
            scale = 0;
         else if (p("--simplex"))
            method = 0;
         else if (p("--interior"))
            method = 1;
         else if (p("-o") || p("--output"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No solution output file specified");
               return 1;
            }
            if (out_sol != NULL)
            {  print("Only one solution output file allowed");
               return 1;
            }
            out_sol = argv[k];
         }
         else if (p("-w") || p("--write"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No output solution file specifed");
               return 1;
            }
            if (display != NULL)
            {  print("Only one output solution file allowed");
               return 1;
            }
            out_res = argv[k];
         }
         else if (p("--bounds"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No sensitivity bounds output file specified");
               return 1;
            }
            if (out_bnds != NULL)
            {  print("Only one sensitivity bounds output file allowed");
               return 1;
            }
            out_bnds = argv[k];
         }
         else if (p("--tmlim"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No time limit specified");
               return 1;
            }
            if (str2int(argv[k], &tmlim) || tmlim < 0)
            {  print("Invalid time limit `%s'", argv[k]);
               return 1;
            }
         }
         else if (p("--memlim"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No memory limit specified");
               return 1;
            }
            if (str2int(argv[k], &memlim) || memlim < 0)
            {  print("Invalid memory limit `%s'", argv[k]);
               return 1;
            }
         }
         else if (p("--check"))
            check = 1;
         else if (p("--name"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem name specified");
               return 1;
            }
            if (newname != NULL)
            {  print("Only one problem name allowed");
               return 1;
            }
            newname = argv[k];
         }
         else if (p("--plain"))
            orig = 0;
         else if (p("--orig"))
            orig = 1;
         else if (p("--wmps"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No fixed MPS output file specified");
               return 1;
            }
            if (out_mps != NULL)
            {  print("Only one fixed MPS output file allowed");
               return 1;
            }
            out_mps = argv[k];
         }
         else if (p("--wfreemps"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No free MPS output file specified");
               return 1;
            }
            if (out_mps != NULL)
            {  print("Only one free MPS output file allowed");
               return 1;
            }
            out_freemps = argv[k];
         }
         else if (p("--wcpxlp") || p("--wlpt"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No CPLEX LP output file specified");
               return 1;
            }
            if (out_cpxlp != NULL)
            {  print("Only one CPLEX LP output file allowed");
               return 1;
            }
            out_cpxlp = argv[k];
         }
         else if (p("--wtxt"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem output file specified");
               return 1;
            }
            if (out_txt != NULL)
            {  print("Only one problem output file allowed");
               return 1;
            }
            out_txt = argv[k];
         }
         else if (p("--wpb"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem output file specified");
               return 1;
            }
            if (out_pb != NULL)
            {  print("Only one OPB output file allowed");
               return 1;
            }
            out_pb = argv[k];
         }
         else if (p("--wnpb"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem output file specified");
               return 1;
            }
            if (out_npb != NULL)
            {  print("Only one normalized OPB output file allowed");
               return 1;
            }
            out_npb = argv[k];
         }
         else if (p("--log"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No log file specified");
               return 1;
            }
            if (log_file != NULL)
            {  print("Only one log file allowed");
               return 1;
            }
            log_file = argv[k];
         }
         else if (p("--wglp"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem output file specified");
               return 1;
            }
            if (out_glp != NULL)
            {  print("Only one problem output file allowed");
               return 1;
            }
            out_glp = argv[k];
         }
         else if (p("-h") || p("--help"))
         {  display_help(argv[0]);
            return -1;
         }
         else if (p("-v") || p("--version"))
         {  display_version();
            return -1;
         }
         else if (p("--luf"))
            bf_type = 1;
         else if (p("--cbg"))
            bf_type = 2;
         else if (p("--cgr"))
            bf_type = 3;
#if 0 /* 12/XI-2006 */
         else if (p("--primal"))
            primal = 1;
#endif
         else if (p("--std"))
            basis = 0, presol = 0;
         else if (p("--adv"))
            basis = 1, presol = 0;
         else if (p("--bib"))
            basis = 2, presol = 0;
         else if (p("--bas"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No basis input file specifed");
               return 1;
            }
            if (in_bas != NULL)
            {  print("Only one basis input file allowed");
               return 1;
            }
            in_bas = argv[k];
         }
         else if (p("--steep"))
            price = 1;
         else if (p("--nosteep"))
            price = 0;
         else if (p("--relax"))
            relax = 1;
         else if (p("--norelax"))
            relax = 0;
         else if (p("--presol"))
            presol = 1;
         else if (p("--nopresol"))
            presol = 0;
         else if (p("--exact"))
            exact = 1, presol = 0;
         else if (p("--xcheck"))
            xcheck = 1;
         else if (p("--wbas"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No basis output file specified");
               return 1;
            }
            if (out_bas != NULL)
            {  print("Only one basis output file allowed");
               return 1;
            }
            out_bas = argv[k];
         }
         else if (p("--nomip"))
            nomip = 1;
         else if (p("--first"))
            branch = 0;
         else if (p("--last"))
            branch = 1;
         else if (p("--drtom"))
            branch = 2;
         else if (p("--mostf"))
            branch = 3;
         else if (p("--dfs"))
            btrack = 0;
         else if (p("--bfs"))
            btrack = 1;
         else if (p("--bestp"))
            btrack = 2;
         else if (p("--bestb"))
            btrack = 3;
         else if (p("--mipgap"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No relative gap tolerance specified");
               return 1;
            }
            if (str2num(argv[k], &mip_gap) || mip_gap < 0.0)
            {  print("Invalid relative gap tolerance `%s'", argv[k]);
               return 1;
            }
         }
         else if (p("--intopt"))
            intopt = 1;
         else if (p("--binarize"))
            intopt = 1, binarize = 1;
         else if (p("--cover"))
            intopt = 1, use_cuts |= LPX_C_COVER;
         else if (p("--clique"))
            intopt = 1, use_cuts |= LPX_C_CLIQUE;
         else if (p("--gomory"))
            intopt = 1, use_cuts |= LPX_C_GOMORY;
         else if (p("--mir"))
            intopt = 1, use_cuts |= LPX_C_MIR;
         else if (p("--cuts"))
            intopt = 1, use_cuts |= LPX_C_ALL;
#ifdef _GLP_USE_MIPOPT
         else if (p("--mipopt"))
            intopt = 2;
#endif
         else if (argv[k][0] == '-' ||
                 (argv[k][0] == '-' && argv[k][1] == '-'))
         {  print("Invalid option `%s'; try %s --help",
               argv[k], argv[0]);
            return 1;
         }
         else
         {  if (in_file != NULL)
            {  print("Only one input file allowed");
               return 1;
            }
            in_file = argv[k];
         }
      }
      return 0;
}

/*----------------------------------------------------------------------
-- lpx_main - stand-alone LP/MIP solver.
--
-- This main program is called by the control program and manages the
-- solution process. */

int lpx_main(int argc, const char *argv[])
{     LPX *lp = NULL;
      MPL *mpl = NULL;
      int ret;
      glp_ulong start;
      /* parse command line parameters */
      ret = parse_cmdline(argc, argv);
      if (ret < 0)
      {  ret = EXIT_SUCCESS;
         goto done;
      }
      if (ret > 0)
      {  ret = EXIT_FAILURE;
         goto done;
      }
      /* set available memory limit */
      if (memlim >= 0) glp_mem_limit(memlim);
      /* remove all output files specified in the command line */
      if (display != NULL) remove(display);
      if (out_bas != NULL) remove(out_bas);
      if (out_sol != NULL) remove(out_sol);
      if (out_res != NULL) remove(out_res);
      if (out_bnds != NULL) remove(out_bnds);
      if (out_mps != NULL) remove(out_mps);
      if (out_freemps != NULL) remove(out_freemps);
      if (out_cpxlp != NULL) remove(out_cpxlp);
      if (out_txt != NULL) remove(out_txt);
      if (out_glp != NULL) remove(out_glp);
      if (log_file != NULL) remove(log_file);
      /* open hardcopy file, if necessary */
      if (log_file != NULL)
      {  if (lib_open_log(log_file))
         {  print("Unable to create log file");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* read problem data from the input file */
      if (in_file == NULL)
      {  print("No input file specified; try %s --help", argv[0]);
         ret = EXIT_FAILURE;
         goto done;
      }
      switch (format)
      {  case 0:
            lp = lpx_read_mps(in_file);
            if (lp == NULL)
            {  print("MPS file processing error");
               ret = EXIT_FAILURE;
               goto done;
            }
            orig = 1;
            break;
         case 1:
            lp = lpx_read_cpxlp(in_file);
            if (lp == NULL)
            {  print("CPLEX LP file processing error");
               ret = EXIT_FAILURE;
               goto done;
            }
            break;
         case 2:
            /* initialize the translator database */
            mpl = mpl_initialize();
            /* read model section and optional data section */
            ret = mpl_read_model(mpl, (char *)in_file, in_data != NULL);
            if (ret == 4)
err:        {  print("Model processing error");
               ret = EXIT_FAILURE;
               goto done;
            }
            xassert(ret == 1 || ret == 2);
            /* read data section, if necessary */
            if (in_data != NULL)
            {  xassert(ret == 1);
               ret = mpl_read_data(mpl, (char *)in_data);
               if (ret == 4) goto err;
               xassert(ret == 2);
            }
            /* generate model */
            ret = mpl_generate(mpl, (char *)display);
            if (ret == 4) goto err;
            /* extract problem instance */
            lp = lpx_extract_prob(mpl);
            xassert(lp != NULL);
            break;
         case 3:
            lp = lpx_read_prob((char *)in_file);
            if (lp == NULL)
            {  print("GNU LP file processing error");
               ret = EXIT_FAILURE;
               goto done;
            }
            break;
         case 4:
            lp = lpx_read_freemps(in_file);
            if (lp == NULL)
            {  print("MPS file processing error");
               ret = EXIT_FAILURE;
               goto done;
            }
            break;
         default:
            xassert(format != format);
      }
      /* order rows and columns of the constraint matrix */
      lpx_order_matrix(lp);
      /* change problem name (if required) */
      if (newname != NULL) lpx_set_prob_name(lp, newname);
      /* change optimization direction (if required) */
      if (dir != 0) lpx_set_obj_dir(lp, dir);
      /* write problem in fixed MPS format (if required) */
      if (out_mps != NULL)
      {  lpx_set_int_parm(lp, LPX_K_MPSORIG, orig);
         ret = lpx_write_mps(lp, out_mps);
         if (ret != 0)
         {  print("Unable to write problem in fixed MPS format");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem in free MPS format (if required) */
      if (out_freemps != NULL)
      {  ret = lpx_write_freemps(lp, out_freemps);
         if (ret != 0)
         {  print("Unable to write problem in free MPS format");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem in CPLEX LP format (if required) */
      if (out_cpxlp != NULL)
      {  ret = lpx_write_cpxlp(lp, out_cpxlp);
         if (ret != 0)
         {  print("Unable to write problem in CPLEX LP format");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem in plain text format (if required) */
      if (out_txt != NULL)
      {  lpx_set_int_parm(lp, LPX_K_LPTORIG, orig);
         ret = lpx_print_prob(lp, out_txt);
         if (ret != 0)
         {  print("Unable to write problem in plain text format");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem in GNU LP format (if required) */
      if (out_glp != NULL)
      {  ret = lpx_write_prob(lp, (char *)out_glp);
         if (ret != 0)
         {  print("Unable to write problem in GNU LP format");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem in OPB format (if required) */
      if (out_pb != NULL)
      {  ret = lpx_write_pb(lp, out_pb, 0);
         if (ret != 0)
         {  print("Unable to write problem in OPB format");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem in normalized OPB format (if required) */
      if (out_npb != NULL)
      {  ret = lpx_write_pb(lp, out_npb, 1);
         if (ret != 0)
         {  print("Unable to write problem in normalized OPB format");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* if only data check is required, skip computations */
      if (check) goto done;
      /* if solution is provided, read it and skip computations */
      if (in_res != NULL)
      {  switch (method)
         {  case 0:
               if (nomip || lpx_get_class(lp) == LPX_LP)
                  ret = glp_read_sol(lp, in_res);
               else
               {  method = 2;
                  ret = glp_read_mip(lp, in_res);
               }
               break;
            case 1:
               ret = glp_read_ipt(lp, in_res);
               break;
            default:
               xassert(method != method);
         }
         if (ret != 0)
         {  print("Unable to read problem solution");
            ret = EXIT_FAILURE;
            goto done;
         }
         goto ssss;
      }
      /* scale the problem data (if required) */
      if (scale && (!presol || method == 1)) lpx_scale_prob(lp);
      /* build initial LP basis */
      if (method == 0 && !presol && in_bas == NULL)
      {  switch (basis)
         {  case 0:
               lpx_std_basis(lp);
               break;
            case 1:
               if (lpx_get_num_rows(lp) > 0 && lpx_get_num_cols(lp) > 0)
                  lpx_adv_basis(lp);
               break;
            case 2:
               if (lpx_get_num_rows(lp) > 0 && lpx_get_num_cols(lp) > 0)
                  lpx_cpx_basis(lp);
               break;
            default:
               xassert(basis != basis);
         }
      }
      /* or read initial basis from input text file in MPS format */
      if (in_bas != NULL)
      {  if (method != 0)
         {  print("Initial LP basis is useless for interior-point solve"
               "r and therefore ignored");
            goto nobs;
         }
         lpx_set_int_parm(lp, LPX_K_MPSORIG, orig);
         ret = lpx_read_bas(lp, in_bas);
         if (ret != 0)
         {  print("Unable to read initial LP basis");
            ret = EXIT_FAILURE;
            goto done;
         }
         if (presol)
         {  presol = 0;
            print("LP presolver disabled because initial LP basis has b"
               "een provided");
         }
nobs:    ;
      }
      /* set some control parameters, which might be changed in the
         command line */
      lpx_set_int_parm(lp, LPX_K_BFTYPE, bf_type);
      lpx_set_int_parm(lp, LPX_K_PRICE, price);
      if (!relax) lpx_set_real_parm(lp, LPX_K_RELAX, 0.0);
      lpx_set_int_parm(lp, LPX_K_PRESOL, presol);
      lpx_set_int_parm(lp, LPX_K_BRANCH, branch);
      lpx_set_int_parm(lp, LPX_K_BTRACK, btrack);
      lpx_set_real_parm(lp, LPX_K_TMLIM, (double)tmlim);
      lpx_set_int_parm(lp, LPX_K_BINARIZE, binarize);
      lpx_set_int_parm(lp, LPX_K_USECUTS, use_cuts);
      lpx_set_real_parm(lp, LPX_K_MIPGAP, mip_gap);
      /* solve the problem */
      start = xtime();
      switch (method)
      {  case 0:
            if (nomip || lpx_get_class(lp) == LPX_LP)
            {  ret = (!exact ? lpx_simplex(lp) : lpx_exact(lp));
               if (xcheck)
               {  if (!presol || ret == LPX_E_OK)
                     lpx_exact(lp);
                  else
                     print("If you need checking final basis for non-op"
                        "timal solution, use --nopresol");
               }
               if (presol && ret != LPX_E_OK && (out_bas != NULL ||
                  out_sol != NULL))
                  print("If you need actual output for non-optimal solu"
                     "tion, use --nopresol");
            }
            else
            {  method = 2;
               if (!intopt)
               {  ret = (!exact ? lpx_simplex(lp) : lpx_exact(lp));
                  if (xcheck && (!presol || ret == LPX_E_OK))
                     lpx_exact(lp);
                  lpx_integer(lp);
               }
               else if (intopt == 1)
                  lpx_intopt(lp);
#ifdef _GLP_USE_MIPOPT
               else
                  glp_mipopt(lp);
#endif
            }
            break;
         case 1:
            if (nomip || lpx_get_class(lp) == LPX_LP)
               lpx_interior(lp);
            else
            {  print("Interior-point method is not able to solve MIP pr"
                  "oblem; use --simplex");
               ret = EXIT_FAILURE;
               goto done;
            }
            break;
         default:
            xassert(method != method);
      }
      /* display statistics */
      print("Time used:   %.1f secs", xdifftime(xtime(), start));
      {  glp_ulong tpeak;
         char buf[50];
         lib_mem_usage(NULL, NULL, NULL, &tpeak);
         print("Memory used: %.1f Mb (%s bytes)",
            (4294967296.0 * tpeak.hi + tpeak.lo) / 1048576.0,
            ultoa(tpeak, buf, 10));
      }
ssss: if (mpl != NULL && mpl_has_solve_stmt(mpl))
      {  int n, j, round;
         /* store the solution to the translator database */
         n = lpx_get_num_cols(lp);
         round = lpx_get_int_parm(lp, LPX_K_ROUND);
         lpx_set_int_parm(lp, LPX_K_ROUND, 1);
         switch (method)
         {  case 0:
               for (j = 1; j <= n; j++)
                  mpl_put_col_value(mpl, j, lpx_get_col_prim(lp, j));
               break;
            case 1:
               for (j = 1; j <= n; j++)
                  mpl_put_col_value(mpl, j, lpx_ipt_col_prim(lp, j));
               break;
            case 2:
               for (j = 1; j <= n; j++)
                  mpl_put_col_value(mpl, j, lpx_mip_col_val(lp, j));
               break;
            default:
               xassert(method != method);
         }
         lpx_set_int_parm(lp, LPX_K_ROUND, round);
         /* perform postsolving */
         ret = mpl_postsolve(mpl);
         if (ret == 4)
         {  print("Model postsolving error");
            ret = EXIT_FAILURE;
            goto done;
         }
         xassert(ret == 3);
      }
      /* write final LP basis (if required) */
      if (out_bas != NULL)
      {  lpx_set_int_parm(lp, LPX_K_MPSORIG, orig);
         ret = lpx_write_bas(lp, out_bas);
         if (ret != 0)
         {  print("Unable to write final LP basis");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem solution found by the solver (if required) */
      if (out_sol != NULL)
      {  switch (method)
         {  case 0:
               ret = lpx_print_sol(lp, out_sol);
               break;
            case 1:
               ret = lpx_print_ips(lp, out_sol);
               break;
            case 2:
               ret = lpx_print_mip(lp, out_sol);
               break;
            default:
               xassert(method != method);
         }
         if (ret != 0)
         {  print("Unable to write problem solution");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      if (out_res != NULL)
      {  switch (method)
         {  case 0:
               ret = glp_write_sol(lp, out_res);
               break;
            case 1:
               ret = glp_write_ipt(lp, out_res);
               break;
            case 2:
               ret = glp_write_mip(lp, out_res);
               break;
            default:
               xassert(method != method);
         }
         if (ret != 0)
         {  print("Unable to write problem solution");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write sensitivity bounds information (if required) */
      if (out_bnds != NULL)
      {  if (method != 0)
         {  print("Cannot write sensitivity bounds information for inte"
               "rior-point or MIP solution");
            ret = EXIT_FAILURE;
            goto done;
         }
         ret = lpx_print_sens_bnds(lp, out_bnds);
         if (ret != 0)
         {  print("Unable to write sensitivity bounds information");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* all seems to be ok */
      ret = EXIT_SUCCESS;
done: /* delete the problem object */
      if (lp != NULL) lpx_delete_prob(lp);
      /* if the translator database exists, destroy it */
      if (mpl != NULL) mpl_terminate(mpl);
      xassert(gmp_pool_count() == 0);
      gmp_free_mem();
      /* close the hardcopy file */
      if (log_file != NULL) lib_close_log();
      /* check that no memory blocks are still allocated */
      {  int count;
         glp_ulong total;
         lib_mem_usage(&count, NULL, &total, NULL);
         xassert(count == 0);
         xassert(total.lo == 0 && total.hi == 0);
      }
      /* free the library environment */
      lib_free_env();
      /* return to the control program */
      return ret;
}

/* eof */
