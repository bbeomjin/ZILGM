#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
 *   The following name(s) appear with different usages
 *     e.g., with different numbers of arguments:
 *
 *         penGLM
 *
 *           This needs to be resolved in the tables and any declarations.
 *           */

/* FIXME: 
 *    Check these declarations against the C/Fortran source code.
 *    */

/* .Fortran calls */
extern void F77_NAME(checkconvergence)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(deveval)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(glmlink)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(loglikfor)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(outloop)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(penglm)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pred)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zeval)(void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
	    {"checkconvergence", (DL_FUNC) &F77_NAME(checkconvergence),  7},
	    {"deveval",          (DL_FUNC) &F77_NAME(deveval),           7},
	    {"glmlink",          (DL_FUNC) &F77_NAME(glmlink),           6},
	    {"loglikfor",        (DL_FUNC) &F77_NAME(loglikfor),         7},
	    {"outloop",          (DL_FUNC) &F77_NAME(outloop),          34},
	    {"penglm",           (DL_FUNC) &F77_NAME(penglm),            7},
	    {"pred",             (DL_FUNC) &F77_NAME(pred),              9},
	    {"zeval",            (DL_FUNC) &F77_NAME(zeval),             7},
					        {NULL, NULL, 0}
};

void R_init_mpath(DllInfo *dll)
{
	    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
	        R_useDynamicSymbols(dll, FALSE);
}