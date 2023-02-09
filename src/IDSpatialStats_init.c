#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void get_pi_a_clustsurvey_hh_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_a_clustsurvey_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_a_typed_clustsurvey_hh_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_a_typed_clustsurvey_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_b_clustsurvey_hh_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_b_clustsurvey_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_b_typed_clustsurvey_hh_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_b_typed_clustsurvey_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_clustsurvey(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_clustsurvey_hh_wts_window(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_clustsurvey_window(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_clustsurvey_wts_window(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_hh_typed(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_typed(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_typed_survey(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_pi_typed_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_clustsurvey(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_clustsurvey_hh_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_clustsurvey_hh_wts_window(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_clustsurvey_window(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_clustsurvey_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_clustsurvey_wts_window(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_hh_typed(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_typed(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_typed_clustsurvey_hh_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_typed_clustsurvey_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_typed_survey(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_tau_typed_wts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_theta_typed(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_theta_typed_survey(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP get_pi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_tau(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_theta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"get_pi_a_clustsurvey_hh_wts",       (DL_FUNC) &get_pi_a_clustsurvey_hh_wts,       16},
    {"get_pi_a_clustsurvey_wts",          (DL_FUNC) &get_pi_a_clustsurvey_wts,          14},
    {"get_pi_a_typed_clustsurvey_hh_wts", (DL_FUNC) &get_pi_a_typed_clustsurvey_hh_wts, 20},
    {"get_pi_a_typed_clustsurvey_wts",    (DL_FUNC) &get_pi_a_typed_clustsurvey_wts,    18},
    {"get_pi_b_clustsurvey_hh_wts",       (DL_FUNC) &get_pi_b_clustsurvey_hh_wts,       16},
    {"get_pi_b_clustsurvey_wts",          (DL_FUNC) &get_pi_b_clustsurvey_wts,          14},
    {"get_pi_b_typed_clustsurvey_hh_wts", (DL_FUNC) &get_pi_b_typed_clustsurvey_hh_wts, 20},
    {"get_pi_b_typed_clustsurvey_wts",    (DL_FUNC) &get_pi_b_typed_clustsurvey_wts,    18},
    {"get_pi_clustsurvey",                (DL_FUNC) &get_pi_clustsurvey,                10},
    {"get_pi_clustsurvey_hh_wts_window",  (DL_FUNC) &get_pi_clustsurvey_hh_wts_window,  16},
    {"get_pi_clustsurvey_window",         (DL_FUNC) &get_pi_clustsurvey_window,         11},
    {"get_pi_clustsurvey_wts_window",     (DL_FUNC) &get_pi_clustsurvey_wts_window,     14},
    {"get_pi_hh_typed",                   (DL_FUNC) &get_pi_hh_typed,                   12},
    {"get_pi_typed",                      (DL_FUNC) &get_pi_typed,                      11},
    {"get_pi_typed_survey",               (DL_FUNC) &get_pi_typed_survey,               11},
    {"get_pi_typed_wts",                  (DL_FUNC) &get_pi_typed_wts,                  12},
    {"get_tau_clustsurvey",               (DL_FUNC) &get_tau_clustsurvey,               10},
    {"get_tau_clustsurvey_hh_wts",        (DL_FUNC) &get_tau_clustsurvey_hh_wts,        18},
    {"get_tau_clustsurvey_hh_wts_window", (DL_FUNC) &get_tau_clustsurvey_hh_wts_window, 16},
    {"get_tau_clustsurvey_window",        (DL_FUNC) &get_tau_clustsurvey_window,        11},
    {"get_tau_clustsurvey_wts",           (DL_FUNC) &get_tau_clustsurvey_wts,           16},
    {"get_tau_clustsurvey_wts_window",    (DL_FUNC) &get_tau_clustsurvey_wts_window,    14},
    {"get_tau_hh_typed",                  (DL_FUNC) &get_tau_hh_typed,                  13},
    {"get_tau_typed",                     (DL_FUNC) &get_tau_typed,                     12},
    {"get_tau_typed_clustsurvey_hh_wts",  (DL_FUNC) &get_tau_typed_clustsurvey_hh_wts,  22},
    {"get_tau_typed_clustsurvey_wts",     (DL_FUNC) &get_tau_typed_clustsurvey_wts,     20},
    {"get_tau_typed_survey",              (DL_FUNC) &get_tau_typed_survey,              12},
    {"get_tau_typed_wts",                 (DL_FUNC) &get_tau_typed_wts,                 12},
    {"get_theta_typed",                   (DL_FUNC) &get_theta_typed,                   11},
    {"get_theta_typed_survey",            (DL_FUNC) &get_theta_typed_survey,            11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"get_pi",    (DL_FUNC) &get_pi,    7},
    {"get_tau",   (DL_FUNC) &get_tau,   8},
    {"get_theta", (DL_FUNC) &get_theta, 7},
    {NULL, NULL, 0}
};

void R_init_IDSpatialStats(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}