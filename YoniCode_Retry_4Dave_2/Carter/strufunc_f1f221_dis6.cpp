#include "strufunc_f1f221_dis6.h"

// a C++ wrapper for the fortran subroutines in strufunc_f1f221_dis6.f
extern "C"
{
    // f1sfun
    double f1sfun_(double *x, double *q2);

    // f2sfun
    double f2sfun_(double *x, double *q2);

    // g1sfun
    double g1sfun_(double *x, double *q2);

    // g2sfun
    double g2sfun_(double *x, double *q2);
}

double f1sfun(double *x, double *q2) {
    return f1sfun_(x, q2);
}

double f2sfun(double *x, double *q2) {
    return f2sfun_(x, q2);
}

double g1sfun(double *x, double *q2) {
    return g1sfun_(x, q2);
}

double g2sfun(double *x, double *q2) {
    return g2sfun_(x, q2);
}
