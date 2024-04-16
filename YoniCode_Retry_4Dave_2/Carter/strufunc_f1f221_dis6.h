#ifndef strufunc_f1f221_dis6_CWRAPPER_H
#define strufunc_f1f221_dis6_CWRAPPER_H

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

#endif
