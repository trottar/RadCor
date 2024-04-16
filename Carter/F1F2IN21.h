#ifndef F1F2IN21_CWRAPPER_H
#define F1F2IN21_CWRAPPER_H

// a C++ wrapper for the fortran subroutines in strufunc_f1f221_dis6.f
extern "C"
{
    // F1F2IN21
    void F1F2IN21(double Z, double A, double QSQ, double WSQ, double *F1, double *F2);
}

#endif
