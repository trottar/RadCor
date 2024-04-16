#include "F1F2IN21.h"

// a C++ wrapper for the fortran subroutines in strufunc_f1f221_dis6.f
extern "C"
{
    // F1F2IN21
    void F1F2IN21_(double Z, double A, double QSQ, double WSQ, double *F1, double *F2);
}

void F1F2IN21(double Z, double A, double QSQ, double WSQ, double F1, double F2) {
    F1F2IN21_(Z, A, QSQ, WSQ, &F1, &F2);
}
