#include "strufunc_f1f221_dis6.h"
#include "F1F2IN21.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

const double deg2rad = 0.0174533;
const double Mp = 0.938;

int main() {
// void xs_gen_dis6(double Ebeam = 10.38 /*GeV*/, double theta = 30 /*deg*/) {
    // double Ep = 1.0; // GeV
    // double dEp = 0.1;
    // double theta_rad = theta*deg2rad;
    // double Q2 = 2*Ebeam*Ep*
    double F1 = 0, F2 = 0;
    double Q2, x;
    Q2 = 1.016;
    x = 0.058;
    cout << "Q^2 = " << Q2 << "GeV^2 \t" << "x = " << x << endl;
    cout << "F1: " << f1sfun_(&x, &Q2) << endl;
    cout << "F2: " << f2sfun_(&x, &Q2) << endl;
    // F1F2IN21(1.0, 1.0, 5.0, 2.0, F1, F2);
    // cout << "F1F2IN21: " << F1 << "\t" << F2 << endl;
}
