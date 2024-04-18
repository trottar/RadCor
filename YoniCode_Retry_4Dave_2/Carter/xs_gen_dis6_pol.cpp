#include "strufunc_f1f221_dis6.h"
#include "F1F2IN21.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double deg2rad = 0.0174533;   // degrees to radians conversion
const double Mp = 0.938;            // proton mass in [GeV]
const double M3he = 2.8094141619;     // 3He mass in [GeV]
const double alpha = 0.0072973525693;   // fine structure constant
const double A3he = 3.0;
const double nbarn = 0.389e6;           // barn: (1 GeV)**-2 = 0.389e-3 barn


double sqr(double x) {
    return pow(x, 2.0);
}

double q2_calc(double E, double Ep, double theta) {
    return 2*E*Ep*(1 - cos(theta*deg2rad));
}

double x_calc(double E, double Ep, double theta) {
    double Q2 = q2_calc(E, Ep, theta);
    double nu = E - Ep;
    double denom = 2*Mp*nu;
    return Q2/denom;
}

double w_calc(double E, double Ep, double theta) {
    double Q2 = 0, nu = 0, W2 = 0;
    Q2 = q2_calc(E, Ep, theta);
    nu = E - Ep;
    W2 = sqr(Mp) + 2*Mp*nu - Q2;
    return sqrt(W2);
}

// Calculates Mott XS in units [GeV^-2]
double xs_mott(double E /*[GeV]*/, double theta /*degrees*/) {
    double num = 0, denom = 0;
    num = sqr(alpha)*sqr(cos(theta*deg2rad/2.0));
    denom = 4*sqr(E)*sqr(sin(theta*deg2rad/2.0));
    return num/denom;
}

// Calculates unpolarized XS from form factors in units []
// double xs_calc(double E, double Ep, double theta, double Q2, double x) {
//     double nu = 0, first_term = 0, second_term = 0, f1 = 0, f2 = 0, xsMott = 0;
//     nu = E - Ep;
//     f1 = A3he*f1sfun_(&x, &Q2);
//     f2 = A3he*f2sfun_(&x, &Q2);
//     first_term = f2/nu;
//     second_term = 2*f1*sqr(tan(theta*deg2rad/2.0))/Mp;
//     xsMott = xs_mott(E, theta)*nbarn;   // convert to units of nbarn
//     return xsMott*(first_term + second_term);
// }

// Calculates unpolarized XS from form factors in units []
double xs_calc(double E, double Ep, double theta, double Q2, double x) {
    double nu = 0, first_term = 0, second_term = 0, f1 = 0, f2 = 0, xsMott = 0;
    nu = E - Ep;
    f1 = A3he*f1sfun_(&x, &Q2);
    f2 = A3he*f2sfun_(&x, &Q2);
    first_term = f2/nu;
    second_term = 2*f1*sqr(tan(theta*deg2rad/2.0))/Mp;
    xsMott = xs_mott(E, theta)*nbarn;   // convert to units of nbarn
    return xsMott*(first_term + second_term);
}

double xsLong_calc(double E, double Ep, double theta, double Q2, double x) {
    double nu = 0, first_term = 0, second_term = 0, outer = 0, g1 = 0, g2 = 0;
    nu = E - Ep;
    g1 = A3he*g1sfun_(&x, &Q2);
    g2 = A3he*g2sfun_(&x, &Q2);
    outer = 4*sqr(alpha)*Ep/(nu*E*Q2)/Mp;
    first_term = E + Ep*cos(theta*deg2rad);
    first_term *= g1;
    second_term = 2.0*Mp*x*g2;
    return outer*(first_term - second_term)*nbarn;
}

double xsTrans_calc(double E, double Ep, double theta, double Q2, double x) {
    double nu = 0, first_term = 0, second_term = 0, outer = 0, g1 = 0, g2 = 0;
    nu = E - Ep;
    g1 = A3he*g1sfun_(&x, &Q2);
    g2 = A3he*g2sfun_(&x, &Q2);
    outer = 4*sqr(alpha*Ep)/(nu*E*Q2)*sin(theta*deg2rad)/Mp;
    first_term = g1;
    second_term = 2.0*E/nu;
    second_term *= g2;
    return outer*(first_term + second_term)*nbarn;
}

int main() {
// void xs_gen_dis6(double Ebeam = 10.38 /*GeV*/, double theta = 30 /*deg*/) {
    double E = 10.38;   // GeV
    int Elabel = int(E*1000);
    double dE = 1.0;
    double Ep = 0.0;    // GeV
    double dEp = 0.001;
    double theta = 30.0;
    int thetalabel = int(theta);
    int num_Ebeams = 5;
    
    double F1 = 0.0, F2 = 0.0;  // These are form factors per nucleus, so need to multiply form factors per nucleon by A
    double Q2 = 0.0, x = 0.0;
    double xsLong = 0.0, xsTrans = 0.0, nu = 0.0;  


    E = 5.73;
    Ep = 1.31;
    theta = 34.9;
    Q2 = 2.708;
    x = 0.327;
    cout << "Q^2 = " << Q2 << "GeV^2 \t" << "x = " << x << endl;
    // cout << "F1: " << f1sfun_(&x, &Q2) << endl;
    // cout << "F2: " << f2sfun_(&x, &Q2) << endl;
    // cout << "G1: " << g1sfun_(&x, &Q2) << endl;
    // cout << "G2: " << g2sfun_(&x, &Q2) << endl;
    cout << "sig_unpol: " << xs_calc(E,Ep,theta,Q2,x) << endl;
    cout << "sig_par: " << xsLong_calc(E,Ep,theta,Q2,x) << endl;
    cout << "sig_perp: " << xsTrans_calc(E,Ep,theta,Q2,x) << endl;
    cout << "A_par: " << xsLong_calc(E,Ep,theta,Q2,x)/(2*xs_calc(E,Ep,theta,Q2,x)) << endl;
    cout << "A_perp: " << xsTrans_calc(E,Ep,theta,Q2,x)/(2*xs_calc(E,Ep,theta,Q2,x)) << endl;

    E = 5.73;
    Ep = 1.71;
    theta = 34.88;
    Q2 = 3.515;
    x = 0.466;
    cout << "Q^2 = " << Q2 << "GeV^2 \t" << "x = " << x << endl;
    // cout << "F1: " << f1sfun_(&x, &Q2) << endl;
    // cout << "F2: " << f2sfun_(&x, &Q2) << endl;
    // cout << "G1: " << g1sfun_(&x, &Q2) << endl;
    // cout << "G2: " << g2sfun_(&x, &Q2) << endl;
    cout << "sig_unpol: " << xs_calc(E,Ep,theta,Q2,x) << endl;
    cout << "sig_par: " << xsLong_calc(E,Ep,theta,Q2,x) << endl;
    cout << "sig_perp: " << xsTrans_calc(E,Ep,theta,Q2,x) << endl;
    cout << "A_par: " << xsLong_calc(E,Ep,theta,Q2,x)/(2*xs_calc(E,Ep,theta,Q2,x)) << endl;
    cout << "A_perp: " << xsTrans_calc(E,Ep,theta,Q2,x)/(2*xs_calc(E,Ep,theta,Q2,x)) << endl;

    E = 5.73;
    Ep = 1.44;
    theta = 44.94;
    Q2 = 4.831;
    x = 0.601;
    cout << "Q^2 = " << Q2 << "GeV^2 \t" << "x = " << x << endl;
    // cout << "F1: " << f1sfun_(&x, &Q2) << endl;
    // cout << "F2: " << f2sfun_(&x, &Q2) << endl;
    // cout << "G1: " << g1sfun_(&x, &Q2) << endl;
    // cout << "G2: " << g2sfun_(&x, &Q2) << endl;
    cout << "sig_unpol: " << xs_calc(E,Ep,theta,Q2,x) << endl;
    cout << "sig_par: " << xsLong_calc(E,Ep,theta,Q2,x) << endl;
    cout << "sig_perp: " << xsTrans_calc(E,Ep,theta,Q2,x) << endl;
    cout << "A_par: " << xsLong_calc(E,Ep,theta,Q2,x)/(2*xs_calc(E,Ep,theta,Q2,x)) << endl;
    cout << "A_perp: " << xsTrans_calc(E,Ep,theta,Q2,x)/(2*xs_calc(E,Ep,theta,Q2,x)) << endl;

    return 0;

    ofstream outputLongFile, outputTransFile;
    string outputLongFileName, outputTransFileName, inputLongFilePath, inputTransFilePath;
    // outputFile.open("E10380_30deg.dat");

    int countL = 0, countT = 0;

    for (int i=0; i < num_Ebeams; i++) {
        countL = 0;
        countT = 0;
        if (i != 0) { 
            E -= dE;
            Elabel -= int(dE*1000);
        }
        outputLongFileName = "Data/Long/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_long.dat";
        outputLongFile.open(outputLongFileName);
        outputTransFileName = "Data/Trans/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_trans.dat";
        outputTransFile.open(outputTransFileName);
        inputLongFilePath = "/home/carterhedinger/Zheng/RadCor/CAnalyzer-master/example/data/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_long.dat";
        inputTransFilePath = "/home/carterhedinger/Zheng/RadCor/CAnalyzer-master/example/data/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_trans.dat";

        Ep = dEp;
        Q2 = q2_calc(E, Ep, theta);
        x = x_calc(E, Ep, theta);
        F1 = A3he*f1sfun_(&x, &Q2);
        F2 = A3he*f2sfun_(&x, &Q2);
        while(Ep < E) {
            // cout << x << "\t" << Q2 << "\t" << Ep << endl;
            if (x < 1.0) {
                xsLong = xsLong_calc(E, Ep, theta, Q2, x);
                xsTrans = xsTrans_calc(E, Ep, theta, Q2, x);
                nu = (E - Ep)*1000;
                if (xsLong > 0) { 
                    outputLongFile << nu << "\t" << xsLong << "\t" << 0.0 << "\t" << 0.0 << endl;
                    countL++;
                }
                if (xsTrans > 0) {
                    outputTransFile << nu << "\t" << xsTrans << "\t" << 0.0 << "\t" << 0.0 << endl;
                    countT++;
                }
            }
            Ep += dEp;
            Q2 = q2_calc(E, Ep, theta);
            x = x_calc(E, Ep, theta);
            F1 = A3he*f1sfun_(&x, &Q2);
            F2 = A3he*f2sfun_(&x, &Q2);
        }
        outputLongFile.close();
        cout << countL << " data points saved in" << outputLongFileName << endl;
        outputTransFile.close();
        cout << countT << " data points saved in" << outputTransFileName << endl;

        ifstream sourceLongFile(outputLongFileName, ios::binary);
        ofstream copyLongFile(inputLongFilePath, ios::binary);
        if (!copyLongFile) {
            cout << "Error opening " << inputLongFilePath << endl;
        }
        copyLongFile << sourceLongFile.rdbuf();
        sourceLongFile.close();
        copyLongFile.close();
        cout << "Data copied successfully to " << inputLongFilePath << endl;


        ifstream sourceTransFile(outputTransFileName, ios::binary);
        ofstream copyTransFile(inputTransFilePath, ios::binary);
        if (!copyTransFile) {
            cout << "Error opening " << inputTransFilePath << endl;
        }
        copyTransFile << sourceTransFile.rdbuf();
        sourceTransFile.close();
        copyTransFile.close();
        cout << "Data copied successfully to " << inputTransFilePath << endl;

    }

}
