#include "strufunc_jam.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <iomanip>
#include <random>

using namespace std;

// const double deg2rad = 0.0174533;   // degrees to radians conversion
// const double Mp = 0.938;            // proton mass in [GeV]
// const double M3he = 2.8094141619;     // 3He mass in [GeV]
// const double alpha = 0.0072973525693;   // fine structure constant
// const double A3he = 3.0;
// const double nbarn = 0.389e6;           // barn: (1 GeV)**-2 = 0.389e-3 barn
// // const double maxPercDiff = 1.5;
// // const double minPercDiff = 4/maxPercDiff;


// double sqr(double x) {
//     return pow(x, 2.0);
// }

// // double perc_diff(double x, double y) {
// //     if (x==0 || y==0) return 1e20;
// //     double z = x/y;
// //     return fabs(1 - z)/(fabs(1 + z)/2);
// // }

// double q2_calc(double E, double Ep, double theta) {
//     return 2*E*Ep*(1 - cos(theta*deg2rad));
// }

// double x_calc(double E, double Ep, double theta) {
//     double Q2 = q2_calc(E, Ep, theta);
//     double nu = E - Ep;
//     double denom = 2*Mp*nu;
//     return Q2/denom;
// }

// double w_calc(double E, double Ep, double theta) {
//     double Q2 = 0, nu = 0, W2 = 0;
//     Q2 = q2_calc(E, Ep, theta);
//     nu = E - Ep;
//     W2 = sqr(Mp) + 2*Mp*nu - Q2;
//     return sqrt(W2);
// }

// double nu_calc(double x, double q2) {
//     double inv = x*2*Mp/q2;
//     return 1/inv;
// }

// // Calculates Mott XS in units [GeV^-2]
// // Latex - \frac{d^2\sigma_{Mott}}{d\Omega} = \frac{\alpha^2 cos^2\frac{\theta}{2}}{4E^2sin^4\frac{\theta}{2}}
// double xs_mott(double E /*[GeV]*/, double theta /*degrees*/) {
//     double num = 0, denom = 0;
//     num = sqr(alpha)*sqr(cos(theta*deg2rad/2.0));
//     denom = 4*sqr(E)*sqr(sqr(sin(theta*deg2rad/2.0)));
//     return num/denom;
// }

// // Calculates unpolarized XS from structure functions in units [nb/MeV]
// // Latex - \frac{d^2\sigma}{d\Omega dE} = \frac{d^2\sigma_{Mott}}{d\Omega}\cdot\left(\frac{F_2(x,Q^2)}{\nu} + \frac{2}{M}tan^2\frac{\theta}{2}F_1(x,Q^2)\right)
// double xs_calc(double E, double Ep, double theta, double F1, double F2) {
//     double nu = 0, first_term = 0, second_term = 0, xsMott = 0;
//     nu = E - Ep;
//     first_term = F2/nu;
//     second_term = 2*F1*sqr(tan(theta*deg2rad/2.0))/Mp;
//     xsMott = xs_mott(E, theta)*nbarn;   // convert to units of nbarn
//     return xsMott*(first_term + second_term);
// }

// // Calculates difference between parallel and antiparallel longitudinal (parallel) XS from structure functions
// // Latex - \frac{d^2\sigma_{\parallel}}{d\Omega dE} = \frac{4\alpha^2E'}{\nu EQ^2M}\cdot\left[ (E + E'cos\theta)g_1(x,Q^2) + 2Mxg_2(x,Q^2) \right]
// double xsLong_calc(double E, double Ep, double theta, double Q2, double x, double g1, double g2) {
//     double nu = 0, first_term = 0, second_term = 0, outer = 0;
//     nu = E - Ep;
//     outer = 4*sqr(alpha)*Ep/(nu*E*Q2)/Mp;
//     first_term = E + Ep*cos(theta*deg2rad);
//     first_term *= g1;
//     second_term = 2.0*Mp*x*g2;
//     return outer*(first_term - second_term)*nbarn;
// }

// // Calculates transverse (perpendicular) XS from structure functions
// // Latex - \frac{d^2\sigma_{\perp}}{d\Omega dE} = \frac{4\alpha^2E'^2}{\nu EQ^2Msin\theta}\cdot\left( g_1(x,Q^2) + \frac{2E}{\nu}g_2(x,Q^2) \right)
// double xsTrans_calc(double E, double Ep, double theta, double Q2, double g1, double g2) {
//     double nu = 0, first_term = 0, second_term = 0, outer = 0;
//     nu = E - Ep;
//     outer = 4*sqr(alpha*Ep)/(nu*E*Q2)*sin(theta*deg2rad)/Mp;
//     first_term = g1;
//     second_term = 2.0*E/nu;
//     second_term *= g2;
//     return outer*(first_term + second_term)*nbarn;
// }

int main() {
// void xs_gen_dis6(double Ebeam = 10.38 /*GeV*/, double theta = 30 /*deg*/) {
    // double E = 10.38;   // GeV
    // int Elabel = int(E*1000);
    // double dE = 1.0;
    // double Ep = 0.0;    // GeV
    // double dEp = 0.001;
    // vector<double> theta_values = {11.0, 18.0, 30.0};
    // double theta = theta_values[0];
    // int thetalabel = int(theta);
    // int num_Ebeams = 5;
    
    // double F1 = 0.0, F2 = 0.0, g1 = 0.0, g2 = 0.0;
    // Point pt;
    // double Q2 = 0.0, x = 0.0, w = 0.0;
    // double xsUnpol = 0.0, xsLong = 0.0, xsTrans = 0.0, nu = 0.0;
    // double xsUnpol_prev = xsUnpol, xsLong_prev = xsLong, xsTrans_prev = xsTrans;

    map<double, vector<Point>> dataPoints;
    vector<double> q2_keys;
    getDataPoints(dataPoints, q2_keys);

    double q2min = q2_keys[0];
    double q2max = q2_keys[q2_keys.size()-1];
    // double xmin = dataPoints[q2min][0].x;
    // double xmax = dataPoints[q2max][dataPoints[q2max].size()-1].x;

    // cout << dataPoints.size() << " " << q2_keys.size() << endl;


    // E = 5.73;
    // Ep = 1.31;
    // theta = 34.9;
    // Q2 = 2.708;
    // x = 0.327;
    // cout << "Q^2 = " << Q2 << "GeV^2 \t" << "x = " << x << endl;
    // pt.setEqual(interpolate2DLinear(x, Q2, dataPoints, q2_keys));
    // cout << "F1: " << pt.F1 << endl;
    // cout << "F2: " << pt.F2 << endl;
    // cout << "G1: " << pt.g1 << endl;
    // cout << "G2: " << pt.g2 << endl;
    // cout << "sig_unpol: " << xs_calc(E,Ep,theta,pt.F1,pt.F2) << endl;
    // cout << "sig_par: " << xsLong_calc(E,Ep,theta,Q2,x,pt.g1,pt.g2) << endl;
    // cout << "sig_perp: " << xsTrans_calc(E,Ep,theta,Q2,pt.g1,pt.g2) << endl;
    // cout << "A_par: " << xsLong_calc(E,Ep,theta,Q2,x,pt.g1,pt.g2)/(2*xs_calc(E,Ep,theta,pt.F1,pt.F2)) << endl;
    // // \mathrm{A}_{\parallel} = \frac{\Delta\frac{d^2\sigma_{\parallel}}{d\Omega dE}}{2 \cdot \frac{d^2\sigma}{d\Omega dE}}
    // cout << "A_perp: " << xsTrans_calc(E,Ep,theta,Q2,pt.g1,pt.g2)/(2*xs_calc(E,Ep,theta,pt.F1,pt.F2)) << endl;
    // // \mathrm{A}_{\perp} = \frac{\Delta\frac{d^2\sigma_{\perp}}{d\Omega dE}}{2 \cdot \frac{d^2\sigma}{d\Omega dE}}

    // E = 5.73;
    // Ep = 1.71;
    // theta = 34.88;
    // Q2 = 3.515;
    // x = 0.466;
    // cout << "Q^2 = " << Q2 << "GeV^2 \t" << "x = " << x << endl;
    // pt.setEqual(interpolate2DLinear(x, Q2, dataPoints, q2_keys));
    // cout << "F1: " << pt.F1 << endl;
    // cout << "F2: " << pt.F2 << endl;
    // cout << "G1: " << pt.g1 << endl;
    // cout << "G2: " << pt.g2 << endl;
    // cout << "sig_unpol: " << xs_calc(E,Ep,theta,pt.F1,pt.F2) << endl;
    // cout << "sig_par: " << xsLong_calc(E,Ep,theta,Q2,x,pt.g1,pt.g2) << endl;
    // cout << "sig_perp: " << xsTrans_calc(E,Ep,theta,Q2,pt.g1,pt.g2) << endl;
    // cout << "A_par: " << xsLong_calc(E,Ep,theta,Q2,x,pt.g1,pt.g2)/(2*xs_calc(E,Ep,theta,pt.F1,pt.F2)) << endl;
    // cout << "A_perp: " << xsTrans_calc(E,Ep,theta,Q2,pt.g1,pt.g2)/(2*xs_calc(E,Ep,theta,pt.F1,pt.F2)) << endl;

    // E = 5.73;
    // Ep = 1.44;
    // theta = 44.94;
    // Q2 = 4.831;
    // x = 0.601;
    // cout << "Q^2 = " << Q2 << "GeV^2 \t" << "x = " << x << endl;
    // pt.setEqual(interpolate2DLinear(x, Q2, dataPoints, q2_keys));
    // cout << "F1: " << pt.F1 << endl;
    // cout << "F2: " << pt.F2 << endl;
    // cout << "G1: " << pt.g1 << endl;
    // cout << "G2: " << pt.g2 << endl;
    // cout << "sig_unpol: " << xs_calc(E,Ep,theta,pt.F1,pt.F2) << endl;
    // cout << "sig_par: " << xsLong_calc(E,Ep,theta,Q2,x,pt.g1,pt.g2) << endl;
    // cout << "sig_perp: " << xsTrans_calc(E,Ep,theta,Q2,pt.g1,pt.g2) << endl;
    // cout << "A_par: " << xsLong_calc(E,Ep,theta,Q2,x,pt.g1,pt.g2)/(2*xs_calc(E,Ep,theta,pt.F1,pt.F2)) << endl;
    // cout << "A_perp: " << xsTrans_calc(E,Ep,theta,Q2,pt.g1,pt.g2)/(2*xs_calc(E,Ep,theta,pt.F1,pt.F2)) << endl;

    // return 0;

    ofstream outputFile;
    string outputFileName = "../../data_strufunc_jam.dat";
    outputFile.open(outputFileName);

    int count = 0;

    std::random_device rd;
    std::mt19937 gen(rd()); // Mersenne Twister PRNG engine
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double random_number = distribution(gen);

    double dQ2 = 0.001, Q2 = q2min, dx = 0.00001, x, q2low, q2high, xlow, xhigh;
    Point pt;
    cout << endl;
    // Q2 Loop
    if (false) {
        Q2 = (q2max-q2min)*random_number + q2min;
        q2low = 0;
        for (auto q2_key : q2_keys) {
            if (q2_key <= Q2) {
                q2low = q2_key;
                continue;
            }
            break;
        }
        x = dataPoints[q2low][0].x;
        cout << "Q^2 = " << Q2 << " GeV^2" << endl;
        while (x <= dataPoints[q2low][dataPoints[q2low].size()-1].x) {
            cout << std::fixed << std::setprecision(2) << 
            ((x-dataPoints[q2low][0].x)/(dataPoints[q2low][dataPoints[q2low].size()-1].x-dataPoints[q2low][0].x))*100 
            << "%" << "\r";
            cout.flush();
            pt = interpolate2DLinear(x, Q2, dataPoints, q2_keys);
            outputFile << x << "\t" << Q2 << "\t" << pt.F1 << "\t" << pt.F2 << "\t" << pt.g1 << "\t" << pt.g2 << endl;
            x += dx;
        }
        Q2 += dQ2;
    }

    if (true) {
        x = random_number;
        Q2 = q2min;
        cout << "x = " << x << endl;
        while (Q2 < q2max) {
            Q2 += dQ2;
            // cout << Q2 << endl;
            for (auto q2_key : q2_keys) {
                if (q2_key <= Q2) q2low = q2_key;
                if (q2_key > Q2) {
                    q2high = q2_key;
                    break;
                }
            }
            xlow = max(dataPoints[q2low][0].x,dataPoints[q2high][0].x);
            xhigh = min(dataPoints[q2low][dataPoints[q2low].size()-1].x,dataPoints[q2high][dataPoints[q2high].size()-1].x);
            if (x < xlow || x > xhigh) continue;
            cout << std::fixed << std::setprecision(2) << ((Q2-q2min)/(q2max-q2min))*100 << "%" << "\r";
            cout.flush();
            pt = interpolate2DLinear(x, Q2, dataPoints, q2_keys);
            outputFile << x << "\t" << Q2 << "\t" << pt.F1 << "\t" << pt.F2 << "\t" << pt.g1 << "\t" << pt.g2 << endl;
        }

    }

    outputFile.close();
    cout << "Data saved in " << outputFileName << endl;


    // // Theta Loop
    // for (double theta_value : theta_values) {
    //     theta = theta_value;
    //     thetalabel = int(theta);
    //     E = 10.38;
    //     Elabel = int(E*1000);
    //     cout << endl << "Theta: " << fixed << setprecision(0) << theta << " degrees" << endl;
    //     // Beam Energy Loop
    //     for (int i=0; i < num_Ebeams; i++) {
    //         count = 0;
    //         countL = 0;
    //         countT = 0;
    //         if (i != 0) { 
    //             E -= dE;
    //             Elabel -= int(dE*1000);
    //         }
    //         outputFileName = "Data/Unpol/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg.dat";
    //         outputFile.open(outputFileName);
    //         outputLongFileName = "Data/Long/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_long.dat";
    //         outputLongFile.open(outputLongFileName);
    //         outputTransFileName = "Data/Trans/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_trans.dat";
    //         outputTransFile.open(outputTransFileName);
    //         inputFilePath = "../../CAnalyzer-master/example/data/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg.dat";
    //         inputLongFilePath = "../../CAnalyzer-master/example/data/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_long.dat";
    //         inputTransFilePath = "../../CAnalyzer-master/example/data/E"+to_string(Elabel)+"_"+to_string(thetalabel)+"deg_trans.dat";
    //         Ep = dEp;
    //         Q2 = q2_calc(E, Ep, theta);
    //         x = x_calc(E, Ep, theta);
    //         w = w_calc(E, Ep, theta);
    //         cout << endl << "Ebeam: " << fixed << setprecision(2) << E << " GeV" << endl;
    //         while(Ep < E) {
    //             // cout << x << "\t" << Q2 << "\t" << Ep << endl;
    //             cout << std::fixed << std::setprecision(4) << Ep/E*100 << "%" << "\r";
    //             cout.flush();
    //             // if (x < 1.0 && w > Mp) {
    //             double q2low = 0;
    //             for (auto q2_key : q2_keys) {
    //                 if (q2_key <= Q2) {
    //                     q2low = q2_key;
    //                     continue;
    //                 }
    //                 break;
    //             }
    //             if (Q2 >= q2min && Q2 <= q2max && x <= dataPoints[q2low][dataPoints[q2low].size()-1].x && w >= Mp) {
    //                 pt = interpolate2DLinear(x, Q2, dataPoints, q2_keys);
    //                 xsUnpol = xs_calc(E, Ep, theta, pt.F1, pt.F2);
    //                 xsLong = xsLong_calc(E, Ep, theta, Q2, x, pt.g1, pt.g2);
    //                 xsTrans = xsTrans_calc(E, Ep, theta, Q2, pt.g1, pt.g2);
    //                 bool jumpUnpol = false, jumpLong = false, jumpTrans = false;
    //                 if (!(xsUnpol_prev == 0 && xsUnpol != 0)) jumpUnpol = fabs(xsUnpol) > 5e2; //jumpUnpol = (fabs(xsUnpol_prev - xsUnpol) > 1e5*min(fabs(xsUnpol_prev),fabs(xsUnpol)));
    //                 if (!(xsLong_prev == 0 && xsLong != 0)) jumpLong = fabs(xsLong) > 5e2; //jumpLong = (fabs(xsLong_prev - xsLong) > 1e5*min(fabs(xsLong_prev),fabs(xsLong)));
    //                 if (!(xsTrans_prev == 0 && xsTrans != 0)) jumpTrans = fabs(xsTrans) > 5e2; //jumpTrans = (fabs(xsTrans_prev - xsTrans) > 1e5*min(fabs(xsTrans_prev),fabs(xsTrans)));

    //                 // nu = (E - Ep)*1000;     // input in MeV instead of GeV
    //                 // if (!isnan(xsUnpol) && !jumpUnpol) {
    //                 //     outputFile << nu << "\t" << xsUnpol << "\t" << 0.0 << "\t" << 0.0 << endl;
    //                 //     count++;
    //                 //     xsUnpol_prev = xsUnpol;                        
    //                 // }
    //                 // if (!isnan(xsLong) && !jumpLong) {
    //                 //     outputLongFile << nu << "\t" << xsLong << "\t" << 0.0 << "\t" << 0.0 << endl;
    //                 //     countL++;
    //                 //     xsLong_prev = xsLong;
    //                 // }
    //                 // if (!isnan(xsTrans) && !jumpTrans) {
    //                 //     outputTransFile << nu << "\t" << xsTrans << "\t" << 0.0 << "\t" << 0.0 << endl;
    //                 //     countT++;
    //                 //     xsTrans_prev = xsTrans;
    //                 // }

    //                 bool jump = (jumpUnpol || jumpLong || jumpTrans);
    //                 if (!(isnan(xsUnpol) || isnan(xsLong) || isnan(xsTrans) || xsUnpol <= 0) && !jump) {
    //                 // if (true) {
    //                     nu = (E - Ep)*1000;     // input in MeV instead of GeV
    //                     // if (xsUnpol > 0) {
    //                         outputFile << nu << "\t" << xsUnpol << "\t" << 0.0 << "\t" << 0.0 << endl;
    //                         count++;
    //                     // }
    //                     // if (xsLong > 0) { 
    //                         outputLongFile << nu << "\t" << xsLong << "\t" << 0.0 << "\t" << 0.0 << endl;
    //                         countL++;
    //                     // }
    //                     // if (xsTrans > 0) {
    //                         outputTransFile << nu << "\t" << xsTrans << "\t" << 0.0 << "\t" << 0.0 << endl;
    //                         countT++;
    //                     // }
    //                     xsUnpol_prev = xsUnpol;
    //                     xsLong_prev = xsLong;
    //                     xsTrans_prev = xsTrans;
    //                 }
    //             }
    //             Ep += dEp;
    //             Q2 = q2_calc(E, Ep, theta);
    //             x = x_calc(E, Ep, theta);
    //             w = w_calc(E, Ep, theta);
    //         }
    //         outputFile.close();
    //         cout << count << " data points saved in " << outputFileName << endl;
    //         outputLongFile.close();
    //         cout << countL << " data points saved in" << outputLongFileName << endl;
    //         outputTransFile.close();
    //         cout << countT << " data points saved in" << outputTransFileName << endl;

    //         ifstream sourceFile(outputFileName, ios::binary);
    //         ofstream copyFile(inputFilePath, ios::binary);
    //         if (!copyFile) {
    //             cout << "Error opening " << inputFilePath << endl;
    //         }
    //         copyFile << sourceFile.rdbuf();
    //         sourceFile.close();
    //         copyFile.close();
    //         cout << "Data copied successfully to " << inputFilePath << endl;

    //         ifstream sourceLongFile(outputLongFileName, ios::binary);
    //         ofstream copyLongFile(inputLongFilePath, ios::binary);
    //         if (!copyLongFile) {
    //             cout << "Error opening " << inputLongFilePath << endl;
    //         }
    //         copyLongFile << sourceLongFile.rdbuf();
    //         sourceLongFile.close();
    //         copyLongFile.close();
    //         cout << "Data copied successfully to " << inputLongFilePath << endl;

    //         ifstream sourceTransFile(outputTransFileName, ios::binary);
    //         ofstream copyTransFile(inputTransFilePath, ios::binary);
    //         if (!copyTransFile) {
    //             cout << "Error opening " << inputTransFilePath << endl;
    //         }
    //         copyTransFile << sourceTransFile.rdbuf();
    //         sourceTransFile.close();
    //         copyTransFile.close();
    //         cout << "Data copied successfully to " << inputTransFilePath << endl;
    //     }
    // }
}
