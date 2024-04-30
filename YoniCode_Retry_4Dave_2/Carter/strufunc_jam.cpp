#include "strufunc_jam.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <sstream>

using namespace std;

// map<double, vector<Point>> dataPoints;
// vector<double> q2_keys;

void getDataPoints(map<double, vector<Point>> &dataPoints, vector<double> &q2_keys) {
    q2_keys.clear();
    dataPoints.clear();

    string jamFileName = "table_3He_JAM_smeared_no_QE_ipol1IA14_SF23_AC11.csv";
    // string jamFileName = "table_3He_JAM_smeared_inelastic_plus_QE_ipol1IA14_SF23_AC11.csv";
    ifstream jamFile(jamFileName);
    string line;
    while (getline(jamFile, line)) {
        stringstream ss(line);
        Point pt;
        char comma;
        double other;

        // used for "table_3He_JAM_smeared_no_QE_ipol1IA14_SF23_AC11.csv"
        ss >> pt.q2 >> comma >> pt.x;
        for (int i = 0; i < 10; i++) ss >> comma >> other;
        ss >> comma >> pt.F2;
        for (int i = 0; i < 3; i++) ss >> comma >> other;
        ss >> comma >> pt.g1 >> comma >> pt.g2 >> comma >> pt.F1;
        
        // // used for "table_3He_JAM_smeared_inelastic_plus_QE_ipol1IA14_SF23_AC11.csv"
        // ss >> pt.q2 >> comma >> pt.x >> comma >> pt.F2 >> comma >> other >> comma >> pt.F1 >> comma;
        // ss >> pt.g1 >> comma >> pt.g2;

        dataPoints[pt.q2].push_back(pt);
        bool q2inkeys = false;
        for (auto key : q2_keys) {
            if (key == pt.q2) {
                q2inkeys = true;
                break;
            }
        }
        if (!q2inkeys) q2_keys.push_back(pt.q2);
    }
    sort(q2_keys.begin(),q2_keys.end());
    jamFile.close();
    return;
}

Point interpolate2DLinear(double x, double q2, map<double, vector<Point>> &dataPoints, vector<double> q2_keys) {
    Point pt;

    double q2min = 0, q2max = 0, xmin1 = 0, xmax1 = 0, xmin2 = 0, xmax2 = 0, F1 = 0, F2 = 0, g1 = 0, g2 = 0;
    int ixmin1 = 0, ixmax1 = 0, ixmin2 = 0, ixmax2 = 0;
    for (auto q2_key : q2_keys) {
        if (q2_key <= q2) q2min = q2_key;
        if (q2_key > q2) { 
            q2max = q2_key;
            break;
        }
    }
    for (int i = 0; i < dataPoints[q2min].size(); i++) {
        if (dataPoints[q2min][i].x <= x) ixmin1 = i;
        if (dataPoints[q2min][i].x > x) {
            ixmax1 = i;
            break;
        }
    }
    for (int i = 0; i < dataPoints[q2max].size(); i++) {
        if (dataPoints[q2max][i].x <= x) ixmin2 = i;
        if (dataPoints[q2max][i].x > x) {
            ixmax2 = i;
            break;
        }
    }

    xmin1 = dataPoints[q2min][ixmin1].x;
    xmax1 = dataPoints[q2min][ixmax1].x;
    xmin2 = dataPoints[q2max][ixmin2].x;
    xmax2 = dataPoints[q2max][ixmax2].x;

    F1 = (q2max - q2)/(q2max - q2min) * ((xmax2 - x)/(xmax2 - xmin2)*dataPoints[q2min][ixmin1].F1 + (x - xmin2)/(xmax2 - xmin2)*dataPoints[q2min][ixmax1].F1) 
         + (q2 - q2min)/(q2max - q2min) * ((xmax1 - x)/(xmax1 - xmin1)*dataPoints[q2max][ixmin1].F1 + (x - xmin1)/(xmax1 - xmin1)*dataPoints[q2max][ixmax1].F1);

    F2 = (q2max - q2)/(q2max - q2min) * ((xmax2 - x)/(xmax2 - xmin2)*dataPoints[q2min][ixmin1].F2 + (x - xmin2)/(xmax2 - xmin2)*dataPoints[q2min][ixmax1].F2) 
         + (q2 - q2min)/(q2max - q2min) * ((xmax1 - x)/(xmax1 - xmin1)*dataPoints[q2max][ixmin1].F2 + (x - xmin1)/(xmax1 - xmin1)*dataPoints[q2max][ixmax1].F2);

    g1 = (q2max - q2)/(q2max - q2min) * ((xmax2 - x)/(xmax2 - xmin2)*dataPoints[q2min][ixmin1].g1 + (x - xmin2)/(xmax2 - xmin2)*dataPoints[q2min][ixmax1].g1) 
         + (q2 - q2min)/(q2max - q2min) * ((xmax1 - x)/(xmax1 - xmin1)*dataPoints[q2max][ixmin1].g1 + (x - xmin1)/(xmax1 - xmin1)*dataPoints[q2max][ixmax1].g1);

    g2 = (q2max - q2)/(q2max - q2min) * ((xmax2 - x)/(xmax2 - xmin2)*dataPoints[q2min][ixmin1].g2 + (x - xmin2)/(xmax2 - xmin2)*dataPoints[q2min][ixmax1].g2) 
         + (q2 - q2min)/(q2max - q2min) * ((xmax1 - x)/(xmax1 - xmin1)*dataPoints[q2max][ixmin1].g2 + (x - xmin1)/(xmax1 - xmin1)*dataPoints[q2max][ixmax1].g2);

    pt.setVals(x, q2, F1, F2, g1, g2);
    return pt;
}


