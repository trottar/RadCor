#ifndef strufunc_jam_H
#define strufunc_jam_H
#include <vector>
#include <map>

using namespace std;

class Point {
    public:
    double x, q2, F1, F2, g1, g2;
    Point() : x(0), q2(0), F1(0), F2(0), g1(0), g2(0) {};
    Point(double x, double q2, double F1, double F2, double g1, double g2) : x(x), q2(q2), F1(F1), F2(F2), g1(g1), g2(g2) {};
    void setVals(double x, double q2, double F1, double F2, double g1, double g2) {
        this->x = x;
        this->q2 = q2;
        this->F1 = F1;
        this->F2 = F2;
        this->g1 = g1;
        this->g2 = g2;
    }
    void setEqual(Point pt) {
        this->setVals(pt.x, pt.q2, pt.F1, pt.F2, pt.g1, pt.g2);
    }

};

void getDataPoints(map<double, vector<Point>> &dataPoints, vector<double> &q2_keys);
Point interpolate2DLinear(double x, double q2, map<double, vector<Point>> &dataPoints, vector<double>);

#endif
