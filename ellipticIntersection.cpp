#define _USE_MATH_DEFINES
#include <algorithm>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

constexpr double ZERO = 1e-6;
constexpr double Y_ACCURACY_THRESHOLD = 8e-4;
constexpr double SK_ACCURACY_THRESHOLD = 1e-6;

vector<pair<double, double>> getEllipse(double &A, double &B, double &C, double &D, double &E, double &F);
vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> getArcsPairs(vector<pair<double, double>> extreme_points0, vector<pair<double, double>> extreme_points1);
vector<pair<pair<double, double>, pair<double, double>>> getArcs(vector<pair<double, double>> extreme_points);
bool isEllipse(double A, double B, double C, double D, double E, double F);
vector<pair<double, double>> calcExtremePoint(double A, double B, double C, double D, double E, double F);
vector<double> solveQuadraticEquation(double a, double b, double c);
vector<pair<double, double>> sortExtremePoint(vector<pair<double, double>> extreme_points);
bool isObviouslynotintersect(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pairs);
void calcIntersection(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double A0, double B0, double C0, double D0, double E0, double F0,vector<pair<double, double>> extreme_points0, double A1, double B1, double C1, double D1, double E1, double F1, vector<pair<double, double>> extreme_points1);
double getFx(double x, pair<pair<double, double>, pair<double, double>> arc, double A, double B, double C, double D, double E, double F);
int getArcNo(pair<pair<double, double>, pair<double, double>> arc);
double getAngle(double x, double y, double left, pair<pair<double, double>, pair<double, double>> arc, double A, double B, double C, double D, double E, double F);
vector<double> getRectangle(pair<pair<double, double>, pair<double, double>> arc);
void addIntersection(vector<pair<double, double>> &intersections, pair<double, double> point);
void XLessThan0(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double x0, double x1, double delta0, double delta1, double A0, double B0, double C0, double D0, double E0, double F0, double A1, double B1, double C1, double D1, double E1, double F1);
double xkWhenSkEqual0(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, double x0, double x1, double a0, double b0, double A0, double B0, double C0, double D0, double E0, double F0, vector<pair<double, double>> extreme_points0, double A1, double B1, double C1, double D1, double E1, double F1, vector<pair<double, double>> extreme_points1);

int main()
{
    double A1, B1, C1, D1, E1, F1, A2, B2, C2, D2, E2, F2;
    vector<pair<double, double>> extreme_points0, extreme_points1, intersections;
    vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> arcs_pairs;
    extreme_points0 = getEllipse(A1, B1, C1, D1, E1, F1);
    extreme_points1 = getEllipse(A2, B2, C2, D2, E2, F2);
    arcs_pairs = getArcsPairs(extreme_points0, extreme_points1);
    for (auto s : arcs_pairs) {
        if (isObviouslynotintersect(s)) {
            continue;
        }
        calcIntersection(s, intersections, A1, B1, C1, D1, E1, F1, extreme_points0, A2, B2, C2, D2, E2, F2, extreme_points1);
    }
    cout << "两椭圆的交点坐标为\n";
    for (auto s : intersections) {
        cout << "(" << (abs(s.first) <= ZERO ? 0 : s.first) << ", " << (abs(s.second) <= ZERO ? 0 : s.second) << ")\n";
    }
    system("pause");
}

vector<pair<double, double>> getEllipse(double &A, double &B, double &C, double &D, double &E, double &F)
{
    static int cont = 1;
    vector<pair<double, double>> extreme_points;
    cout << "请输入第" << cont << "个椭圆的方程Ax^2+By^2+Cxy+Dx+Ey+F=0中的参数，以空格为分隔\n";
    while (true) {
        cin >> A >> B >> C >> D >> E >> F;
        cin.sync();
        if (cin.rdstate() != 0) {
            cout << "输入异常，请重新输入\n";
            cin.clear();
            continue;
        }
        if (!isEllipse(A, B, C, D, E, F)) {
            cout << "参数无法构成椭圆，请重新输入\n";
            continue;
        }
        break;
    }
    extreme_points = calcExtremePoint(A, B, C, D, E, F);
    extreme_points = sortExtremePoint(extreme_points);
    cont++;
    return extreme_points;
}

vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> getArcsPairs(vector<pair<double, double>> extreme_points0, vector<pair<double, double>> extreme_points1)
{
    vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> arcs_pairs;
    vector<pair<pair<double, double>, pair<double, double>>> arcs0, arcs1;
    arcs0 = getArcs(extreme_points0);
    arcs1 = getArcs(extreme_points1);
    for (auto iter0 = arcs0.begin(); iter0 != arcs0.end(); iter0++) {
        for (auto iter1 = arcs1.begin(); iter1 != arcs1.end(); iter1++) {
            arcs_pairs.push_back(make_pair(*iter0, *iter1));
        }
    }
    return arcs_pairs;
}

vector<pair<pair<double, double>, pair<double, double>>> getArcs(vector<pair<double, double>> extreme_points)
{
    vector<pair<pair<double, double>, pair<double, double>>> arcs;
    for (auto iter = extreme_points.begin(); iter != extreme_points.end() - 1; iter++) {
        arcs.push_back(make_pair(*iter, *(iter + 1)));
    }
    arcs.push_back(make_pair(extreme_points.back(), extreme_points.at(0)));
    return arcs;
}

bool isEllipse(double A, double B, double C, double D, double E, double F)
{
    double P = D / 2, Q = E / 2, I1, I2, I3;
    swap(B, C);
    B /= 2;
    I1 = A + C;
    I2 = A * C - B * B;
    I3 = A * C * F + B * Q * P + P * B * Q - P * C * P - B * B * F - A * Q * Q;
    if (I2 > 0 && abs(I3) > 0 && I1 * I3 < 0) {
        return true;
    } else {
        return false;
    }
}

vector<pair<double, double>> calcExtremePoint(double A, double B, double C, double D, double E, double F)
{
    vector<double> solution = solveQuadraticEquation(4 * A * A * B - A * C * C, 4 * A * B * D - 2 * A * C * E, B * D * D - C * D * E + C * C * F);
    vector<pair<double, double>> extreme_points;
    if (solution.size() == 2) {
        extreme_points.push_back(make_pair(solution.at(0), -(2 * A * solution.at(0) + D) / C));
        extreme_points.push_back(make_pair(solution.at(1), -(2 * A * solution.at(1) + D) / C));
    } else if (solution.size() == 1) {
        extreme_points.push_back(make_pair(solution.at(0), sqrt((-4 * A * B * solution.at(0) * solution.at(0) - 4 * B * D * solution.at(0) - 4 * B * F + E * E) / (4 * B * B)) - E / (2 * B)));
        extreme_points.push_back(make_pair(solution.at(0), -sqrt((-4 * A * B * solution.at(0) * solution.at(0) - 4 * B * D * solution.at(0) - 4 * B * F + E * E) / (4 * B * B)) - E / (2 * B)));
    }
    solution = solveQuadraticEquation(4 * A * B - C * C, 4 * B * D - 2 * C * E, 4 * B * F - E * E);
    extreme_points.push_back(make_pair(solution.at(0), -(C * solution.at(0) + E) / (2 * B)));
    extreme_points.push_back(make_pair(solution.at(1), -(C * solution.at(1) + E) / (2 * B)));
    return extreme_points;
}

vector<double> solveQuadraticEquation(double a, double b, double c)
{
    vector<double> result;
    double delta = b * b - 4 * a * c;
    if (delta > 0) {
        result.push_back((-b + sqrt(delta)) / (2 * a));
        result.push_back((-b - sqrt(delta)) / (2 * a));
        sort(result.begin(), result.end(), [](double a, double b) { return a > b; });
    } else if (abs(delta) <= ZERO) {
        result.push_back(-b / (2 * a));
    }
    return result;
}

vector<pair<double, double>> sortExtremePoint(vector<pair<double, double>> extreme_points)
{
    iter_swap(extreme_points.begin(), max_element(extreme_points.begin(), extreme_points.end(), [](pair<double, double> a, pair<double, double> b) { return a.second < b.second; }));
    iter_swap(extreme_points.begin() + 1, max_element(extreme_points.begin() + 1, extreme_points.end(), [](pair<double, double> a, pair<double, double> b) { return a.first < b.first; }));
    iter_swap(extreme_points.begin() + 2, min_element(extreme_points.begin() + 2, extreme_points.end(), [](pair<double, double> a, pair<double, double> b) { return a.second < b.second; }));
    return extreme_points;
}

bool isObviouslynotintersect(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pairs)
{
    vector<double> rectangle0, rectangle1;
    rectangle0 = getRectangle(arcs_pairs.first);
    rectangle1 = getRectangle(arcs_pairs.second);
    if (rectangle0.at(0) < rectangle1.at(2) || rectangle0.at(1) < rectangle1.at(3) || rectangle0.at(2) > rectangle1.at(0) || rectangle0.at(3) > rectangle1.at(1)) {
        return true;
    } else {
        return false;
    }
}

void calcIntersection(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double A0, double B0, double C0, double D0, double E0, double F0, vector<pair<double, double>> extreme_points0, double A1, double B1, double C1, double D1, double E1, double F1, vector<pair<double, double>> extreme_points1)
{
    double x0, x1, FP0, FP1, FQ0, FQ1, delta0, delta1, X, xk, FPk, FQk, deltak, a0, a1, b0, b1;
    x0 = max(min(arcs_pair.first.first.first, arcs_pair.first.second.first), min(arcs_pair.second.first.first, arcs_pair.second.second.first));
    x1 = min(max(arcs_pair.first.first.first, arcs_pair.first.second.first), max(arcs_pair.second.first.first, arcs_pair.second.second.first));
    FP0 = getFx(x0, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQ0 = getFx(x0, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    FP1 = getFx(x1, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQ1 = getFx(x1, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    delta0 = FP0 - FQ0;
    delta1 = FP1 - FQ1;
    if (abs(delta0) <= ZERO) {
        addIntersection(intersections, make_pair(x0, FP0));
        return;
    } else if (abs(delta1) <= ZERO) {
        addIntersection(intersections, make_pair(x1, FP1));
        return;
    }
    X = delta0 * delta1;
    if (X < 0) {
        XLessThan0(arcs_pair, intersections, x0, x1, delta0, delta1, A0, B0, C0, D0, E0, F0, A1, B1, C1, D1, E1, F1);
        return;
    }
    a0 = getAngle(x0, FP0, extreme_points0.at(3).first, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    a1 = getAngle(x1, FP1, extreme_points0.at(3).first, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    b0 = getAngle(x0, FQ0, extreme_points1.at(3).first, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    b1 = getAngle(x1, FQ1, extreme_points1.at(3).first, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    if ((a0 - b0) * (a1 - b1) > 0) {
        return;
    }
    xk = xkWhenSkEqual0(arcs_pair, x0, x1, a0, b0, A0, B0, C0, D0, E0, F0, extreme_points0, A1, B1, C1, D1, E1, F1, extreme_points1);
    FPk = getFx(xk, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQk = getFx(xk, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    deltak = FPk - FQk;
    if (abs(deltak) <= ZERO) {
        addIntersection(intersections, make_pair(xk, FPk));
        return;
    }
    if (delta0 * deltak > 0) {
        return;
    }
    XLessThan0(arcs_pair, intersections, x0, xk, delta0, deltak, A0, B0, C0, D0, E0, F0, A1, B1, C1, D1, E1, F1);
    XLessThan0(arcs_pair, intersections, xk, x1, deltak, delta1, A0, B0, C0, D0, E0, F0, A1, B1, C1, D1, E1, F1);
}

double getFx(double x, pair<pair<double, double>, pair<double, double>> arc, double A, double B, double C, double D, double E, double F)
{
    int arc_no = getArcNo(arc);
    vector<double> solution;
    solution = solveQuadraticEquation(B, C * x + E, A * x * x + D * x + F);
    if (solution.size() == 2) {
        if (arc_no == 1 || arc_no == 4) {
            return solution.at(0);
        } else {
            return solution.at(1);
        }
    } else {
        return solution.at(0);
    }
}

int getArcNo(pair<pair<double, double>, pair<double, double>> arc)
{
    if (arc.first.first < arc.second.first && arc.first.second > arc.second.second) {
        return 1;
    } else if (arc.first.first > arc.second.first && arc.first.second > arc.second.second) {
        return 2;
    } else if (arc.first.first > arc.second.first && arc.first.second < arc.second.second) {
        return 3;
    } else if (arc.first.first < arc.second.first && arc.first.second < arc.second.second) {
        return 4;
    }
    return -1;
}

double getAngle(double x, double y, double left, pair<pair<double, double>, pair<double, double>> arc, double A, double B, double C, double D, double E, double F)
{
    int arc_no = getArcNo(arc);
    double derivative, radian;
    if (abs(C * x + 2 * B * y + E) <= ZERO) {
        if (abs(x - left) <= ZERO) {
            if (arc_no == 3) {
                return -90;
            } else {
                return 90;
            }
        } else {
            if (arc_no == 1) {
                return -90;
            } else {
                return 90;
            }
        }
    }
    derivative = (-2 * A * x - C * y - D) / (C * x + 2 * B * y + E);
    radian = atan(derivative);
    return radian * 180 / M_PI;
}

vector<double> getRectangle(pair<pair<double, double>, pair<double, double>> arc)
{
    vector<double> rectangle;
    rectangle.push_back(max(arc.first.second, arc.second.second));
    rectangle.push_back(max(arc.first.first, arc.second.first));
    rectangle.push_back(min(arc.first.second, arc.second.second));
    rectangle.push_back(min(arc.first.first, arc.second.first));
    return rectangle;
}

void addIntersection(vector<pair<double, double>> &intersections, pair<double, double> point)
{
    if (find(intersections.begin(), intersections.end(), point) == intersections.end()) {
        intersections.push_back(point);
    }
}

void XLessThan0(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double x0, double x1, double delta0, double delta1, double A0, double B0, double C0, double D0, double E0, double F0, double A1, double B1, double C1, double D1, double E1, double F1)
{
    double xk, FPk, FQk, deltak;
    xk = (x0 + x1) / 2;
    FPk = getFx(xk, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQk = getFx(xk, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    deltak = FPk - FQk;
    if (abs(deltak) <= Y_ACCURACY_THRESHOLD) {
        addIntersection(intersections, make_pair(xk, FPk));
        return;
    } else if (delta0 * deltak < 0) {
        x1 = xk;
        delta1 = deltak;
    } else {
        x0 = xk;
        delta0 = deltak;
    }
    XLessThan0(arcs_pair, intersections, x0, x1, delta0, delta1, A0, B0, C0, D0, E0, F0, A1, B1, C1, D1, E1, F1);
}

double xkWhenSkEqual0(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, double x0, double x1, double a0, double b0, double A0, double B0, double C0, double D0, double E0, double F0, vector<pair<double, double>> extreme_points0, double A1, double B1, double C1, double D1, double E1, double F1, vector<pair<double, double>> extreme_points1)
{
    double xk, FPk, FQk, ak, bk, sk, s0;
    xk = (x0 + x1) / 2;
    FPk = getFx(xk, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQk = getFx(xk, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    ak = getAngle(xk, FPk, extreme_points0.at(3).first, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    bk = getAngle(xk, FQk, extreme_points1.at(3).first, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    s0 = a0 - b0;
    sk = ak - bk;
    if (abs(sk) <= SK_ACCURACY_THRESHOLD) {
        return xk;
    }
    if (s0 * sk < 0) {
        x1 = xk;
    } else {
        x0 = xk;
        a0 = ak;
        b0 = bk;
    }
    return xkWhenSkEqual0(arcs_pair, x0, x1, a0, b0 , A0, B0, C0, D0, E0, F0, extreme_points0, A1, B1, C1, D1, E1, F1, extreme_points1);
}
