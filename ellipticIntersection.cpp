#include <algorithm>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

constexpr double ZERO = 1e-6;

vector<pair<double, double>> getEllipse(double &, double &, double &, double &, double &, double &);
vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> getArcsPairs(vector<pair<double, double>>, vector<pair<double, double>>);
vector<pair<pair<double, double>, pair<double, double>>> getArcs(vector<pair<double, double>>);
bool isEllipse(double, double, double, double, double, double);
vector<pair<double, double>> calcExtremePoint(double, double, double, double, double, double);
vector<double> solveQuadraticEquation(double, double, double);
vector<pair<double, double>> sortExtremePoint(vector<pair<double, double>>);
bool isObviouslynotintersect(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>);
void calcIntersection(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>, vector<pair<double, double>> &, double, double, double, double, double, double, double, double, double, double, double, double);
double getFx(double, double, double, double, double, double, double, double, double);
vector<double> getRectangle(pair<pair<double, double>, pair<double, double>>);
void addIntersection(vector<pair<double, double>> &, pair<double, double>);

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
        calcIntersection(s, intersections, A1, B1, C1, D1, E1, F1, A2, B2, C2, D2, E2, F2);
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
    if (I2 > ZERO && abs(I3) > ZERO && I1 * I3 < -ZERO) {
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
    if (delta > ZERO) {
        result.push_back((-b + sqrt(delta)) / (2 * a));
        result.push_back((-b - sqrt(delta)) / (2 * a));
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

void calcIntersection(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double A0, double B0, double C0, double D0, double E0, double F0, double A1, double B1, double C1, double D1, double E1, double F1)
{
    double x0, x1, FP0, FP1, FQ0, FQ1, delta0, delta1;
    x0 = max(min(arcs_pair.first.first.first, arcs_pair.first.second.first), min(arcs_pair.second.first.first, arcs_pair.second.second.first));
    x1 = min(max(arcs_pair.first.first.first, arcs_pair.first.second.first), max(arcs_pair.second.first.first, arcs_pair.second.second.first));
    FP0 = getFx(x0, arcs_pair.first.first.second, arcs_pair.first.second.second, A0, B0, C0, D0, E0, F0);
    FQ0 = getFx(x0, arcs_pair.second.first.second, arcs_pair.second.second.second, A1, B1, C1, D1, E1, F1);
    FP1 = getFx(x1, arcs_pair.first.first.second, arcs_pair.first.second.second, A0, B0, C0, D0, E0, F0);
    FQ1 = getFx(x1, arcs_pair.second.first.second, arcs_pair.second.second.second, A1, B1, C1, D1, E1, F1);
    delta0 = FP0 - FQ0;
    delta1 = FP1 - FQ1;
    if (abs(delta0) <= ZERO) {
        addIntersection(intersections, make_pair(x0, FP0));
    } else if (abs(delta1) <= ZERO) {
        addIntersection(intersections, make_pair(x1, FP1));
    }
}

double getFx(double x, double ymin, double ymax, double A, double B, double C, double D, double E, double F)
{
    vector<double> solution;
    solution = solveQuadraticEquation(B, C * x + E, A * x * x + D * x + F);
    if (ymin > ymax) {
        swap(ymin, ymax);
    }
    if (solution.size() == 2) {
        if (solution.at(0) >= ymin && solution.at(1) <= ymax) {
            return solution.at(0);
        } else {
            return solution.at(1);
        }
    } else {
        return solution.at(0);
    }
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
