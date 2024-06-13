#include <algorithm>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

vector<pair<double, double>> getEllipse(double &, double &, double &, double &, double &, double &);
vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> getArcsPairs(vector<pair<double, double>>, vector<pair<double, double>>);
vector<pair<pair<double, double>, pair<double, double>>> getArcs(vector<pair<double, double>>);
bool isEllipse(double, double, double, double, double, double);
vector<pair<double, double>> calcExtremePoint(double, double, double, double, double, double);
vector<double> solveQuadraticEquation(double, double, double);
vector<pair<double, double>> sortExtremePoint(vector<pair<double, double>>);
bool isObviouslynotintersect(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>);
vector<double> getRectangle(pair<pair<double, double>, pair<double, double>>);

int main()
{
    double A1, B1, C1, D1, E1, F1, A2, B2, C2, D2, E2, F2;
    vector<pair<double, double>> extreme_points0, extreme_points1;
    vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> arcs_pairs;
    extreme_points0 = getEllipse(A1, B1, C1, D1, E1, F1);
    extreme_points1 = getEllipse(A2, B2, C2, D2, E2, F2);
    arcs_pairs = getArcsPairs(extreme_points0, extreme_points1);
    cout << "两个椭圆构成的弧对为\n";
    for (auto s : arcs_pairs) {
        if (isObviouslynotintersect(s)) {
            continue;
        }
        cout << "[(" << s.first.first.first << ", " << s.first.first.second << "), (" << s.first.second.first << ", " << s.first.second.second << ")], [(" << s.second.first.first << ", " << s.second.first.second << "), (" << s.second.second.first << ", " << s.second.second.second<< ")]\n";
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
    if (I2 > 1e-6 && abs(I3) > 1e-6 && I1 * I3 < -1e-6) {
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
    if (delta > 1e-6) {
        result.push_back((-b + sqrt(delta)) / (2 * a));
        result.push_back((-b - sqrt(delta)) / (2 * a));
    } else if (abs(delta) <= 1e-6) {
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

vector<double> getRectangle(pair<pair<double, double>, pair<double, double>> arc)
{
    vector<double> rectangle;
    rectangle.push_back(max(arc.first.second, arc.second.second));
    rectangle.push_back(max(arc.first.first, arc.second.first));
    rectangle.push_back(min(arc.first.second, arc.second.second));
    rectangle.push_back(min(arc.first.first, arc.second.first));
    return rectangle;
}
