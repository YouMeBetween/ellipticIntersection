#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

bool isEllipse(double, double, double, double, double, double);
vector<pair<double, double>> calcExtremePoint(double, double, double, double, double, double);
vector<double> solveQuadraticEquation(double, double, double);

int main()
{
    double a, b, c, d, e, f;
    cout << "请输入方程Ax^2+By^2+Cxy+Dx+Ey+F=0中的参数，以空格为分隔\n";
    cin >> a >> b >> c >> d >> e >> f;
    if (!isEllipse(a, b, c, d, e, f)) {
        cout << "参数无法构成椭圆\n";
        system("pause");
        return -1;
    }
    vector<pair<double, double>> extreme_points0 = calcExtremePoint(a, b, c, d, e, f);
    for (auto s : extreme_points0) {
        cout << ".(" << s.first << ", " << s.second << ")" << endl;
    }
    system("pause");
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
