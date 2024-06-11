#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

vector<pair<double, double>> calcExtremePoint(double, double, double, double, double, double);
vector<double> solveQuadraticEquation(double, double, double);

int main()
{
    vector<pair<double, double>> extreme_points0 = calcExtremePoint(3.0, 2.0, 0.0, 5.0, 1.0, -4.0);
    for (auto s : extreme_points0) {
        cout << ".(" << s.first << ", " << s.second << ")" << endl;
    }
    system("pause");
}

vector<pair<double, double>> calcExtremePoint(double A, double B, double C, double D, double E, double F)
{
    vector<double> solution = solveQuadraticEquation(4 * A * A * B - A * C * C, 4 * A * B * D - 2 * A * C * E, B * D * D - C * D * E + C * C * F);
    pair<double, double> result;
    vector<pair<double, double>> extreme_points;
    if (solution.size() == 2) {
        result.first = solution.at(0);
        result.second = -(2 * A * result.first + D) / C;
        extreme_points.push_back(result);
        result.first = solution.at(1);
        result.second = -(2 * A * result.first + D) / C;
        extreme_points.push_back(result);
    } else if (solution.size() == 1) {
        result.first = solution.at(0);
        result.second = sqrt((-4 * A * B * result.first * result.first - 4 * B * D * result.first - 4 * B * F + E * E) / (4 * B * B)) - E / (2 * B);
        extreme_points.push_back(result);
        result.second = -sqrt((-4 * A * B * result.first * result.first - 4 * B * D * result.first - 4 * B * F + E * E) / (4 * B * B)) - E / (2 * B);
        extreme_points.push_back(result);
    }
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
