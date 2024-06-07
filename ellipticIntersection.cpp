#include <iostream>
#include <ostream>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

int main()
{
    double a1 = 3.0, b1 = 2.0, c1 = 4.8, d1 = 5.0, e1 = 1.0, f1 = -4.0;
    double A = 4 * pow(a1, 2) * b1 - a1 * pow(c1, 2);
    double B = 4 * a1 * b1 * d1 - 2 * a1 * c1 * e1;
    double C = b1 * pow(d1, 2) - c1 * d1 * e1 + pow(c1, 2) * f1;
    double delta = pow(B, 2) - 4 * A * C;
    vector<double> x1, y1;
    if (delta > 0) {
        x1.push_back((-B + sqrt(delta)) / (2 * A));
        x1.push_back((-B - sqrt(delta)) / (2 * A));
    } else if (delta == 0) {
        x1.push_back(-B / (2 * A));
    }
    for (auto s : x1) {
        if (x1.size() > 1) {
            y1.push_back(-(2 * a1 * s + d1) / c1);
        } else {
            y1.push_back(sqrt((-4 * a1 * b1 * pow(s, 2) - 4 * b1 * d1 * s - 4 * b1 * f1 + pow(e1, 2)) / (4 * pow(b1, 2))) - e1 / (2 * b1));
            y1.push_back(-sqrt((-4 * a1 * b1 * pow(s, 2) - 4 * b1 * d1 * s - 4 * b1 * f1 + pow(e1, 2)) / (4 * pow(b1, 2))) - e1 / (2 * b1));
        }
    }
    for (int i = 0; i != y1.size(); i++) {
        if (x1.size() > 1) {
            cout << "(" << x1.at(i) << ", " << y1.at(i) << ")" << endl;
        } else {
            cout << "(" << x1.at(0) << ", " << y1.at(i) << ")" << endl;
        }
    }
    system("pause");
}
