#define _USE_MATH_DEFINES
#include <algorithm>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

constexpr double ZERO = 1e-6;                   /* 判断变量是否等于0的精度阈值 */
constexpr double Y_ACCURACY_THRESHOLD = 1e-9;   /* 交点y轴方向的精度阈值 */
constexpr double SK_ACCURACY_THRESHOLD = 1e-9;  /* 判断ak和bk是否相等的精度阈值 */

/**
 * @brief 从控制台读取椭圆信息并获取四个极值点坐标(按顺时针排序)
 * @param[out] A 椭圆一般方程中的参数A
 * @param[out] B 椭圆一般方程中的参数B
 * @param[out] C 椭圆一般方程中的参数C
 * @param[out] D 椭圆一般方程中的参数D
 * @param[out] E 椭圆一般方程中的参数E
 * @param[out] F 椭圆一般方程中的参数F
 * @return vector<pair<double, double>> 椭圆的四个极值点坐标
 */
vector<pair<double, double>> getEllipse(double &A, double &B, double &C, double &D, double &E, double &F);
/**
 * @brief 获取两个椭圆的四个弧两两组合成的16个弧对
 * @param[in] extreme_points0 第一个椭圆的四个极值点
 * @param[in] extreme_points1 第二个椭圆的四个极值点
 * @return vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> 16个弧对
 */
vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> getArcsPairs(vector<pair<double, double>> extreme_points0, vector<pair<double, double>> extreme_points1);
/**
 * @brief 将椭圆分成四个单调弧
 * @param[in] extreme_points 椭圆的四个极值点
 * @return vector<pair<pair<double, double>, pair<double, double>>> 椭圆的四个单调弧
 */
vector<pair<pair<double, double>, pair<double, double>>> getArcs(vector<pair<double, double>> extreme_points);
/**
 * @brief 判断一个圆锥曲线的一般方程是否是椭圆
 * @param[in] A 圆锥曲线一般方程中的参数A
 * @param[in] B 圆锥曲线一般方程中的参数B
 * @param[in] C 圆锥曲线一般方程中的参数C
 * @param[in] D 圆锥曲线一般方程中的参数D
 * @param[in] E 圆锥曲线一般方程中的参数E
 * @param[in] F 圆锥曲线一般方程中的参数F
 * @return bool 圆锥曲线是否是椭圆
 *      -true   圆锥曲线是椭圆
 *      -false  圆锥曲线不是椭圆
 */
bool isEllipse(double A, double B, double C, double D, double E, double F);
/**
 * @brief 根据椭圆的一般方程计算四个极值点(未排序)
 * @param[in] A 椭圆一般方程中的参数A
 * @param[in] B 椭圆一般方程中的参数B
 * @param[in] C 椭圆一般方程中的参数C
 * @param[in] D 椭圆一般方程中的参数D
 * @param[in] E 椭圆一般方程中的参数E
 * @param[in] F 椭圆一般方程中的参数F
 * return vector<pair<double, double>> 椭圆的四个极值点坐标(未排序)
 */
vector<pair<double, double>> calcExtremePoint(double A, double B, double C, double D, double E, double F);
/**
 * @brief 解一元二次方程
 * @param[in] a 一元二次方程一般形式的参数a
 * @param[in] b 一元二次方程一般形式的参数b
 * @param[in] c 一元二次方程一般形式的参数c
 * @return vector<double> 一元二次方程的解
 */
vector<double> solveQuadraticEquation(double a, double b, double c);
/**
 * @brief 将椭圆的四个极值点按顺时针方向从最上方的点开始排序
 * @param[in] extreme_points 椭圆的四个极值点
 * @return vector<pair<double, double>> 排序后的四个极值点
 */
vector<pair<double, double>> sortExtremePoint(vector<pair<double, double>> extreme_points);
/**
 * @brief 判断一个弧对是否明显不相交, 即所在矩形是否重叠
 * @param[in] arcs_pairs 弧对
 * @return bool 弧对是否明显不相交
 *      -true   弧对明显不相交
 *      -false  弧对可能相交
 */
bool isObviouslyNotIntersect(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pairs);
/**
 * @brief 计算一组弧对交点
 * @param[in] arcs_pair 弧对
 * @param[out] intersections 结果集
 * @param[in] A0 第一个椭圆的一般方程的参数A
 * @param[in] B0 第一个椭圆的一般方程的参数B
 * @param[in] C0 第一个椭圆的一般方程的参数C
 * @param[in] D0 第一个椭圆的一般方程的参数D
 * @param[in] E0 第一个椭圆的一般方程的参数E
 * @param[in] F0 第一个椭圆的一般方程的参数F
 * @param[in] extreme_points0 第一个椭圆的四个极值点
 * @param[in] A1 第二个椭圆的一般方程的参数A
 * @param[in] B1 第二个椭圆的一般方程的参数B
 * @param[in] C1 第二个椭圆的一般方程的参数C
 * @param[in] D1 第二个椭圆的一般方程的参数D
 * @param[in] E1 第二个椭圆的一般方程的参数E
 * @param[in] F1 第二个椭圆的一般方程的参数F
 * @param[in] extreme_points1 第二个椭圆的四个极值点
 */
void calcIntersection(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double A0, double B0, double C0, double D0, double E0, double F0,vector<pair<double, double>> extreme_points0, double A1, double B1, double C1, double D1, double E1, double F1, vector<pair<double, double>> extreme_points1);
/**
 * @brief 计算一段椭圆弧在某一横坐标处对应的纵坐标
 * @param[in] x 横坐标
 * @param[in] arc 椭圆弧
 * @param[in] A 椭圆一般方程的参数A
 * @param[in] B 椭圆一般方程的参数B
 * @param[in] C 椭圆一般方程的参数C
 * @param[in] D 椭圆一般方程的参数D
 * @param[in] E 椭圆一般方程的参数E
 * @param[in] F 椭圆一般方程的参数F
 * @return double 计算得到的纵坐标
 */
double getFx(double x, pair<pair<double, double>, pair<double, double>> arc, double A, double B, double C, double D, double E, double F);
/**
 * @brief 计算一个椭圆弧编号(按顺时针方向从右上方的椭圆弧开始排序)
 * @param[in] 椭圆弧
 * @return int 编号
 */
int getArcNo(pair<pair<double, double>, pair<double, double>> arc);
/**
 * @brief 计算一段椭圆弧在某一点的切矢与x轴的夹角
 * @param[in] x 需计算点的横坐标
 * @param[in] y 需计算点的纵坐标
 * @param[in] left 椭圆左极点横坐标
 * @param[in] arc 椭圆弧
 * @param[in] A 椭圆一般方程的参数A
 * @param[in] B 椭圆一般方程的参数B
 * @param[in] C 椭圆一般方程的参数C
 * @param[in] D 椭圆一般方程的参数D
 * @param[in] E 椭圆一般方程的参数E
 * @param[in] F 椭圆一般方程的参数F
 * return double 计算得到的夹角
 */
double getAngle(double x, double y, double left, pair<pair<double, double>, pair<double, double>> arc, double A, double B, double C, double D, double E, double F);
/**
 * @brief 计算一段椭圆弧所在的最小水平矩形的右上和左下角的坐标
 * @param[in] arc 椭圆弧
 * @return vector<double> 计算得到的坐标(按右上纵坐标、右上横坐标、左下纵坐标、左下横坐标排列)
 */
vector<double> getRectangle(pair<pair<double, double>, pair<double, double>> arc);
/**
 * @brief 将某一点放入计算结果集中
 * @param[out] intersections 结果集
 * @param[in] point 点坐标
 */
void addIntersection(vector<pair<double, double>> &intersections, pair<double, double> point);
/**
 * @brief 计算当 X < 0 时两弧的交点坐标
 * @param[in] arcs_pair 弧对
 * @param[out] intersections 结果集
 * @param[in] x0 两椭圆弧横轴交集的最左侧横坐标
 * @param[in] x1 两椭圆弧横轴交集的最右侧横坐标
 * @param[in] delta0 x0处的delta值
 * @param[in] delta1 x1处的delta值
 * @param[in] A0 第一个椭圆一般方程的参数A
 * @param[in] B0 第一个椭圆一般方程的参数B
 * @param[in] C0 第一个椭圆一般方程的参数C
 * @param[in] D0 第一个椭圆一般方程的参数D
 * @param[in] E0 第一个椭圆一般方程的参数E
 * @param[in] F0 第一个椭圆一般方程的参数F
 * @param[in] A1 第二个椭圆一般方程的参数A
 * @param[in] B1 第二个椭圆一般方程的参数B
 * @param[in] C1 第二个椭圆一般方程的参数C
 * @param[in] D1 第二个椭圆一般方程的参数D
 * @param[in] E1 第二个椭圆一般方程的参数E
 * @param[in] F1 第二个椭圆一般方程的参数F
 */
void XLessThan0(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double x0, double x1, double delta0, double delta1, double A0, double B0, double C0, double D0, double E0, double F0, double A1, double B1, double C1, double D1, double E1, double F1);
/**
 * @brief 计算两椭圆弧在切矢相同位置的横坐标
 * @param[in] arcs_pair 弧对
 * @param[in] x0 两椭圆弧横轴交集的最左侧横坐标
 * @param[in] x1 两椭圆弧横轴交集的最右侧横坐标
 * @param[in] a0 第一个椭圆弧在x0处的切矢方向
 * @param[in] b0 第二个椭圆弧在x0处的切矢方向
 * @param[in] A0 第一个椭圆一般方程的参数A
 * @param[in] B0 第一个椭圆一般方程的参数B
 * @param[in] C0 第一个椭圆一般方程的参数C
 * @param[in] D0 第一个椭圆一般方程的参数D
 * @param[in] E0 第一个椭圆一般方程的参数E
 * @param[in] F0 第一个椭圆一般方程的参数F
 * @param[in] extreme_points0 第一个椭圆的极值点坐标
 * @param[in] A1 第二个椭圆一般方程的参数A
 * @param[in] B1 第二个椭圆一般方程的参数B
 * @param[in] C1 第二个椭圆一般方程的参数C
 * @param[in] D1 第二个椭圆一般方程的参数D
 * @param[in] E1 第二个椭圆一般方程的参数E
 * @param[in] F1 第二个椭圆一般方程的参数F
 * @param[in] extreme_points1 第二个椭圆的极值点坐标
 */
double xkWhenSkEqual0(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, double x0, double x1, double a0, double b0, double A0, double B0, double C0, double D0, double E0, double F0, vector<pair<double, double>> extreme_points0, double A1, double B1, double C1, double D1, double E1, double F1, vector<pair<double, double>> extreme_points1);

int main()
{
    double A1, B1, C1, D1, E1, F1, A2, B2, C2, D2, E2, F2;  /* 两个椭圆一般方程的参数 */
    vector<pair<double, double>> extreme_points0, extreme_points1, intersections;   /* 两个椭圆的极值点坐标以及交点坐标 */
    vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> arcs_pairs;    /* 两个椭圆的四个单调弧两两组合构成的16个弧对 */
    /* 读取椭圆参数并分别计算两个椭圆的极值点坐标 */
    extreme_points0 = getEllipse(A1, B1, C1, D1, E1, F1);
    extreme_points1 = getEllipse(A2, B2, C2, D2, E2, F2);
    /* 获取两个椭圆构成的16个弧对 */
    arcs_pairs = getArcsPairs(extreme_points0, extreme_points1);
    /* 计算交点, 每轮循环计算一个弧对 */
    for (auto s : arcs_pairs) {
        /* 如果两椭圆弧之间显然不相交, 即两椭圆弧所在最小水平矩形不相交, 则直接进入下一轮循环 */
        if (isObviouslyNotIntersect(s)) {
            continue;
        }
        /* 计算两椭圆弧交点, 并将结果记录到intersections中 */
        calcIntersection(s, intersections, A1, B1, C1, D1, E1, F1, extreme_points0, A2, B2, C2, D2, E2, F2, extreme_points1);
    }
    /* 输出结果 */
    if (intersections.empty()) {
        cout << "两椭圆没有交点\n";
    } else {
        cout << "两椭圆的交点坐标为\n";
        for (auto s : intersections) {
            cout << "(" << (abs(s.first) <= ZERO ? 0 : s.first) << ", " << (abs(s.second) <= ZERO ? 0 : s.second) << ")\n";
        }
    }
    system("pause");
}

vector<pair<double, double>> getEllipse(double &A, double &B, double &C, double &D, double &E, double &F)
{
    static int cont = 1;        /* 当前输入的椭圆序号 */
    double h, k, a, b, theta;   /* 椭圆圆心横坐标, 椭圆圆心纵坐标, 长轴长, 短轴长, 长轴与x轴夹角 */
    bool need_retry = false;    /* 是否需要重新输入 */
    vector<pair<double, double>> extreme_points;    /* 椭圆的几点坐标 */
    cout << "请输入第" << cont << "个椭圆的椭圆心横坐标、椭圆心纵坐标、长轴长、短轴长、长轴与x轴夹角，以空格为分隔\n";
    while (true) {
        cin >> h >> k >> a >> b >> theta;
        /* 在每次读取完成后清空缓冲区, 避免读取到上一次输入的多余参数 */
        cin.sync();
        /* 如果输入非法内容, 例如输入了字母, 则清除cin的异常标志并要求重新输入 */
        if (cin.rdstate() != 0) {
            cout << "输入异常，请重新输入\n";
            cin.clear();
            continue;
        }
        /* 如果输入的参数不符合要求, 则进行相应提示并要求重新输入 */
        if (a < 0) {
            cout << "输入的长轴长小于0\n";
            need_retry = true;
        }
        if (b < 0) {
            cout << "输入的短轴长小于0\n";
            need_retry = true;
        }
        if (a < b) {
            cout << "输入的长轴长小于短轴\n";
            need_retry = true;
        }
        if (theta < 0 || theta >= 180) {
            cout << "输入的长轴与x轴夹角不在[0, 180)的范围内\n";
            need_retry = true;
        }
        if (need_retry) {
            cout << "请重新输入\n";
            need_retry = false;
            continue;
        }
        /* 根据输入的参数计算椭圆的一般方程 */
        a /= 2;
        b /= 2;
        theta = theta * M_PI / 180;
        A = (cos(theta) * cos(theta)) / (a * a) + (sin(theta) * sin(theta)) / (b * b);
        B = 2 * cos(theta) * sin(theta) * (1 / (a * a) - 1 / (b * b));
        C = (sin(theta) * sin(theta)) / (a * a) + (cos(theta) * cos(theta)) / (b * b);
        D = -(2 * A * h + B * k);
        E = -(2 * C * k + B * h);
        F = A * h * h + B * h * k + C * k * k - 1;
        swap(B, C);
        /* 如果计算得到的一般方程无法构成椭圆, 则进行提示并要求重新输入 */
        if (!isEllipse(A, B, C, D, E, F)) {
            cout << "参数无法构成椭圆，请重新输入\n";
            continue;
        }
        /* 输入完成后跳出输入循环 */
        break;
    }
    /* 计算椭圆的极值点坐标 */
    extreme_points = calcExtremePoint(A, B, C, D, E, F);
    extreme_points = sortExtremePoint(extreme_points);
    cont++;
    return extreme_points;
}

vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> getArcsPairs(vector<pair<double, double>> extreme_points0, vector<pair<double, double>> extreme_points1)
{
    vector<pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>>> arcs_pairs;    /* 弧对的集合 */
    vector<pair<pair<double, double>, pair<double, double>>> arcs0, arcs1;  /* 两个椭圆各自的四个单调弧 */
    /* 分别获取两个椭圆的四个单调弧 */
    arcs0 = getArcs(extreme_points0);
    arcs1 = getArcs(extreme_points1);
    /* 将两个椭圆的四个单调弧两两组合成16个弧对 */
    for (auto iter0 = arcs0.begin(); iter0 != arcs0.end(); iter0++) {
        for (auto iter1 = arcs1.begin(); iter1 != arcs1.end(); iter1++) {
            arcs_pairs.push_back(make_pair(*iter0, *iter1));
        }
    }
    return arcs_pairs;
}

vector<pair<pair<double, double>, pair<double, double>>> getArcs(vector<pair<double, double>> extreme_points)
{
    vector<pair<pair<double, double>, pair<double, double>>> arcs;  /* 弧的集合 */
    /* 按顺时针顺序, 从右上角的弧开始将弧放入弧的集合中 */
    for (auto iter = extreme_points.begin(); iter != extreme_points.end() - 1; iter++) {
        arcs.push_back(make_pair(*iter, *(iter + 1)));
    }
    arcs.push_back(make_pair(extreme_points.back(), extreme_points.at(0)));
    return arcs;
}

bool isEllipse(double A, double B, double C, double D, double E, double F)
{
    double P = D / 2, Q = E / 2, I1, I2, I3;    /* 圆锥曲线判别式需要用到的参数 */
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
    vector<double> solution = solveQuadraticEquation(4 * A * A * B - A * C * C, 4 * A * B * D - 2 * A * C * E, B * D * D - C * D * E + C * C * F);  /* 通过数学计算得到的极值点的横坐标 */
    vector<pair<double, double>> extreme_points;    /* 极值点坐标 */
    /* 对于上下极值点, 存在横坐标相同和不相同两种情况, 分别计算纵坐标 */
    if (solution.size() == 2) {
        extreme_points.push_back(make_pair(solution.at(0), -(2 * A * solution.at(0) + D) / C));
        extreme_points.push_back(make_pair(solution.at(1), -(2 * A * solution.at(1) + D) / C));
    } else if (solution.size() == 1) {
        extreme_points.push_back(make_pair(solution.at(0), sqrt((-4 * A * B * solution.at(0) * solution.at(0) - 4 * B * D * solution.at(0) - 4 * B * F + E * E) / (4 * B * B)) - E / (2 * B)));
        extreme_points.push_back(make_pair(solution.at(0), -sqrt((-4 * A * B * solution.at(0) * solution.at(0) - 4 * B * D * solution.at(0) - 4 * B * F + E * E) / (4 * B * B)) - E / (2 * B)));
    }
    /* 左右极值点的横坐标必然不同 */
    solution = solveQuadraticEquation(4 * A * B - C * C, 4 * B * D - 2 * C * E, 4 * B * F - E * E);
    extreme_points.push_back(make_pair(solution.at(0), -(C * solution.at(0) + E) / (2 * B)));
    extreme_points.push_back(make_pair(solution.at(1), -(C * solution.at(1) + E) / (2 * B)));
    return extreme_points;
}

vector<double> solveQuadraticEquation(double a, double b, double c)
{
    vector<double> result;              /* 一元二次方程的解 */
    double delta = b * b - 4 * a * c;   /* 一元二次方程的delta */
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
    /* 将纵坐标最大的极值点(上极值点)放到第一位 */
    iter_swap(extreme_points.begin(), max_element(extreme_points.begin(), extreme_points.end(), [](pair<double, double> a, pair<double, double> b) { return a.second < b.second; }));
    /* 将横坐标最大的极值点(右极值点)放到第二位 */
    iter_swap(extreme_points.begin() + 1, max_element(extreme_points.begin() + 1, extreme_points.end(), [](pair<double, double> a, pair<double, double> b) { return a.first < b.first; }));
    /* 将纵坐标最小的极值点(下极值点)放到第三位 */
    iter_swap(extreme_points.begin() + 2, min_element(extreme_points.begin() + 2, extreme_points.end(), [](pair<double, double> a, pair<double, double> b) { return a.second < b.second; }));
    return extreme_points;
}

bool isObviouslyNotIntersect(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pairs)
{
    vector<double> rectangle0, rectangle1;  /* 两个椭圆所在的矩形区域, 按照上边界纵坐标、右边界横坐标、下边界纵坐标、左边界横坐标顺序排列 */
    rectangle0 = getRectangle(arcs_pairs.first);
    rectangle1 = getRectangle(arcs_pairs.second);
    /* 根据两个矩形的边界坐标判断两个矩形是否相交 */
    if (rectangle0.at(0) < rectangle1.at(2) || rectangle0.at(1) < rectangle1.at(3) || rectangle0.at(2) > rectangle1.at(0) || rectangle0.at(3) > rectangle1.at(1)) {
        return true;
    } else {
        return false;
    }
}

void calcIntersection(pair<pair<pair<double, double>, pair<double, double>>, pair<pair<double, double>, pair<double, double>>> arcs_pair, vector<pair<double, double>> &intersections, double A0, double B0, double C0, double D0, double E0, double F0, vector<pair<double, double>> extreme_points0, double A1, double B1, double C1, double D1, double E1, double F1, vector<pair<double, double>> extreme_points1)
{
    /**
     * x0       两椭圆弧横轴交集的最左侧横坐标
     * x1       两椭圆弧横轴交集的最右侧横坐标
     * FP0      第一个椭圆弧在x0处的纵坐标
     * FP1      第一个椭圆弧在x1处的纵坐标
     * FQ0      第二个椭圆弧在x0处的纵坐标
     * FQ1      第二个椭圆弧在x1处的纵坐标
     * delta0   FP0 - FQ1
     * delta1   FP1 - FQ1
     * X        delta0 * delta1
     * xk       两椭圆弧切矢方向相同的位置的横坐标
     * FPk      第一个椭圆弧在xk处的纵坐标
     * FQk      第二个椭圆弧在xk处的纵坐标
     * deltak   FPk - FQk
     * a0       第一个椭圆弧在x0处的切矢与x轴的夹角
     * a1       第一个椭圆弧在x1处的切矢与x轴的夹角
     * b0       第二个椭圆弧在x0处的切矢与x轴的夹角
     * b1       第二个椭圆弧在x1处的切矢与x轴的夹角
     */
    double x0, x1, FP0, FP1, FQ0, FQ1, delta0, delta1, X, xk, FPk, FQk, deltak, a0, a1, b0, b1;
    x0 = max(min(arcs_pair.first.first.first, arcs_pair.first.second.first), min(arcs_pair.second.first.first, arcs_pair.second.second.first));
    x1 = min(max(arcs_pair.first.first.first, arcs_pair.first.second.first), max(arcs_pair.second.first.first, arcs_pair.second.second.first));
    FP0 = getFx(x0, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQ0 = getFx(x0, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    FP1 = getFx(x1, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQ1 = getFx(x1, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    delta0 = FP0 - FQ0;
    delta1 = FP1 - FQ1;
    /* 若delta0或delta1等于0, 则两椭圆弧在端点处相交 */
    if (abs(delta0) <= ZERO) {
        addIntersection(intersections, make_pair(x0, FP0));
        return;
    } else if (abs(delta1) <= ZERO) {
        addIntersection(intersections, make_pair(x1, FP1));
        return;
    }
    X = delta0 * delta1;
    /* 若 X < 0, 则两椭圆弧有一个交点 */
    if (X < 0) {
        XLessThan0(arcs_pair, intersections, x0, x1, delta0, delta1, A0, B0, C0, D0, E0, F0, A1, B1, C1, D1, E1, F1);
        return;
    }
    a0 = getAngle(x0, FP0, extreme_points0.at(3).first, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    a1 = getAngle(x1, FP1, extreme_points0.at(3).first, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    b0 = getAngle(x0, FQ0, extreme_points1.at(3).first, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    b1 = getAngle(x1, FQ1, extreme_points1.at(3).first, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    /* 若 (a0 - b0) * (a1 - b1) > 0, 则两椭圆弧没有交点 */
    if ((a0 - b0) * (a1 - b1) > 0) {
        return;
    }
    xk = xkWhenSkEqual0(arcs_pair, x0, x1, a0, b0, A0, B0, C0, D0, E0, F0, extreme_points0, A1, B1, C1, D1, E1, F1, extreme_points1);
    FPk = getFx(xk, arcs_pair.first, A0, B0, C0, D0, E0, F0);
    FQk = getFx(xk, arcs_pair.second, A1, B1, C1, D1, E1, F1);
    deltak = FPk - FQk;
    /* 若deltak等于0, 则两椭圆弧在(xk, FPk)处相交 */
    if (abs(deltak) <= ZERO) {
        addIntersection(intersections, make_pair(xk, FPk));
        return;
    }
    /* 若 delta0 * deltak > 0, 则两椭圆弧没有交点 */
    if (delta0 * deltak > 0) {
        return;
    }
    /* 若 delta0 * deltak < 0, 则两椭圆弧在(x0, xk)和(xk, x0)区间内各有一个交点 */
    XLessThan0(arcs_pair, intersections, x0, xk, delta0, deltak, A0, B0, C0, D0, E0, F0, A1, B1, C1, D1, E1, F1);
    XLessThan0(arcs_pair, intersections, xk, x1, deltak, delta1, A0, B0, C0, D0, E0, F0, A1, B1, C1, D1, E1, F1);
}

double getFx(double x, pair<pair<double, double>, pair<double, double>> arc, double A, double B, double C, double D, double E, double F)
{
    int arc_no = getArcNo(arc); /* 椭圆弧的编号 */
    vector<double> solution;    /* 椭圆在横坐标x处的纵坐标 */
    solution = solveQuadraticEquation(B, C * x + E, A * x * x + D * x + F);
    /* 如果是第一或第四段弧, 就采用较大的解, 如果是第二或第三段弧, 就采用较小的解 */
    if (solution.size() == 2) {
        if (arc_no == 1 || arc_no == 4) {
            return solution.at(0);
        } else {
            return solution.at(1);
        }
    } else {    /* 如果只有一个解, 说明x处于椭圆的左右极点, 直接使用解作为结果 */
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
    if (find_if(intersections.begin(), intersections.end(), [point](pair<double, double> stored_point) { return abs(stored_point.first - point.first) <= ZERO && abs(stored_point.second - point.second) <= ZERO; }) == intersections.end()) {
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
