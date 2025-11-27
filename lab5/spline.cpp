#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;


struct SplineTuple {
    double a, b, c, d, x;
};

double func(double x) {
    return log10(log10(x));
}

// Побудова сплайна
vector<SplineTuple> buildSpline(const vector<double>& x, const vector<double>& y, int n) {
    vector<SplineTuple> splines(n);
    
    // Ініціалізація a_i = y_i
    for (int i = 0; i < n; ++i) {
        splines[i].x = x[i];
        splines[i].a = y[i];
    }
    
    // Крок 1: Обчислення h_i
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }
    
    // Крок 2: Прямий хід методу прогонки (для знаходження c_i)
    // Система лінійних рівнянь для c_i має тридіагональну матрицю
    vector<double> alpha(n - 1);
    vector<double> beta(n - 1);
    
    // Граничні умови для природного сплайна: c_0 = 0, c_n = 0
    alpha[0] = 0.0;
    beta[0] = 0.0;
    
    for (int i = 1; i < n - 1; ++i) {
        double A = h[i - 1];
        double C = 2.0 * (h[i - 1] + h[i]);
        double B = h[i];
        double F = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        
        double z = (A * alpha[i - 1] + C);
        alpha[i] = -B / z;
        beta[i] = (F - A * beta[i - 1]) / z;
    }
    
    // Крок 3: Зворотний хід (знаходимо c_i, потім b_i та d_i)
    splines[n - 1].c = 0.0; // Гранична умова
    
    for (int i = n - 2; i >= 0; --i) {
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
    }
    
    for (int i = n - 2; i >= 0; --i) {
        double hi = h[i];
        splines[i].d = (splines[i + 1].c - splines[i].c) / (3.0 * hi);
        splines[i].b = (y[i + 1] - y[i]) / hi - (hi * (splines[i + 1].c + 2.0 * splines[i].c)) / 3.0;
    }
    
    return splines;
}

int main() {
    cout << fixed << setprecision(5);

    double start = 1.2;
    double end = 4.0;   
    int n;

    cout << "Variant 4: f(x) = lg(lg(x))" << endl;
    cout << "Enter the number of interpolation nodes (>= 2): ";
    cin >> n;

    if (n < 2) {
        cout << "Error: Need at least 2 points." << endl;
        return 1;
    }

    // Генерація вузлів
    vector<double> x(n), y(n);
    double step = (end - start) / (n - 1);

    cout << "\n--- Interpolation Nodes ---" << endl;
    cout << "X\t\tY" << endl;
    for (int i = 0; i < n; ++i) {
        x[i] = start + i * step;
        y[i] = func(x[i]);
        cout << x[i] << "\t" << y[i] << endl;
    }

    // Побудова сплайна
    vector<SplineTuple> splines = buildSpline(x, y, n);

    // Вивід результатів у консоль (аналітичний вигляд поліномів)
    cout << "\n--- Cubic Spline Coefficients ---" << endl;
    cout << "Interval\t\ta_i\t\tb_i\t\tc_i\t\td_i" << endl;
    for (int i = 0; i < n - 1; ++i) {
        cout << "[" << splines[i].x << "; " << x[i+1] << "]\t"
             << splines[i].a << "\t" << splines[i].b << "\t" 
             << splines[i].c << "\t" << splines[i].d << endl;
    }

    // Генерація даних для графіку (файл data.txt)
    ofstream outFile("data.txt");
    outFile << "X Original Spline" << endl;
    
    double plotStep = 0.05;
    for (double xi = start; xi <= end; xi += plotStep) {
        SplineTuple s;
        if (xi >= x[n-1]) s = splines[n-2]; 
        else {
            for (int i = 0; i < n - 1; ++i) {
                if (xi >= x[i] && xi < x[i+1]) {
                    s = splines[i];
                    break;
                }
            }
        }
        
        double delta = xi - s.x;
        double splineVal = s.a + s.b * delta + s.c * pow(delta, 2) + s.d * pow(delta, 3);
        double originalVal = func(xi);
        
        outFile << xi << " " << originalVal << " " << splineVal << endl;
    }
    outFile.close();
    
    cout << "\nData for plotting saved to 'data.txt'." << endl;

    return 0;
}
