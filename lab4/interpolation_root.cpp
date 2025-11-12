#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double f(double x) {
    return x * x * x - 4 * x * x - 4 * x + 13 - sin(x);
}

double direct_interpolation_root(double x0, double y0, double x1, double y1) {
    return x0 - y0 * (x1 - x0) / (y1 - y0);
}

double inverse_lagrange_root(const std::vector<double>& x, const std::vector<double>& y, int n) {
    double root_approximation = 0.0;

    for (int k = 0; k < n; ++k) {
        double l_k_at_zero = 1.0;
        for (int j = 0; j < n; ++j) {
            if (k == j) continue;
            l_k_at_zero *= y[j] / (y[j] - y[k]);
        }
        root_approximation += x[k] * l_k_at_zero;
    }
    return root_approximation;
}

int main() {
    const int n = 10;
    const double a = -2.0;
    const double b = -1.0;

    std::vector<double> x_nodes(n);
    std::vector<double> y_values(n);

    std::cout << std::fixed << std::setprecision(8);

    std::cout << "--- 1. Обчислення вузлiв Чебишова та значень f(x) ---" << std::endl;
    std::cout << "Iнтервал [a, b] = [" << a << ", " << b << "], n = " << n << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "| k |     x_k (Вузол)     |     y_k = f(x_k) (Значення) |" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    for (int k = 0; k < n; ++k) {
        double t_k = cos((2.0 * k + 1.0) * M_PI / (2.0 * n));
        x_nodes[k] = (b - a) / 2.0 * t_k + (a + b) / 2.0;
        y_values[k] = f(x_nodes[k]);

        std::cout << "| " << std::setw(1) << k << " | "
                  << std::setw(19) << x_nodes[k] << " | "
                  << std::setw(27) << y_values[k] << " |" << std::endl;
    }
    std::cout << "--------------------------------------------------------" << std::endl;

    std::cout << "\n--- 2. Метод 1: Пряма iнтерполяцiя (Метод хорд) ---" << std::endl;
    
    int k_left = -1, k_right = -1;
    for (int i = 0; i < n - 1; ++i) {
        if (y_values[i] * y_values[i + 1] < 0) {
            k_left = i;
            k_right = i + 1;
            break;
        }
    }

    if (k_left != -1) {
        double x0 = x_nodes[k_left];
        double y0 = y_values[k_left];
        double x1 = x_nodes[k_right];
        double y1 = y_values[k_right];

        std::cout << "Знайдено iнтервал зi змiною знаку мiж вузлами k=" << k_left << " та k=" << k_right << ":" << std::endl;
        std::cout << "  (x" << k_left << ", y" << k_left << ") = (" << x0 << ", " << y0 << ")" << std::endl;
        std::cout << "  (x" << k_right << ", y" << k_right << ") = (" << x1 << ", " << y1 << ")" << std::endl;

        double root_direct = direct_interpolation_root(x0, y0, x1, y1);
        std::cout << "\nНаближений корiнь (пряма iнтерполяцiя): " << root_direct << std::endl;
    } else {
        std::cout << "Не вдалося знайти сусiднi вузли зi змiною знаку." << std::endl;
    }

    std::cout << "\n--- 3. Метод 2: Обернена iнтерполяцiя (Полiном Лагранжа) ---" << std::endl;
    std::cout << "Використовуються всi " << n << " вузлiв для побудови G(y) та обчислення G(0)." << std::endl;
    
    double root_inverse = inverse_lagrange_root(x_nodes, y_values, n);
    
    std::cout << "\nНаближений корiнь (обернена iнтерполяцiя): " << root_inverse << std::endl;

    return 0;
}
