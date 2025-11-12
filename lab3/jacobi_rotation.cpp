#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

void printMatrix(const Matrix& A) {
    for (const auto& row : A) {
        for (double val : row) {
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << val;
        }
        std::cout << std::endl;
    }
}

void printVector(const Vector& v) {
    for (double val : v) {
        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << val;
    }
    std::cout << std::endl;
}

void findMaxOffDiagonal(const Matrix& A, int& p, int& q) {
    double maxVal = 0.0;
    int n = A.size();
    p = 0; // Ініціалізуємо, щоб уникнути невизначеної поведінки
    q = 1; // якщо матриця вже діагональна
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(A[i][j]) > maxVal) {
                maxVal = std::abs(A[i][j]);
                p = i;
                q = j;
            }
        }
    }
}

void jacobiRotation(Matrix& A, Matrix& V, double epsilon = 1e-10, int maxIterations = 1000) {
    int n = A.size();
    
    V.assign(n, Vector(n, 0.0));
    for (int i = 0; i < n; ++i) {
        V[i][i] = 1.0;
    }

    std::cout << "\n--- Ітерація 0 (Початковий стан) ---\n";
    std::cout << "Матриця A:\n";
    printMatrix(A);
    std::cout << "Матриця V:\n";
    printMatrix(V);

    for (int iter = 0; iter < maxIterations; ++iter) {
        int p, q;
        findMaxOffDiagonal(A, p, q);

        double apq = A[p][q];
        
        if (std::abs(apq) < epsilon) {
            std::cout << "\nАлгоритм збігся за " << iter << " ітерацій.\n";
            break;
        }

        double app = A[p][p];
        double aqq = A[q][q];
        
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0) ? 1.0 / (tau + std::sqrt(1.0 + tau * tau)) 
                            : -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
        double c = 1.0 / std::sqrt(1.0 + t * t);
        double s = t * c;

        Matrix A_old = A; // Зберігаємо стару матрицю A для розрахунків
        A[p][p] = app - t * apq;
        A[q][q] = aqq + t * apq;
        A[p][q] = 0.0;
        A[q][p] = 0.0;

        for (int i = 0; i < n; ++i) {
            if (i != p && i != q) {
                double aip = A_old[i][p];
                double aiq = A_old[i][q];
                A[i][p] = c * aip - s * aiq;
                A[p][i] = A[i][p];
                A[i][q] = s * aip + c * aiq;
                A[q][i] = A[i][q];
            }
        }
        
        Matrix V_old = V; // Зберігаємо стару матрицю V для розрахунків
        for (int i = 0; i < n; ++i) {
            double vip = V_old[i][p];
            double viq = V_old[i][q];
            V[i][p] = c * vip - s * viq;
            V[i][q] = s * vip + c * viq;
        }

        std::cout << "\n--- Ітерація " << iter + 1 << " (обертання p=" << p << ", q=" << q << ") ---\n";
        std::cout << "Матриця A:\n";
        printMatrix(A);
        std::cout << "Матриця V:\n";
        printMatrix(V);

        if (iter == maxIterations - 1) {
            std::cout << "\nДосягнуто ліміту ітерацій (" << maxIterations << ").\n";
        }
    }
}

int main() {
    Matrix A = {
        {2.3, 3.5, 1.4},
        {3.5, 0.4, 0.6},
        {1.4, 0.6, 1.3}
    };
    
    int n = A.size();
    Matrix V; 

    jacobiRotation(A, V);

    std::cout << "\n=============================================\n";
    std::cout << "           КІНЦЕВИЙ РЕЗУЛЬТАТ\n";
    std::cout << "=============================================\n";

    std::cout << "\nФінальна матриця A (власні значення на діагоналі):\n";
    printMatrix(A);

    Vector eigenvalues(n);
    for (int i = 0; i < n; ++i) {
        eigenvalues[i] = A[i][i];
    }

    std::cout << "\nВласні значення (lambda):\n";
    printVector(eigenvalues);

    std::cout << "\nФінальна матриця власних векторів V (стовпці - вектори):\n";
    printMatrix(V);

    return 0;
}
