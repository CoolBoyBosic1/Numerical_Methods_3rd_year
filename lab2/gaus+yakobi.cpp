#include <bits/stdc++.h> 
using namespace std; 

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

class SLAE
{ 
private:
    Matrix A;
    Vector B;
    int n;

    double vector_norm_Linf(const Vector& v){
        double max_val = 0.0;
        for (int i = 0; i < n; i++) {
            max_val = max(max_val, fabs(v[i]));
        }
        return max_val;
    }

    double matrix_norm_Linf(const Matrix& M){
        double max_sum = 0.0;
        for (int i = 0; i < n; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < n; j++) {
                row_sum += fabs(M[i][j]);
            }
            max_sum = max(max_sum, row_sum);
        }
        return max_sum;
    }

    Matrix get_inverse_matrix(){
        Matrix Aug = A;
        Matrix Inv(n, Vector(n, 0.0));

        for(int i=0; i<n; ++i) {
            Inv[i][i] = 1.0;
        }

        for (int i = 0; i < n; i++) {
            int pivot_row = i;
            for (int j = i + 1; j < n; j++) {
                if (fabs(Aug[j][i]) > fabs(Aug[pivot_row][i])) {
                    pivot_row = j;
                }
            }

            if (pivot_row != i) {
                swap(Aug[i], Aug[pivot_row]);
                swap(Inv[i], Inv[pivot_row]);
            }

            if (fabs(Aug[i][i]) < 1e-12) {
                cout << "помилка: Матриця сингулярна, оберненої не існує!!!" << endl;
                return Matrix(n, Vector(n, 0.0));
            }

            double div = Aug[i][i];
            for (int j = 0; j < n; j++) {
                Aug[i][j] /= div;
                Inv[i][j] /= div;
            }

            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double factor = Aug[j][i];
                    for (int k = 0; k < n; k++) {
                        Aug[j][k] -= factor * Aug[i][k];
                        Inv[j][k] -= factor * Inv[i][k];
                    }
                }
            }
        }
        return Inv;
    }

public:
    SLAE(const Matrix& matrix, const Vector& vec){
        A = matrix;
        B = vec;
        n = B.size();
    }

    Vector gausian_method(){
        Matrix Aug = A;
        Vector b_temp = B;

        for(int i=0; i<n; ++i) {
            Aug[i].push_back(b_temp[i]);
        }

        for (int i = 0; i < n; i++) {
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (fabs(Aug[j][i]) > fabs(Aug[pivot][i])) {
                    pivot = j;
                }
            }
            swap(Aug[i], Aug[pivot]);

            double div = Aug[i][i];
            if (fabs(div) < 1e-12) {
                 cout << "помилка: Система має або безліч, або 0 розв'язків." << endl;
                 return Vector(n, 0.0);
            }
            
            for (int j = i; j <= n; j++) { 
                Aug[i][j] /= div;
            }

            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double factor = Aug[j][i];
                    for (int k = i; k <= n; k++) {
                        Aug[j][k] -= factor * Aug[i][k];
                    }
                }
            }
        }

        Vector X(n);
        for(int i=0; i<n; ++i) {
            X[i] = Aug[i][n];
        }
        return X;
    }

    void calculate_determinant(){
        Matrix tempA = A;
        double det = 1.0;
        int swaps = 0;

        for (int i = 0; i < n; i++) {
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (fabs(tempA[j][i]) > fabs(tempA[pivot][i])) {
                    pivot = j;
                }
            }

            if (pivot != i) {
                swap(tempA[i], tempA[pivot]);
                swaps++;
            }

            if (fabs(tempA[i][i]) < 1e-12) {
                det = 0.0;
                break;
            }

            det *= tempA[i][i]; 

            for (int j = i + 1; j < n; j++) {
                double factor = tempA[j][i] / tempA[i][i];
                for (int k = i; k < n; k++) {
                    tempA[j][k] -= factor * tempA[i][k];
                }
            }
        }
        
        if (swaps % 2 != 0) {
            det = -det;
        }
        
        cout << "Визначник: " << det << endl;
    }

    void calculate_condition_number(){
        cout << "Число обумовленості (L-inf норма)" << endl;
        
        double normA = matrix_norm_Linf(A);
        Matrix A_inv = get_inverse_matrix();
        double normA_inv = matrix_norm_Linf(A_inv);
        
        cout << "Норма A: " << normA << endl;
        cout << "Норма A_inv: " << normA_inv << endl;
        
        if(normA_inv == 0.0){
            cout << "Число обумовленості: (Матриця сингулярна)" << endl;
        } else {
            cout << "Число обумовленості: " << normA * normA_inv << endl;
        }
    }

    bool check_jacobi_convergence(){
        for (int i = 0; i < n; i++) {
            double diag = fabs(A[i][i]);
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += fabs(A[i][j]);
                }
            }
            if (diag <= sum) {
                cout << "Достатня умова збіжності Якобі НЕ виконана (рядок " << i << ")" << endl;
                return false;
            }
        }
        cout << "Достатня умова збіжності Якобі виконана." << endl;
        return true; 
    }

    Vector Yakobi(double epsilon = 0.0001){
        cout << "Метод Якобі" << endl;
        cout << "(Для умови припинення використовується L-inf норма)" << endl;

        bool converges = check_jacobi_convergence();
        if(!converges){
             cout << "УВАГА: Метод може не збігатися." << endl;
        }

        Vector X(n, 0.0); 
        Vector X_new(n, 0.0); 
        double norm = epsilon + 1.0;
        int iter = 0;

        while (norm > epsilon) {
            if (iter > 1000) {
               cout << "Перевищено ліміт ітерацій (1000)!" << endl;
               break;
            }
            
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        sum += A[i][j] * X[j];
                    }
                }
                X_new[i] = (B[i] - sum) / A[i][i];
            }

            Vector diff(n);
            for(int i=0; i<n; ++i) {
                diff[i] = X_new[i] - X[i];
            }
            
            norm = vector_norm_Linf(diff);
            
            if (iter == 0) {
                cout << "Ітерація 1: Норма = " << norm << endl;
            }

            X = X_new;
            iter++;
        }
        
        cout << "Остання ітерація (" << iter << "): Норма = " << norm << endl;
        return X;
    }

    void print_vector(const Vector& v, string name){
         cout << name << ": [ ";
         for(auto val : v) {
             cout << val << " ";
         }
         cout << "]" << endl;
    }
}; 

int main()
{ 
    
    int n;
    cout << "Введіть розмірність системи (n): ";
    cin >> n;

    Matrix A(n, Vector(n));
    Vector B(n);

    cout << "Введіть матрицю A (" << n << "x" << n << ") по рядках:" << endl;
    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            cin >> A[i][j];
        }
    }

    cout << "Введіть вектор B (" << n << "x1):" << endl;
    for(int i=0; i<n; ++i){
        cin >> B[i];
    }
    
    cout << "\n Система прийнята" << endl;
    
    cout << fixed << setprecision(6);

    SLAE system(A, B);

    cout << "\n1.Метод Гауса-Жордана" << endl;
    Vector X_gauss = system.gausian_method();
    system.print_vector(X_gauss, "Розв'язок (Гаус)");
    
    cout << "\n2.Визначник" << endl;
    system.calculate_determinant();
    
    cout << "\n3.Число обумовленості" << endl;
    system.calculate_condition_number();
    
    cout << "\n4.Метод Якобі" << endl;
    Vector X_jacobi = system.Yakobi(0.0001);
    system.print_vector(X_jacobi, "Розв'язок (Якобі)");

    return 0;
}
