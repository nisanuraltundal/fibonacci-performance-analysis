#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

typedef vector<vector<long long>> MatrixF;
MatrixF multiplyF(MatrixF A, MatrixF B) {
    MatrixF C(2, vector<long long>(2));
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];
    return C;
}

MatrixF powerF(MatrixF A, int p) {
    MatrixF res = { {1, 0}, {0, 1} };
    while (p > 0) {
        if (p & 1) res = multiplyF(res, A);
        A = multiplyF(A, A);
        p >>= 1;
    }
    return res;
}

long long fibonacciMatrix(int n) {
    if (n <= 1) return n;
    MatrixF T = { {1, 1}, {1, 0} };
    T = powerF (T, n - 1);
    return T[0][0];
}

long long fibonacciIterative(int n) {
    if (n <= 1) return n;
    long long prev = 0, curr = 1;
    for (int i = 2; i <= n; i++) {
        long long next = prev + curr;
        prev = curr;
        curr = next;
    }
    return curr;
}

int fibonacciRecursive(int n) {
    if (n <= 1) return n;
    return fibonacciRecursive(n - 1) + fibonacciRecursive(n - 2);
}


// Standard (Naive) Matrix Multiplication
void standardMultiply(const vector<vector<double> >& A, const vector<vector<double> >& B,vector<vector<double>>& C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < N; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

// Tiled (Blocked) Matrix Multiplication
void tiledMultiply(const vector<vector<double> >& A, const vector<vector<double> >& B,vector<vector<double>>& C, int N, int BLOCK_SIZE) {
    for (int i = 0; i < N; i += BLOCK_SIZE) {
        for (int j = 0; j < N; j += BLOCK_SIZE) {
            for (int k = 0; k < N; k += BLOCK_SIZE) {
                // Computing sub-blocks
                for (int ii = i; ii < i + BLOCK_SIZE && ii < N; ii++) {
                    for (int jj = j; jj < j + BLOCK_SIZE && jj < N; jj++) {
                        double sum = C[ii][jj];
                        for (int kk = k; kk < k + BLOCK_SIZE && kk < N; kk++) {
                            sum += A[ii][kk] * B[kk][jj];
                        }
                        C[ii][jj] = sum;
                    }
                }
            }
        }
    }
}


// ANA TEST MERKEZÝ 
int main() {
    int n_fib = 40;
    int n_matrix = 1000;
    int N1024 = 1024;
    int BSIZE = 32;
    string flag = "Od"; // for flags

    

    // IKJ MATRIX TEST
    vector<int> Aikj(n_matrix * n_matrix, 1), Bikj(n_matrix * n_matrix, 1), Cikj(n_matrix * n_matrix, 0);
    auto s_ikj= chrono::high_resolution_clock::now();
    for (int i = 0; i < n_matrix; i++) {
        for (int k = 0; k < n_matrix; k++) {
            int a_val = Aikj[i * n_matrix + k];
            for (int j = 0; j < n_matrix; j++) {
                Cikj[i * n_matrix + j] += a_val * Bikj[k * n_matrix + j];
            }
        }
    }
    auto e_ikj = chrono::high_resolution_clock::now();
    chrono::duration<double> d_ikj = e_ikj - s_ikj;

    // STANDART MATRÝX
    vector<vector<double>> A1024(N1024, vector<double>(N1024, 1.1)), B1024(N1024, vector<double>(N1024, 2.2)), C_std(N1024, vector<double>(N1024, 0.0));
    auto s_std1024 = chrono::high_resolution_clock::now();
    standardMultiply(A1024, B1024, C_std, N1024);
    auto e_std1024 = chrono::high_resolution_clock::now();
    chrono::duration<double> d_std1024 = e_std1024 - s_std1024;

    // TILED MATRÝS
    vector<vector<double>> C_tile(N1024, vector<double>(N1024, 0.0));
    auto s_tile1024 = chrono::high_resolution_clock::now();
    tiledMultiply(A1024, B1024, C_tile, N1024, BSIZE);
    auto e_tile1024 = chrono::high_resolution_clock::now();
    chrono::duration<double> d_tile1024 = e_tile1024 - s_tile1024;

    // FIBONACCI RECURSIVE
    auto s_fib_rec = chrono::high_resolution_clock::now();
    fibonacciRecursive(n_fib);
    auto e_fib_rec = chrono::high_resolution_clock::now();
    chrono::duration<double> d_fib_rec = e_fib_rec - s_fib_rec;

    // FIBONACCI ITERATIVE
    auto s_fib_itr = chrono::high_resolution_clock::now();
    fibonacciIterative(n_fib);
    auto e_fib_itr = chrono::high_resolution_clock::now();
    chrono::duration<double> d_fib_itr = e_fib_itr - s_fib_itr;

    // FIBONACCI MATRIX 
    auto s_fib_mat = chrono::high_resolution_clock::now();
    fibonacciMatrix(n_fib);
    auto e_fib_mat = chrono::high_resolution_clock::now();
    chrono::duration<double> d_fib_mat = e_fib_mat - s_fib_mat;

    // EXCEL KAYIT
    ofstream excel("deney_sonuclari.csv", ios::app);
    if (excel.is_open()) {
        excel << "IKJ_Matrix;" << n_matrix << ";" << flag << ";" << d_ikj.count() << endl;
        excel << "Matrix_Standart_1024;" << N1024 << ";" << flag << ";" << d_std1024.count() << endl;
        excel << "Matrix_Tiled_1024;" << N1024 << ";" << flag << ";" << d_tile1024.count() << endl;
        excel << "Fib_Recursive;" << n_fib << ";" << flag << ";" << d_fib_rec.count() << endl;
        excel << "Fib_Iterative;" << n_fib << ";" << flag << ";" << d_fib_itr.count() << endl;
        excel << "Fib_Matrix_Method;" << n_fib << ";" << flag << ";" << d_fib_mat.count() << endl;
        excel.close();
        cout << "The test is completed successfully! (" << flag << ")!" << endl;
    }
    return 0;
}