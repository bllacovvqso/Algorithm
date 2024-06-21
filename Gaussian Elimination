연립방정식을 풀기 위해 문자의 계수를 행렬로 모으고 답을 포함하는 REF(Row Echelon Form matrix)로 변환하는 풀이.
```cpp
#include <iostream>
#include <vector>
using namespace std;

using matrix = vector<vector<double>>;

void print_matrix(matrix &a) {
    for (int i = 0; i < a.size(); i++, puts(""))
        for (int q = 0; q < a[i].size(); q++) 
            cout << a[i][q] << ' ';
}

namespace GaussianElimination {
    int n; // n is the size of matrix M
    matrix M; // M is system of linear equations

    void run() {
        /**********Gaussian_Elimination***********/
        for (int i = 0; i < n - 1; i++) {
            if (!M[i][i])
                for (int q = i + 1; q < n; q++)
                    if (M[q][i]) {
                        for (int cur = -1; cur++ < n; ) 
                            swap(M[i][cur], M[q][cur]);
                        i--;
                        break;
                    }
            else {
                double v = M[i][i];
                for (int q = -1; q++ < n; ) 
                    M[i][q] /= v;
                for (int j = i; ++j < n; ) {
                    double multiple = M[j][i];
                    for (int q = -1; q++ < n; ) 
                        M[j][q] -= multiple * M[i][q];
                }
            }
        }
        double multiple = M[n - 1][n - 1];
        for (int q = -1; q++ < n; ) 
            M[n - 1][q] /= multiple;
        
        /**************Guass_Jordan***************/
        for (int i = n; --i > 0; )
            for (int j = i; j--; ) {
                double multiple = M[j][i];
                for (int q = i - 1; q++ < n; ) 
                    M[j][q] -= multiple * M[i][q];
            }
    }
}

int main() {
    using namespace GaussianElimination;
    cout << "Enter the number of variables >> ";
    cin >> n;
    M = vector<vector<double>>(n, vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        cout << "Enter the coefficients of the " << i << "th equation\n>> ";
        for (int q = 0; q <= n; q++) 
            cin >> M[i][q];
    }
    cout << "Matrix before running Gaussian Elimination\n";
    print_matrix(M);
    run();
    cout << "Matrix after running Gaussian Elimination\n";
    print_matrix(M);
}

```
