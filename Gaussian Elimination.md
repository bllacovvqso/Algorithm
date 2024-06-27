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

## 0. Gaussian Elimination
 선형 시스템의 해를 구할 때는 가우스 소거법이 유용하다. 연립일차방정식을 행사다리꼴 행렬(REF, Row Echelon Form)로 변환하여 나타내는 것인데, 아래 예시를 참고하자.
 
 위와 같이 $Ax=B$라는 선형시스템을 계수행렬 $A$와 상수항 벡터 $b$로 나누고 이를 결합시킨 첨가행렬, Augmented Matrix를 얻을 수 있다. 이때 Elementary Row Operation이라는 연산을 통해 Augmented Matrix를 REF로 조작할 수 있다.
 
 
 그림과 같이 주대각선 좌측 하단의 성분이 모두 0으로 변환된 것이 REF이며 이렇게 정리돼야 해가 쉽게 구해진다. 
 
 
 ## 1. Consistency of Roots
 한편, 선형시스템의 해의 존재성은 pivot을 통해 알 수 있다. pivot은 REF의 각 행에서 처음으로 나타나는 0이 아닌 성분을 지칭한다. 

 RREF(Reduced REF)
