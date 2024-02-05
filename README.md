# 각종 최적화
##### 1. fastio

##### 2. #include <iomanip>와 long double
#include <iomanip>와 long double을 이용하면

    cout << fixed << setprecision(100) << p;
  
꼴로 작성하여 고정소숫점 표기법(fixed)으로 소숫점 아래 100째자리까지(set precision) 출력할 수 있다.

##### 3. Gaussian Elimination
연립방정식을 풀기 위해 문자의 계수를 행렬로 모으고 답을 포함하는 REF(Row Echelon Form matrix)로 변환하는 풀이이다.

    #include <iostream>
    #include <vector>
    using namespace std;
    using matrix=vector<vector<double>>;
    void print_matrix(matrix &a){
        for(int i=0;i<a.size();i++,puts(""))for(int q=0;q<a[i].size();q++)cout<<a[i][q]<<' ';
    }
    namespace GaussianElimination{
        int n; // n is the size of matrix M
        matrix M; // M is system of linear equations
        void run(){
            /**********Gaussian_Elimination***********/
            for(int i=0;i<n-1;i++){
                if(!M[i][i])for(int q=i+1;q<n;q++)if(M[q][i]){
                        for(int cur=-1;cur++<n;)swap(M[i][cur],M[q][cur]);
                        i--;break;
                    }
                else{
                    double v=M[i][i];
                    for(int q=-1;q++<n;)M[i][q]/=v;
                    for(int j=i;++j<n;){
                        double multiple=M[j][i];
                        for(int q=-1;q++<n;)M[j][q]-=multiple*M[i][q];
                    }
                }
            }
            double multiple=M[n-1][n-1];
            for(int q=-1;q++<n;)M[n-1][q]/=multiple;
            /**************Guass_Jordan***************/
            for(int i=n;--i>0;)for(int j=i;j--;){
                double multiple=M[j][i];
                for(int q=i-1;q++<n;)M[j][q]-=multiple*M[i][q];
            }
        }
    }
    
    int main(){
        using namespace GaussianElimination;
        cout<<"Enter the number of various >> ";
        cin>>n;
        M=vector<vector<double>>(n,vector<double>(n+1));
        for(int i=0;i<n;i++){
            cout<<"Enter the coefficients of the "<<i<<"th equation\n>>";
            for(int q=0;q<=n;q++)cin>>M[i][q];
        }
        cout<<"Matrix before run Gaussian Elimination\n";
        print_matrix(M);
        run();
        cout<<"Matrix after run Gaussian Elimination\n";
        print_matrix(M);
    }


### BOJ21296

x_N = n + d 로 지표와 가수로 나타내고 주어진 식을 간단하게 바꾸자.
d, d가 0.5 이상이면 1-d인 값이 1.5*998244353/(998244352.98*1e14) 이하여야함을 알 수 있다.
한편 x_N은 곧 N!을 나열했을 때 앞 두자리 수를 지표로 취하고 뒤 수들을 모두 가수로 취한 수이다.
따라서 가수가 매우 작아야 하므로, 문제는 곧 N!=abcdefg ... 로 나타날 때
ab.cdefg... 가 정수에 최대한 가까운 수를 찾는 것이다.

이러한 수는 Log를 취했을 때 가수부가 log0.1, log0.2, log0.3, ..., log9.9가 되므로
문제는 log(N!)의 가수부가 log0.1, log0.2, log0.3, ..., log9.9와 가장 가까운 수 N을 찾는 것이다.

##### 1. 고려사항
1. N 범위가 7부터 1e14이므로 eGPU 가속이 필요하다.
2. fastio 최적화





### BOJ21298
