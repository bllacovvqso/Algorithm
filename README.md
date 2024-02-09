# Gaussian Elimination
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
1. N 범위가면

### BOJ30641
거듭제곱을 여러번 더할 때 그 시간 복잡도를 줄이기 위해선 따로 거듭제곱을 계산하는 함수를 호출하기보다 DP를 해야 한다. 이때 굳이 배열은 필요 없고 아래와 같이 하면 된다.
 ##### 시간 초과 코드
     #include <iostream>
    #define lint long long int

    using namespace std;

    lint sum;
    lint mod = 1000000007;

    lint powmod(lint n)  
    {
        lint res = 1;
    	if (n == 0) return 1;
        for (lint i = 1; i <= n; i++)
        {
            res *= 26; 
            res %= mod;
        }
        
        return res;
    }


    lint slvSum(lint L, lint U)
    {	
    	lint lTmp = (L-1)/2;
    	lint uTmp = (U-1)/2;
    	
    	if(L == U) return powmod(lTmp);
    	else if (lTmp == uTmp) return (2*powmod(lTmp))%mod;
    	
    	if (L%2) sum += 2*powmod(lTmp); 
    	else sum += powmod(lTmp);
    	
    	if (U%2) sum += powmod(uTmp);
    	else sum += 2*powmod(uTmp);
    	
    	if (uTmp-lTmp >= 2)
    	{
    		for (int i = lTmp+1; i <= uTmp-1; i++)
    		{
    			sum += powmod(i);
    			sum %= mod;
    		}
    		sum *= 2;	
    	}
    	
    	return sum%mod;
    }

    int main()
    {
        lint L, U;
        cin >> L >> U;
        
        if (L == 2 && U == 2) cout << "H" << '\n'; 
        else if (L == 1 && U == 1) cout << "H" << '\n';
        else if (L == 2) cout << "H" << '\n';
    	else cout << "A" << '\n';
    	
        cout << slvSum(L, U);
    }



    

이걸 수정해주면 다음과 같다.
    #include <iostream>
    #define lint long long int

    using namespace std;

    lint sum;
    lint mod = 1000000007;



    lint slvSum(lint L, lint U)
    {	
    	lint multi = 1;

    	for (int i = 1; i <= U; ++i)
        {
        	if (i <= 2 && i >= L) sum += 1;
        	else if (i <= 2) continue;
        	else
        	{
        		if (i % 2 == 1) multi = (multi * 26) % mod;
          		if (i >= L) sum = (sum + multi) % mod;	
    		}
        }

    	return sum;
    }

    int main()
    {
        lint L, U;
        cin >> L >> U;

        if (L == 2 && U == 2) cout << "H" << '\n'; 
        else if (L == 1 && U == 1) cout << "H" << '\n';
        else if (L == 2) cout << "H" << '\n';
    	else cout << "A" << '\n';

        cout << slvSum(L, U);
    }



#std:map
BOJ11444와 같은 경우를 보면 dp도 항상 배열이 아닌 적당한 자료구조에 해야 함을 알 수 있다.

    #include <iostream>
    #define lint unsigned long int
    const int max_size = 10000000;
    using namespace std;
    
    int fibo(lint n);
    
    int mod = 1000000007;
    lint n, tmp1, tmp2, tmp3;
    int dp[max_size];
    
    int main()
    {
        dp[0] = 0;
        dp[1] = 1;
        dp[2] = 1;
        for (int i = 1; i <= 10000; i++) {
    
        cin >> n;
        cout << fibo(n);
        }
    }
    int fibo(lint n)
    {
      if (n < max_size && dp[n]!=0) return dp[n]; 
      if (n == 0) return 0;
      if (n%2)
      {
        tmp1 = fibo((n+1)/2)*fibo((n+1)/2);
        tmp2 = fibo((n-1)/2)*fibo((n-1)/2);
        tmp3=(tmp1%mod+tmp2%mod)%mod;
        if(n < max_size)
        {
            dp[n]= tmp3;
            return dp[n];
        } 
        return tmp3;
      }
      else
      {
        tmp1 = 2*fibo(n/2-1)+fibo(n/2);
        tmp2 = (fibo(n/2)*(tmp1%mod))%mod;
        if(n < max_size)
        {
            dp[n] = tmp2;
            return dp[n];
        }
        return tmp2; 
      }
    }

이 문제 특성 상 사용하는 배열에 버려지는 칸이 있으면 TL이 되는데, 사용한 점화식은 f(n/2)를 저장하기 때문에 배열에 버려지는 칸이 너무 많아진다. 따라서 map을 사용해야 한다.
