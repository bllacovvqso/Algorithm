# 각종 최적화
##### 1. fastio

##### 2. #include <iomanip>와 long double
#include <iomanip>와 long double을 이용하면

    cout << fixed << setprecision(100) << p;
  
꼴로 작성하여 고정소숫점 표기법(fixed)으로 소숫점 아래 100째자리까지(set precision) 출력할 수 있다.


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
