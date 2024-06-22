# 1. I/O Optimization
- 표준입출력(getchar, puts, printf, etc)이 cin/cout보다 빠르다.
-  ₩n이 endl보다 빠르다.
- 아래 코드로 C++ 표준 스트림과 C 표준 스트림의 입출력 연산 후 동기화 여부를 설정한다. 
  sync를 false로 두면 사용 버퍼 수가 줄어들어 가속되지만 싱글 쓰레드 환경에서 코드를 작성해야 한다. 이렇게 해도 표준입출력이 더 빠르다.
```cpp
ios_base::sync_with_stdio(bool sync);
```
실전에서 아래처럼 사용할 수 있는데 비추.
```cpp
#define fastio ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
main(){
fastio;
}
```
# 2. 

```cpp
#define _CRT_SECURE_NO_WARNINGS
#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#pragma pack(1);
```
