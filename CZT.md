# 1. FFT
```cpp
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll n;
const double PI=acos(-1);
typedef complex<double> cpx;
void FFT(vector<cpx> &v,cpx w){
    ll m=v.size();
    if(m==1)return;
    vector<cpx> odd(m/2),even(m/2);
    for(ll i=0;i<m;i++)
        if(i%2)odd[i/2]=v[i];
        else even[i/2]=v[i];
    FFT(even,w*w);FFT(odd,w*w);cpx z(1,0);
    for(ll i=0;i<m/2;i++,z*=w){
        v[i]=even[i]+z*odd[i];
        v[i+m/2]=even[i]-z*odd[i];
    }
}
vector<ll> multiply(vector<ll> v,vector<ll> u){
    vector<ll> w(n,0);
    ll m=1;
    while(m<=v.size()||m<=u.size())m<<=1;
    m<<=1;
    v.resize(m);u.resize(m);
    vector<cpx> _v(m),_u(m);
    for(ll i=0;i<m;i++){
        _v[i]=cpx(v[i],0);
        _u[i]=cpx(u[i],0);
    }
    cpx unit(cos(2*PI/m),sin(2*PI/m));
    FFT(_v,unit);
    FFT(_u,unit);
    vector<cpx> _w(m);
    for(ll i=0;i<m;i++)_w[i]=_v[i]*_u[i];
    FFT(_w,cpx(1,0)/unit);
    for(ll i=0;i<m;i++)_w[i]/=cpx(m,0);
    for(ll i=0;i<m;i++)w[i%n]+=round(_w[i].real());
    return w;
}
int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    cin>>n;
    vector<ll> a(n,0),w,c(n,0),b(n,0);
    for(ll i=1;i<=n/2;i++){
        if(2*i!=n){
            a[i*i%n]+=2;
            c[(2*i*i)%n]+=2;
        }else{
            a[i*i%n]++;
            c[(2*i*i)%n]++;
        }
        b[i]=i*i%n;
    }
    w=multiply(a,a);
    ll ans=0,sq=0;
    for(ll i=1;i<=n/2;i++){
        if(2*i!=n){ans+=w[b[i]]*2;sq+=c[b[i]]*2;}
        else{ans+=w[b[i]];sq+=c[b[i]];}
    }
    cout<<((ans+sq+1)/2);
}
```
# 2. FFT in recursive
```cpp
void FFT(vector<lint> &v, lint w){
    lint m = v.size();
    if(m == 1) return;
    
    vector<lint> odd(m/2), even(m/2);
    
    for(lint i=0;i<m;i++)
        if(i%2) odd[i/2] = v[i];
        else even[i/2] = v[i];
        
    FFT(even, w*w); FFT(odd, w*w);
    lint z(1,0);
    for(lint i = 0; i < m/2; i++, z*=w){
        v[i] = even[i] + z*odd[i];
        v[i+m/2] = even[i] - z*odd[i];
    }
}
```
# 3. Convolution with mod, FFT with mod.
```cpp
namespace fft{  
    
    __m256d mult(__m256d a, __m256d b) {
        __m256d c = _mm256_movedup_pd(a);
        __m256d d = _mm256_shuffle_pd(a, a, 15);
        __m256d cb = _mm256_mul_pd(c, b);
        __m256d db = _mm256_mul_pd(d, b);
        __m256d e = _mm256_shuffle_pd(db, db, 5);
        __m256d r = _mm256_addsub_pd(cb, e);
        return r;
    }
    void fft(int n, __m128d a[], bool invert) {
        for (int i = 1, j = 0; i < n; ++i) {
            int bit = n >> 1;
            for (; j >= bit; bit >>= 1) j -= bit;
            j += bit;
            if (i < j) swap(a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1) {
            double ang = 2 * 3.14159265358979 / len * (invert ? -1 : 1);
            __m256d wlen; wlen[0] = cos(ang), wlen[1] = sin(ang);
            for (int i = 0; i < n; i += len) {
                __m256d w; w[0] = 1; w[1] = 0;
                for (int j = 0; j < len / 2; ++j) {
                    w = _mm256_permute2f128_pd(w, w, 0);
                    wlen = _mm256_insertf128_pd(wlen, a[i + j + len / 2], 1);
                    w = mult(w, wlen);
                    __m128d vw = _mm256_extractf128_pd(w, 1);
                    __m128d u = a[i + j];
                    a[i + j] = _mm_add_pd(u, vw);
                    a[i + j + len / 2] = _mm_sub_pd(u, vw);
                }
            }
        }
        if (invert) {
            __m128d inv; inv[0] = inv[1] = 1.0 / n;
            for (int i = 0; i < n; ++i) a[i] = _mm_mul_pd(a[i], inv);
        }
    }

    vector<li> multiply_fft(const vector<li>& v,const vector<li>& w) {
        int n = 2; while (n < v.size() + w.size()) n <<= 1;
        __m128d* fv = new __m128d[n];
        for (int i = 0; i < n; ++i) fv[i][0] = fv[i][1] = 0;
        for (int i = 0; i < v.size(); ++i) fv[i][0] = v[i];
        for (int i = 0; i < w.size(); ++i) fv[i][1] = w[i];
        fft(n, fv, 0); //(a+bi) is stored in FFT
        for (int i = 0; i < n; i += 2) {
            __m256d a;
            a = _mm256_insertf128_pd(a, fv[i], 0);
            a = _mm256_insertf128_pd(a, fv[i + 1], 1);
            a = mult(a, a);
            fv[i] = _mm256_extractf128_pd(a, 0);
            fv[i + 1] = _mm256_extractf128_pd(a, 1);
        }
        fft(n, fv, 1);
        vector<li> ret(n);
        for (int i = 0; i < n; ++i) ret[i] = ((li)round(fv[i][1] / 2))%M;
        delete[] fv;
        return ret;
    }
    vector<li> mod_multiply_fft(const vector<li>& v,const vector<li>& w){
        vector<li> A1(v.size()),A2(v),B1(w.size()),B2(w);
        for(int i=0;i<v.size();i++){
            A1[i]=A2[i]&(131071);
            A2[i]>>=17;
        } // note. A1 + 2^17 * A2
        for(int i=0;i<w.size();i++){
            B1[i]=B2[i]&(131071);
            B2[i]>>=17;
        } // note. B1 + 2^17 * B2
        vector<li> A1B1 = multiply_fft(A1,B1);
        vector<li> A1B2 = multiply_fft(A1,B2);
        vector<li> A2B1 = multiply_fft(A2,B1);
        vector<li> A2B2 = multiply_fft(A2,B2);
        int _sz = v.size()+w.size()-1;
        A1B1.resize(_sz);A1B2.resize(_sz);A2B1.resize(_sz);A2B2.resize(_sz);
        for(int i=0;i<_sz;i++)
            A1B1[i] = (A1B1[i]+((A1B2[i]+A2B1[i])<<17)%M+(A2B2[i]*((1LL<<34)%M))%M)%M;
        return A1B1;
    }
}
```
Another version:
```cpp
using i64 = int64_t;
using u64 = uint64_t;
using li = int64_t;
using base = complex<double>;
i64 N,P;

inline i64 add(i64 a,i64 b){return a+b>=P?a+b-P:a+b;}
inline i64 sub(i64 a,i64 b){return a<b?a-b+P:a-b;}
inline i64 mul(i64 a,i64 b){return a*b%P;}
i64 ipow(i64 a,i64 b){
    i64 res=1;
    for(a%=P;b;b>>=1,a=a*a%P)if(b&1)res=res*a%P;
    return res;
}
inline i64 modInv(i64 a){return ipow(a,P-2);}

void fft(vector<base> &a, bool inv){
   int n = a.size(), j = 0;
   for(int i=1; i<n; i++){
      int bit = (n >> 1);
      while(j >= bit){
         j -= bit;
         bit >>= 1;
      }
      j += bit;
      if(i < j) swap(a[i], a[j]);
   }
   double ang = 2 * acos(-1) / n * (inv ? -1 : 1);
   
   vector<base> roots(n/2);

   for(int i=0; i<n/2; i++)
      roots[i] = base(cos(ang * i), sin(ang * i));

   for(int i=2; i<=n; i<<=1){
      int step = n / i;
      for(int j=0; j<n; j+=i){
         for(int k=0; k<i/2; k++){
            base u = a[j+k], v = a[j+k+i/2] * roots[step * k];
            a[j+k] = u+v;
            a[j+k+i/2] = u-v;
         }
      }
   }
   if(inv) for(int i=0; i<n; i++) a[i] /= n; 
}

vector<i64> multiply(vector<i64> &v, vector<i64> &w){
   int n = 2; while(n < v.size() + w.size()) n <<= 1;
   vector<base> v1(n), v2(n), r1(n), r2(n);
   for(int i=0; i<v.size(); i++)
      v1[i] = base(v[i] >> 15, v[i] & 32767);
   for(int i=0; i<w.size(); i++)
      v2[i] = base(w[i] >> 15, w[i] & 32767);
   fft(v1, 0);
   fft(v2, 0);
   for(int i=0; i<n; i++){
      int j = (i ? (n - i) : i);
      base ans1 = (v1[i] + conj(v1[j])) * base(0.5, 0);
      base ans2 = (v1[i] - conj(v1[j])) * base(0, -0.5);
      base ans3 = (v2[i] + conj(v2[j])) * base(0.5, 0);
      base ans4 = (v2[i] - conj(v2[j])) * base(0, -0.5);
      r1[i] = (ans1 * ans3) + (ans1 * ans4) * base(0, 1);
      r2[i] = (ans2 * ans3) + (ans2 * ans4) * base(0, 1);
   }
   fft(r1, 1);
   fft(r2, 1);
   vector<i64> ret(n);
   for(int i=0; i<n; i++){
      i64 av = ((i64)round(r1[i].real())%P+P)%P;
      i64 bv = ((i64)round(r1[i].imag())%P + (i64)round(r2[i].real())%P + 2*P)%P;
      i64 cv = ((i64)round(r2[i].imag())%P+P)%P;
      ret[i] = add(add(mul(av,(1<<30)%P),mul(bv,1<<15)),cv);
   }
   return ret;
}
```
