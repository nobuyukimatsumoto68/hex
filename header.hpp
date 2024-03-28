#pragma once

#include <random>
#include <stack>
#include <cassert>

#include <sstream>
#include <fstream>

#include <string>
#include <functional>


using Idx = long;
constexpr int TWO = 2;
constexpr int THREE = 3;
constexpr int SIX = 6;


Idx Lx = 3;
Idx Ly = 3;
int nu = 1;
double K = 1.0;
int nparallel = 8;


std::mt19937 gen; // mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<double> d01D(0.0, 1.0); // (1, 6);
std::uniform_int_distribution<int> d01I(0, 1); // (1, 6);
std::uniform_int_distribution<Idx> d0N(0, Lx*Ly); // (1, 6);
void set_gen( const int seed ) { gen.seed(seed); }
double dist01(){ return d01D(gen); }
Idx dist0N(){ return d0N(gen); }
int distpm1(){ return 2*d01I(gen)-1; }


Idx mod(const Idx a, const int b){ return (b +(a%b))%b; }

Idx idx(const Idx x, const Idx y)  { return x + Lx*y; }

void get_xy(Idx& x, Idx& y, const Idx i)  {
  x = i%Lx;
  y = (i-x)/Lx;
}


bool is_site(const Idx x, const Idx y)  {
  const Idx c = mod(x-y, 3);
  bool res = true;
  if(c==1) res = false; // e.g., (1,0)
  return res;
}

bool is_site(const Idx i)  {
  Idx x,y;
  get_xy(x, y, i);
  return is_site(x,y);
}


bool is_link(const Idx x, const Idx y, const int mu)  {
  const Idx c = mod(x-y, 3);
  bool res = false;
  if(c==0 && mu<3) res = true; // e.g., (0,0)
  else if(c==2 && mu>=3) res = true; // e.g., (0,1)
  return res;
}

bool is_link(const Idx i, const int mu)  {
  Idx x,y;
  get_xy(x, y, i);
  return is_link(x,y,mu);
}


int cshift(Idx& xp, Idx& yp, const Idx x, const Idx y, const int mu)  {
  int res = 1;

  if(mu==0){
    xp=mod(x-1,Lx);
    yp=y;
    if(x==0 && nu>=3) res *= -1;
  }
  else if(mu==1){
    xp=mod(x+1,Lx);
    yp=mod(y-1,Ly);
    if(x==Lx-1 && nu>=3) res *= -1;
    if(y==0 && nu/2==1) res *= -1;
  }
  else if(mu==2){
    xp=x;
    yp=mod(y+1,Ly);
    if(y==Ly-1 && nu/2==1) res *= -1;
  }
  else if(mu==3){
    xp=mod(x+1,Lx);
    yp=y;
    if(x==Lx-1 && nu>=3) res *= -1;
  }
  else if(mu==4){
    xp=mod(x-1,Lx);
    yp=mod(y+1,Ly);
    if(x==0 && nu>=3) res *= -1;
    if(y==Ly-1 && nu/2==1) res *= -1;
  }
  else if(mu==5){
    xp=x;
    yp=mod(y-1,Ly);
    if(y==0 && nu/2==1 ) res *= -1;
  }
  else assert(false);
  return res;
}

int cshift(Idx& ip, const Idx i, const int mu)  {
  Idx x, y, xp, yp;
  get_xy(x, y, i);
  const int res = cshift(xp, yp, x, y, mu);
  ip = idx(xp, yp);
  return res;
}



struct Corr{
}




struct Spin {
  Idx N;
  std::vector<int> s;

  int& operator()(const Idx x, const Idx y) { return s[idx(x,y)]; }
  int operator()(const Idx x, const Idx y) const { return s[idx(x,y)]; }

  int& operator[](const Idx i) { return s[i]; }
  int operator[](const Idx i) const { return s[i]; }

  Spin( const int N_ )
    : N(N_)
    , s( N )
  {
    set1();
  }

  void set1() {
    for(Idx i=0; i<Lx*Ly; i++){
      if( is_site(i) ) s[i] = 1;
      else s[i] = 0;
    }
  }

  void random() {
    for(Idx i=0; i<Lx*Ly; i++){
      if( is_site(i) ) s[i] = distpm1();
      else s[i] = 0;
    }
  }

  int ss( const Idx x, const Idx y, const int mu ) const {
    Idx xp, yp;
    const int sign = cshift(xp, yp, x, y, mu);
    return sign * (*this)(x,y) * (*this)(xp,yp);
  }

  double ss_even( const int mu ) const {
    double tot = 0;
    Idx counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_link(x,y,mu) ) continue;
        tot += ss(x,y,mu);
        counter++;
      }
    }
    return tot/counter;
  }


  double ss_corr( const Idx dx, const Idx dy ) const {
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    double res = 0.0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;

        int sign = 1;
        if( x+dx>=Lx && nu>=3 ) sign *= -1;
        if( y+dy>=Ly && nu/2==1 ) sign *= -1;

        res += sign * (*this)(x,y) * (*this)(x+dx,y+dy);
      }}

    return res;
  }


  void ss_corr( std::vector<double>& corr ) const {
    // void ss_corr( Corr& corr ) const {
    assert( corr.size()==N );

#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        corr[idx(dx,dy)] = ss_corr( dx, dy );
      }}
  }


  std::string print() const {
    std::stringstream ss;
    for(int y=Ly-1; y>=0; y--){
      for(int x= 0; x<Lx; x++) {
        ss << std::setw(5) << (*this)(x, y);
      }
      ss << std::endl;
    }
    return ss.str();
  }
};


void heatbath( Spin& s ){
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
  for(Idx i=0; i<Lx*Ly; i++){
    if( !is_site(i) ) continue;

    int senv = 0;
    for(int mu=0; mu<SIX; mu++){
      if( !is_link(i,mu) ) continue;
      Idx j;
      const int sign = cshift( j, i, mu );
      senv += sign * s[j];
    }

    const double p = std::exp(2.0*K*senv);
    const double r = dist01();
    if( r<p/(1.0+p) ) s[i] = 1;
    else s[i] = -1;
  }
}


void wolff( Spin& s ){
  std::vector<bool> is_cluster(Lx*Ly, false);
  std::stack<Idx> stack_idx;

  Idx init = dist0N();
  while( !is_site(init) ) init = dist0N();

  is_cluster[init] = true;
  stack_idx.push(init);

  while( stack_idx.size() != 0 ){

    const Idx p = stack_idx.top();
    stack_idx.pop();
    s[p] = -s[p]; // flip when visited

    for(int mu = 0; mu < SIX; mu++){
      if( !is_link(p,mu) ) continue;
      Idx q;
      const int sign_q = cshift(q, p, mu);
      if( sign_q*s[q] == s[p] || is_cluster[q] ) continue; // s[x]*sR[y]<0 or y in c

      const double r = dist01();
      if( r < exp(-2.0 * K) ) continue; // reject

      is_cluster[q] = true;
      stack_idx.push(q);
    }
  }
}


template<typename T>
struct Obs {
  bool is_initialized = false;
  std::string description;
  int N;
  std::function<T(const Spin&)> f;

  T sum;
  int counter;

  Obs(){}

  Obs
  (
   const std::string& description_,
   const int N_,
   const std::function<T(const Spin&)>& f_
   )
    : is_initialized(true)
    , description(description_)
    , N(N_)
    , f(f_)
    , sum(0.0)
    , counter(0)
  {}

  void init
  (
   const std::string description_,
   const int N_,
   const std::function<T(const Spin&)>& f_
   ){
    is_initialized = true;
    description = description_;
    N = N_;
    f = f_;
    sum = 0.0;
    counter = 0;
  }

  void clear(){
    assert(is_initialized);
    sum = 0.0;
    counter = 0;
  }

  void meas( const Spin& s ) {
    assert(is_initialized);
    sum += f(s);
    counter++;
  }

  double mean() const {
    assert( is_initialized );
    assert( counter==N );
    return sum/N;
  }

};



void write( Obs<double>& obs, const std::string& dir, const int label ){
  const std::string filename = dir + obs.description + "_" + std::to_string(label) + ".dat";
  std::ofstream of( filename, std::ios::out | std::ios::trunc );
  if(!of) assert(false);

  of << std::scientific << std::setprecision(15);
  of << obs.mean();

  obs.clear();
}






// int cshift_minus(Idx& xp, Idx& yp, const Idx x, const Idx y, const int mu)  {
//   int res = 1;

//   if(mu==0){
//     xp=mod(x+1,Lx);
//     yp=y;
//     if(x==Lx-1 && nu>=3) res *= -1;
//   }
//   else if(mu==1){
//     xp=mod(x-1,Lx);
//     yp=mod(y+1,Ly);
//     if(x==0 && nu>=3) res *= -1;
//     if(y==Ly-1 && nu/2==1) res *= -1;
//   }
//   else if(mu==2){
//     xp=x;
//     yp=mod(y-1,Ly);
//     if(y==0 && nu/2==1) res *= -1;
//   }
//   else if(mu==3){
//     xp=mod(x-1,Lx);
//     yp=y;
//     if(x==0 && nu>=3) res *= -1;
//   }
//   else if(mu==4){
//     xp=mod(x+1,Lx);
//     yp=mod(y-1,Ly);
//     if(x==Lx-1 && nu>=3) res *= -1;
//     if(y==0 && nu/2==1) res *= -1;
//   }
//   else if(mu==5){
//     xp=x;
//     yp=mod(y+1,Ly);
//     if(y==Ly-1 && nu/2==1) res *= -1;
//   }
//   else assert(false);
//   return res;
// }

// int cshift_minus(Idx& ip, const Idx i, const int mu)  {
//   Idx x, y, xp, yp;
//   get_xy(x, y, i);
//   const int res = cshift_minus(xp, yp, x, y, mu);
//   ip = idx(xp,yp);
//   return res;
// }
