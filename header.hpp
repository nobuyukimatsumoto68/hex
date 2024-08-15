#pragma once

#include <random>
#include <stack>
#include <cassert>

#include <sstream>
#include <fstream>

#include <string>
#include <vector>
#include <array>
#include <functional>
// #include <any>


using Idx = std::size_t;
// using Idx = unsigned long int;
// using Idx = int;
using Double = long double;
// using int = short;

using namespace std;

constexpr int TWO = 2;
constexpr int THREE = 3;
constexpr int SIX = 6;

const double abs_tautil = 1.2; // 1.0 <= abs
const double arg_tautil = 4.0/9.0 * M_PI; // pi/3.0 <= arg <= pi/2.0 // 3.35
// const double abs_tautil = 1.0; // 1.0 <= abs
// const double arg_tautil = 3.0/9.0 * M_PI; // pi/3.0 <= arg <= pi/2.0 // 3.35
const bool tautil_default = false; // true->below
double tautil1 = 0.3420201433256688;
double tautil2 = 0.9396926207859083 + 0.00001;


int mult = 8;
// int mult = 4;
Idx Lx = 3*1*mult; // 12
Idx Ly = 3*1*mult;

// constexpr Idx Lx = 6*4; // 12
// constexpr Idx Ly = 2*Lx;
constexpr int nparallel = 12; //12


constexpr int nu = 1; // PP, PA, AA, AP
// const Double beta_c = 0.5 * std::log(2.0 + std::sqrt(3.0));
//
// constexpr Double kappa = 2.0/3.0;
// const Double cos6 = std::cos(M_PI/6.0);
// const Double B = cos6 / (1.0 - kappa*kappa * cos6 * cos6);
//
// constexpr Double alat = 1.0/Lx;
// constexpr Double xipsi = std::sqrt( 1.5*std::sqrt(3.0)*alat / (2.0*M_PI) );

#ifndef _OPENMP
// int omp_get_thread_num() { return 0; }
#endif



std::mt19937 gen;
std::uniform_real_distribution<Double> d01D(0.0, 1.0); // (1, 6);
std::uniform_int_distribution<int> d01I(0, 1); // (1, 6);
std::uniform_int_distribution<Idx> d0N(0, Lx*Ly-1); // (1, 6);
void set_gen( const int seed ) {
  std::mt19937 gen0( seed );
  gen.seed( gen0() );
}
Double dist01(){ return d01D(gen); }
Idx dist0N(){ return d0N(gen); }
int distpm1(){ return 2*d01I(gen)-1; }


// ---------------
// GLOBAL FUNCTIONS
// ---------------

double ell0[2]; // = {1.0, 0.0};
double ell1[2]; // = -omega;
double ell2[2]; // -ell2 - ell0

double ell[3]; // {|ell0|, |ell1|, |ell2|}

double ell_star0[2];
double ell_star1[2];
double ell_star2[2];

double e0[2];
double e1[2];
double e2[2];


double kappa[3];
double theta[3];
double cosH[3];
double Beta[3];


void set_tautil(){
  if(!tautil_default){
    std::cout << "setting tautil from abs, arg" << std::endl;
    tautil1 = abs_tautil*std::cos(arg_tautil);
    tautil2 = abs_tautil*std::sin(arg_tautil);
  }
}


void set_ell(){
  ell1[0] = 1.0+tautil1;
  ell1[1] = tautil2;

  ell0[0] = 1.0-2.0*tautil1;
  ell0[1] = -2.0*tautil2;

  ell2[0] = -2.0+tautil1;
  ell2[1] = tautil2;

  ell[0] = std::sqrt( ell0[0]*ell0[0] + ell0[1]*ell0[1] );
  ell[1] = std::sqrt( ell1[0]*ell1[0] + ell1[1]*ell1[1] );
  ell[2] = std::sqrt( ell2[0]*ell2[0] + ell2[1]*ell2[1] );
}


void set_ell_star(){
  const double s = 0.5 * (ell[0]+ell[1]+ell[2]);
  const double area = std::sqrt( s*(s-ell[0])*(s-ell[1])*(s-ell[2]) );

  double coeff;
  coeff = 0.25 * (ell[1]*ell[1] + ell[2]*ell[2] - ell[0]*ell[0]) / area;
  ell_star0[0] =  ell0[1] * coeff;
  ell_star0[1] = -ell0[0] * coeff;

  coeff = 0.25 * (ell[2]*ell[2] + ell[0]*ell[0] - ell[1]*ell[1]) / area;
  ell_star1[0] =  ell1[1] * coeff;
  ell_star1[1] = -ell1[0] * coeff;

  coeff = 0.25 * (ell[0]*ell[0] + ell[1]*ell[1] - ell[2]*ell[2]) / area;
  ell_star2[0] =  ell2[1] * coeff;
  ell_star2[1] = -ell2[0] * coeff;
}


void set_kappa(){
  kappa[0] = 2.0*ell[0] / (ell[0] + ell[1] + ell[2]);
  kappa[1] = 2.0*ell[1] / (ell[0] + ell[1] + ell[2]);
  kappa[2] = 2.0*ell[2] / (ell[0] + ell[1] + ell[2]);
}


void set_e(){
  double len;

  len = std::sqrt( ell_star0[0]*ell_star0[0] + ell_star0[1]*ell_star0[1] );
  e0[0] = ell_star0[0]/len;
  e0[1] = ell_star0[1]/len;

  len = std::sqrt( ell_star1[0]*ell_star1[0] + ell_star1[1]*ell_star1[1] );
  e1[0] = ell_star1[0]/len;
  e1[1] = ell_star1[1]/len;

  len = std::sqrt( ell_star2[0]*ell_star2[0] + ell_star2[1]*ell_star2[1] );
  e2[0] = ell_star2[0]/len;
  e2[1] = ell_star2[1]/len;
}


void set_theta(){
  theta[0] = std::acos( -e1[0]*e2[0]-e1[1]*e2[1] );
  theta[1] = std::acos( -e2[0]*e0[0]-e2[1]*e0[1] );
  theta[2] = std::acos( -e0[0]*e1[0]-e0[1]*e1[1] );
}


void set_cos(){
  cosH[0] = std::cos(0.5*theta[0]);
  cosH[1] = std::cos(0.5*theta[1]);
  cosH[2] = std::cos(0.5*theta[2]);
}


void set_beta(){
  Beta[0] = std::atanh( kappa[0]*cosH[1]*cosH[2]/cosH[0] );
  Beta[1] = std::atanh( kappa[1]*cosH[2]*cosH[0]/cosH[1] );
  Beta[2] = std::atanh( kappa[2]*cosH[0]*cosH[1]/cosH[2] );
}


void set_all(){
  set_tautil();
  set_ell();
  set_ell_star();
  set_kappa();
  set_e();
  set_theta();
  set_cos();
  set_beta();
}





Idx idx(const Idx x, const Idx y)  { return x + Lx*y; }

void get_xy(Idx& x, Idx& y, const Idx i)  {
  x = (i+Lx)%Lx;
  y = (i-x)/Lx;
}

inline int get_char( const Idx x, const Idx y) { return (x-y+Lx*Ly)%3; }

bool is_site(const Idx x, const Idx y)  {
  const Idx c = get_char(x,y);
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
  const Idx c = get_char(x,y);
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


void cshift(Idx& xp, Idx& yp, const Idx x, const Idx y, const int mu)  {
  int dx = -(mu+2)%3+1;
  int dy = -(mu+1)%3+1;
  if(mu>=3){
    dx *= -1;
    dy *= -1;
  }
  xp = (x+dx+Lx)%Lx;
  yp = (y+dy+Ly)%Ly;
}

void cshift(Idx& ip, const Idx i, const int mu)  {
  Idx x, y, xp, yp;
  get_xy(x, y, i);
  cshift(xp, yp, x, y, mu);
  ip = idx(xp, yp);
}



// ----------------
// CLASS DEFINITIONS
// ----------------


struct Spin {
  Idx N;
  std::vector<int> s;

  inline int& operator()(const Idx x, const Idx y) { return s[idx(x,y)]; }
  inline int operator()(const Idx x, const Idx y) const { return s[idx(x,y)]; }

  int& operator[](const Idx i) { return s[i]; }
  int operator[](const Idx i) const { return s[i]; }

  Spin() = delete;

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
    cshift(xp, yp, x, y, mu);
    return (*this)(x,y) * (*this)(xp,yp);
  }

  Double ss_even( const int mu ) const {
    Double tot = 0.0;
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

  Double ss_corr( const Idx dx, const Idx dy ) const {
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        // if( !is_site(x,y) ) continue;
        if( (x-y+Lx*Ly)%3!=0 ) continue;
        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;
        if( !is_site(xp,yp) ) continue;

        res += (*this)(x,y) * (*this)(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> ss_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        corr[idx(dx,dy)] = ss_corr( dx, dy );
      }}

    return corr;
  }


//   Double eps( const Idx x, const Idx y ) const {
//     assert(0<=x && x<Lx);
//     assert(0<=y && y<Ly);
//     assert( is_site(x,y) );

//     Double res = 0.0;
//     for(int mu=0; mu<SIX; mu++){
//       if( !is_link(x,y,mu) ) continue;
//       Idx xp, yp;
//       cshift( xp, yp, x, y, mu );
//       res += (*this)(x,y) * (*this)(xp,yp);
//     }
//     // res *= 0.5 * kappa * B; // @@@

//     return res;
//   }


//   Double eps_1pt() const {
//     Double res = 0.0;

// // #ifdef _OPENMP
// // #pragma omp parallel for collapse(2) num_threads(nparallel)
// // #endif
//     int counter = 0;
//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         if( !is_site(x,y) ) continue;
//         res += eps( x, y );
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }


//   Double epseps_corr( const Idx dx, const Idx dy ) const {
//     assert(0<=dx && dx<Lx);
//     assert(0<=dy && dy<Ly);

//     Double res = 0.0;
//     int counter = 0;

//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         // if( !is_site(x,y) ) continue;
//         if( (x-y+Lx*Ly)%3!=0 ) continue;
//         const Idx xp = (x+dx)%Lx;
//         const Idx yp = (y+dy)%Ly;
//         if( !is_site(xp,yp) ) continue;

//         res += eps(x,y) * eps(xp,yp);
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }


//   std::vector<Double> epseps_corr() const {
//     std::vector<Double> corr(N, 0.0);

// #ifdef _OPENMP
// #pragma omp parallel for collapse(2) num_threads(nparallel)
//     // #pragma omp parallel for num_threads(nparallel) schedule(static)
// #endif
//     for(Idx dx=0; dx<Lx; dx++){
//       for(Idx dy=0; dy<Ly; dy++){
//         corr[idx(dx,dy)] = epseps_corr( dx, dy );
//       }}

//     return corr;
//   }


  Double K( const Idx x, const Idx y, const int mu ) const {
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);
    // assert(0<=mu && mu<3);
    assert( is_link(x,y,mu) );

    Idx xp, yp;
    cshift( xp, yp, x, y, mu );

    // const Double res = B * (*this)(x,y)*(*this)(xp,yp);
    const Double res = (*this)(x,y)*(*this)(xp,yp); // @@@
    return res;
  }


  Double K_1pt( const int mu ) const {
    // std::vector<Double> tmp(nparallel, 0.0);
    // std::vector<int> counter(nparallel, 0);
    int counter = 0;
    Double res = 0.0;

    // #ifdef _OPENMP
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
    // #endif
    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_link(x,y,mu) ) continue;
        res += K( x, y, mu );
        counter++;
        // tmp[omp_get_thread_num()] += K( x, y, mu );
        // counter[omp_get_thread_num()]++;
      }}

    // Double res = 0.0;
    // int tot = 0;
    // for(int i=0; i<nparallel; i++) {
    //   res += tmp[i];
    //   tot += counter[i];
    // }
    res /= counter;
    // res /= tot;
    return res;
  }

  Double KK_corr( const Idx dx, const Idx dy, const int mu, const int nu ) const {
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        // if( !is_site(x,y) ) continue;
        // if( (x-y+Lx*Ly)%3!=0 ) continue;
        if( !is_link(x,y,mu) ) continue;
        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;
        if( !is_link(xp,yp,nu) ) continue;
        // if( !is_site(xp,yp) ) continue;

        res += K(x,y,mu) * K(xp,yp,nu);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> KK_corr(const int mu, const int nu) const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        corr[idx(dx,dy)] = KK_corr( dx, dy, mu, nu );
      }}

    return corr;
  }


  Double K_ss( const Idx dx1, const Idx dy1, const Idx dx2, const Idx dy2, const int mu ) const {
    assert(0<=dx1 && dx1<Lx);
    assert(0<=dy1 && dy1<Ly);
    assert(0<=dx2 && dx2<Lx);
    assert(0<=dy2 && dy2<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        // if( !is_site(x,y) ) continue;
        // if( (x-y+Lx*Ly)%3!=0 ) continue;
        if( !is_link(x,y,mu) ) continue;
        const Idx x1 = (x+dx1)%Lx;
        const Idx y1 = (y+dy1)%Ly;
        const Idx x2 = (x+dx2)%Lx;
        const Idx y2 = (y+dy2)%Ly;
        if( !is_site(x1,y1) || !is_site(x2,y2) ) continue;
        res += K(x,y,mu) * (*this)(x1,y1) * (*this)(x2,y2);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> K_ss_corr( const Idx dx1, const Idx dy1, const int mu ) const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx dx2=0; dx2<Lx; dx2++){
      for(Idx dy2=0; dy2<Ly; dy2++){
        corr[idx(dx2,dy2)] = K_ss( dx1, dy1, dx2, dy2, mu );
      }}
    return corr;
  }



//   Double Txx_even( const Idx x, const Idx y ) const {
//     assert(0<=x && x<Lx);
//     assert(0<=y && y<Ly);
//     assert( get_char(x,y)==0 );

//     Idx xp, yp;
//     Double res = 0.0;
//     int mu;

//     mu = 0;
//     cshift( xp, yp, x, y, mu );
//     res += 2.0*K(x, y, mu);
//     res -= ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     mu = 1;
//     cshift( xp, yp, x, y, mu );
//     res -= K(x, y, mu);
//     res += 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     mu = 2;
//     cshift( xp, yp, x, y, mu );
//     res -= K(x, y, mu);
//     res += 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     res /= 3.0;

//     return res;
//   }

//   Double Txx_odd( const Idx x, const Idx y ) const {
//     Idx xp, yp;
//     Double res = 0.0;
//     int mu;

//     mu = 3;
//     cshift( xp, yp, x, y, mu );
//     res += 2.0*K(x, y, mu);
//     res -= ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     mu = 4;
//     cshift( xp, yp, x, y, mu );
//     res -= K(x, y, mu);
//     res += 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     mu = 5;
//     cshift( xp, yp, x, y, mu );
//     res -= K(x, y, mu);
//     res += 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     res /= 3.0;

//     return res;
//   }

//   Double Txx( const Idx x, const Idx y ) const {
//     assert(0<=x && x<Lx);
//     assert(0<=y && y<Ly);
//     const int c = get_char(x,y);
//     Double res = 0.0;
//     if(c==0) res = Txx_even(x,y);
//     else if(c==2) res = Txx_odd(x,y);
//     else assert(false);
//     return res;
//   }


//   Double Txy_even( const Idx x, const Idx y ) const {
//     assert( get_char(x,y)==0 );

//     Idx xp, yp;
//     Double res = 0.0;
//     int mu;

//     mu = 2;
//     cshift( xp, yp, x, y, mu );
//     res += K(x, y, mu);
//     res -= 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     mu = 1;
//     cshift( xp, yp, x, y, mu );
//     res -= K(x, y, mu);
//     res += 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     res /= std::sqrt(3.0);

//     return res;
//   }


//   Double Txy_odd( const Idx x, const Idx y ) const {
//     assert( get_char(x,y)==2 );

//     Idx xp, yp;
//     Double res = 0.0;
//     int mu;

//     mu = 5;
//     cshift( xp, yp, x, y, mu );
//     res += K(x, y, mu);
//     res -= 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     mu = 4;
//     cshift( xp, yp, x, y, mu );
//     res -= K(x, y, mu);
//     res += 0.5 * ( eps(x,y)+eps(xp,yp) ); // mu deriv

//     res /= std::sqrt(3.0);

//     return res;
//   }

//   Double Txy( const Idx x, const Idx y ) const {
//     assert(0<=x && x<Lx);
//     assert(0<=y && y<Ly);
//     const int c = get_char(x,y);
//     Double res = 0.0;
//     if(c==0) res = Txy_even(x,y);
//     else if(c==2) res = Txy_odd(x,y);
//     else assert(false);
//     return res;
//   }


//   Double Txx_1pt( ) const {
// //     std::vector<Double> tmp(nparallel, 0.0);
// //     std::vector<int> counter(nparallel, 0);

// // #ifdef _OPENMP
// // #pragma omp parallel for num_threads(nparallel)
// // #endif
// //     for(Idx x=0; x<Lx; x++){
// //       for(Idx y=0; y<Ly; y++){
// //         if( !is_link(x,y,mu) ) continue;
// //         tmp[omp_get_thread_num()] += T( x, y, mu );
// //         counter[omp_get_thread_num()]++;
// //       }}

// //     Double res = 0.0;
// //     int tot = 0;
// //     for(int i=0; i<nparallel; i++) {
// //       res += tmp[i];
// //       tot += counter[i];
// //     }
// //     res /= tot;
// //     return res;

//     Double res = 0.0;
//     int counter = 0;

//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         if( !is_site(x,y) ) continue;
//         res += Txx( x, y );
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }


//   Double Txy_1pt( ) const {
// //     std::vector<Double> tmp(nparallel, 0.0);
// //     std::vector<int> counter(nparallel, 0);

// // #ifdef _OPENMP
// // #pragma omp parallel for num_threads(nparallel)
// // #endif
// //     for(Idx x=0; x<Lx; x++){
// //       for(Idx y=0; y<Ly; y++){
// //         if( !is_link(x,y,mu) ) continue;
// //         tmp[omp_get_thread_num()] += T( x, y, mu );
// //         counter[omp_get_thread_num()]++;
// //       }}

// //     Double res = 0.0;
// //     int tot = 0;
// //     for(int i=0; i<nparallel; i++) {
// //       res += tmp[i];
// //       tot += counter[i];
// //     }
// //     res /= tot;
// //     return res;

//     Double res = 0.0;
//     int counter = 0;

//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         if( !is_site(x,y) ) continue;
//         res += Txy( x, y );
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }



//   Double TxxTxx_corr( const Idx dx, const Idx dy ) const {
//     assert(0<=dx && dx<Lx);
//     assert(0<=dy && dy<Ly);

//     Double res = 0.0;
//     int counter = 0;

//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         // if( !is_site(x,y) ) continue;
//         if( (x-y+Lx*Ly)%3!=0 ) continue;
//         const Idx xp = (x+dx)%Lx;
//         const Idx yp = (y+dy)%Ly;
//         if( !is_site(xp,yp) ) continue;
//         res += Txx(x,y) * Txx(xp,yp);
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }

//   Double TxxTxy_corr( const Idx dx, const Idx dy ) const {
//     assert(0<=dx && dx<Lx);
//     assert(0<=dy && dy<Ly);

//     Double res = 0.0;
//     int counter = 0;

//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         // if( !is_site(x,y) ) continue;
//         if( (x-y+Lx*Ly)%3!=0 ) continue;
//         const Idx xp = (x+dx)%Lx;
//         const Idx yp = (y+dy)%Ly;
//         if( !is_site(xp,yp) ) continue;
//         res += Txx(x,y) * Txy(xp,yp);
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }


//   std::vector<Double> TxxTxx_corr() const {
//     std::vector<Double> corr(N, 0.0);

// #ifdef _OPENMP
// #pragma omp parallel for collapse(2) num_threads(nparallel)
//     // #pragma omp parallel for num_threads(nparallel) schedule(static)
// #endif
//     for(Idx dx=0; dx<Lx; dx++){
//       for(Idx dy=0; dy<Ly; dy++){
//         corr[idx(dx,dy)] = TxxTxx_corr( dx, dy );
//       }}
//     return corr;
//   }


//   std::vector<Double> TxxTxy_corr() const {
//     std::vector<Double> corr(N, 0.0);

// #ifdef _OPENMP
// #pragma omp parallel for collapse(2) num_threads(nparallel)
//     // #pragma omp parallel for num_threads(nparallel) schedule(static)
// #endif
//     for(Idx dx=0; dx<Lx; dx++){
//       for(Idx dy=0; dy<Ly; dy++){
//         corr[idx(dx,dy)] = TxxTxy_corr( dx, dy );
//       }}
//     return corr;
//   }


//   Double Txx_ss( const Idx dx1, const Idx dy1, const Idx dx2, const Idx dy2 ) const {
//     assert(0<=dx1 && dx1<Lx);
//     assert(0<=dy1 && dy1<Ly);
//     assert(0<=dx2 && dx2<Lx);
//     assert(0<=dy2 && dy2<Ly);

//     Double res = 0.0;
//     int counter = 0;

//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         // if( !is_site(x,y) ) continue;
//         if( (x-y+Lx*Ly)%3!=0 ) continue;
//         const Idx x1 = (x+dx1)%Lx;
//         const Idx y1 = (y+dy1)%Ly;
//         const Idx x2 = (x+dx2)%Lx;
//         const Idx y2 = (y+dy2)%Ly;
//         if( !is_site(x1,y1) || !is_site(x2,y2) ) continue;
//         res += Txx(x,y) * (*this)(x1,y1) * (*this)(x2,y2);
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }


//   std::vector<Double> Txx_ss_corr( const Idx dx1, const Idx dy1 ) const {
//     std::vector<Double> corr(N, 0.0);

// #ifdef _OPENMP
// #pragma omp parallel for collapse(2) num_threads(nparallel)
//     // #pragma omp parallel for num_threads(nparallel) schedule(static)
// #endif
//     for(Idx dx2=0; dx2<Lx; dx2++){
//       for(Idx dy2=0; dy2<Ly; dy2++){
//         corr[idx(dx2,dy2)] = Txx_ss( dx1, dy1, dx2, dy2 );
//       }}
//     return corr;
//   }


//   Double Txy_ss( const Idx dx1, const Idx dy1, const Idx dx2, const Idx dy2 ) const {
//     assert(0<=dx1 && dx1<Lx);
//     assert(0<=dy1 && dy1<Ly);
//     assert(0<=dx2 && dx2<Lx);
//     assert(0<=dy2 && dy2<Ly);

//     Double res = 0.0;
//     int counter = 0;

//     for(Idx x=0; x<Lx; x++){
//       for(Idx y=0; y<Ly; y++){
//         // if( !is_site(x,y) ) continue;
//         if( (x-y+Lx*Ly)%3!=0 ) continue;
//         const Idx x1 = (x+dx1)%Lx;
//         const Idx y1 = (y+dy1)%Ly;
//         const Idx x2 = (x+dx2)%Lx;
//         const Idx y2 = (y+dy2)%Ly;
//         if( !is_site(x1,y1) || !is_site(x2,y2) ) continue;
//         res += Txy(x,y) * (*this)(x1,y1) * (*this)(x2,y2);
//         counter++;
//       }}

//     res /= counter;
//     return res;
//   }


//   std::vector<Double> Txy_ss_corr( const Idx dx1, const Idx dy1 ) const {
//     std::vector<Double> corr(N, 0.0);

// #ifdef _OPENMP
// #pragma omp parallel for collapse(2) num_threads(nparallel)
//     // #pragma omp parallel for num_threads(nparallel) schedule(static)
// #endif
//     for(Idx dx2=0; dx2<Lx; dx2++){
//       for(Idx dy2=0; dy2<Ly; dy2++){
//         corr[idx(dx2,dy2)] = Txy_ss( dx1, dy1, dx2, dy2 );
//       }}
//     return corr;
//   }



  std::string print() const {
    std::stringstream ss;
    for(Idx y=Ly-1; y>=0; y--){
      for(Idx x= 0; x<Lx; x++) {
        ss << std::setw(5) << (*this)(x, y);
      }
      ss << std::endl;
    }
    return ss.str();
  }



  void ckpoint( const std::string& str ) const {
    std::ofstream of( str, std::ios::out | std::ios::binary | std::ios::trunc);
    if(!of) assert(false);

    int tmp = 0.0;
    for(Idx i=0; i<Lx*Ly; i++){
      tmp = (*this)[i];
      of.write( (char*) &tmp, sizeof(int) );
    }
    of.close();
  }

  void read( const std::string& str ) {
    std::ifstream ifs( str, std::ios::in | std::ios::binary );
    if(!ifs) assert(false);

    int tmp;
    for(Idx i=0; i<Lx*Ly; ++i){
      ifs.read((char*) &tmp, sizeof(int) );
      (*this)[i] = tmp;
    }
    ifs.close();
  }



};


void heatbath( Spin& s ){
  // omp not allowed
  for(Idx i=0; i<Lx*Ly; i++){
    if( !is_site(i) ) continue;

    double henv = 0;
    for(int mu=0; mu<SIX; mu++){
      if( !is_link(i,mu) ) continue;
      Idx j;
      cshift( j, i, mu );
      henv += Beta[mu%THREE]*s[j];
    }

    const Double p = std::exp(2.0*henv);
    const Double r = dist01();
    if( r<p/(1.0+p) ) s[i] = 1;
    else s[i] = -1;
  }
}


void wolff( Spin& s ){
  std::vector<int> is_cluster(Lx*Ly, 0);
  std::stack<Idx> stack_idx;

  Idx init = dist0N();
  while( !is_site(init) ) init = dist0N();

  is_cluster[init] = 1;
  stack_idx.push(init);

  while( stack_idx.size() != 0 ){

    const Idx p = stack_idx.top();
    stack_idx.pop();
    s[p] = -s[p]; // flip when visited

    for(int mu = 0; mu < SIX; mu++){
      if( !is_link(p,mu) ) continue;
      Idx q;
      cshift(q, p, mu);
      if( s[q] == s[p] || is_cluster[q]==1 ) continue; // s[x]*sR[y]<0 or y in c

      const Double r = dist01();
      if( r < std::exp(-2.0 * Beta[mu%THREE]) ) continue; // reject

      is_cluster[q] = 1;
      stack_idx.push(q);
    }
  }
}



struct Scalar {
  Double v;

  Scalar()
    : v(0.0)
  {}

  Scalar( const Double v_ )
    : v(v_)
  {}

  Scalar( const Scalar& other )
    : v(other.v)
  {}

  void clear(){ v = 0.0; }

  Scalar& operator+=(const Scalar& rhs){
    v += rhs.v;
    return *this;
  }

  Scalar& operator+=(const Double& rhs){
    v += rhs;
    return *this;
  }

  Scalar& operator/=(const Double& rhs){
    v /= rhs;
    return *this;
  }

  std::string print() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);
    ss << v;
    return ss.str();
  }

  void print(std::FILE* stream) const {
    fprintf( stream, "%0.15Le\t", v );
  }



};



struct Corr {
  std::vector<Double> v;

  Corr()
    : v(Lx*Ly, 0.0)
  {}

  Corr( const std::vector<Double> v_ )
    : v(v_)
  {}

  Corr( const Corr& other )
    : v(other.v)
  {}

  Double& operator()(const Idx x, const Idx y) { return v[idx(x,y)]; }
  Double operator()(const Idx x, const Idx y) const { return v[idx(x,y)]; }

  void clear(){ for(Idx i=0; i<Lx*Ly; i++) v[i] = 0.0; }

  Corr& operator+=(const Corr& rhs)
  {
    assert( rhs.v.size()==Lx*Ly );
    for(Idx i=0; i<Lx*Ly; i++) v[i] += rhs.v[i];
    return *this;
  }

  Corr& operator+=(const std::vector<Double>& rhs)
  {
    assert( rhs.size()==Lx*Ly );
    for(Idx i=0; i<Lx*Ly; i++) v[i] += rhs[i];
    return *this;
  }

  Corr& operator/=(const Double& rhs)
  {
    for(Idx i=0; i<Lx*Ly; i++) v[i] /= rhs;
    return *this;
  }


  // void print() const {
  //   for(int y=0; y<Ly; y++){
  //     for(int x=0; x<Lx; x++) {
  //       printf( "%0.15e\t", (*this)(x, y) );
  //     }
  //     printf("\n");
  //   }
  // }

  void print(std::FILE* stream) const {
    for(Idx y=0; y<Ly; y++){
      for(Idx x=0; x<Lx; x++) {
        if( !is_site(x,y) ) continue;
        fprintf( stream, "%0.15Le\t", (*this)(x, y) );
      }
      fprintf( stream, "\n");
    }
  }

  std::string print() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);
    for(Idx y=0; y<Ly; y++){
      for(Idx x=0; x<Lx; x++) {
        if( !is_site(x,y) ) continue;
        ss << (*this)(x, y) << " ";
      }
      ss << std::endl;
    }
    return ss.str();
  }

};


template<typename T> // T needs to have: .clear, +=, /= defined
struct Obs {
  std::string description;
  int N;
  std::function<T(const Spin&)> f;

  T sum;
  int counter;

  // Obs() = delete;

  Obs
  (
   const std::string& description_,
   const int N_,
   const std::function<T(const Spin&)>& f_
   )
    : description(description_)
    , N(N_)
    , f(f_)
    , sum()
    , counter(0)
  {}

  void clear(){
    sum.clear();
    counter = 0;
  }

  void meas( const Spin& s ) {
    sum += f(s);
    counter++;
  }

  // T mean() const {
  //   T tmp(sum);
  //   tmp /= counter;
  //   return mean;
  // }

  void write_and_clear( const std::string& dir, const int label ){
    const std::string filename = dir + description + "_" + std::to_string(label) + ".dat";

    // std::ofstream of( filename, std::ios::out | std::ios::trunc );
    // if(!of) assert(false);
    // of << std::scientific << std::setprecision(15);
    // of << sum.print();
    // of.close();

    FILE *stream = fopen(filename.c_str(), "w");
    if (stream == NULL) assert(false);
    std::ofstream of( filename );
    sum /= counter;
    sum.print( stream );
    fclose( stream );

    clear();
  }

};




