/* ============================================
   NM mtsmtnbyk@gmail.com
   ============================================= */
#define Debug 0

/* ---- FIXED ---- */
#define Two 2         // spacetime dimension=2
#define Three 3       // triangular lattice

/* ---- LATTICE ---- */
#define Lx 4
#define Ly 8
#define N Three*Lx*Ly     // 3 * Lx * Ly

/* ---- RANDOM ---- */
#define Seed 1

/* ---- MONTE CARLO ---- */
#define NConf 100
#define Interval 100

#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>


double beta = 1.0;
double hExt = 0.0;


void printSpins(int s[Ly][Lx]);
double heatBath(int s[Ly][Lx], std::mt19937 gen[Ly][Lx]);
double getM(int s[Ly][Lx]);


// random number
std::uniform_int_distribution<> dist; // [a,b] = [0,2147483647]; we check a%2=0, b%2=1.
int uniformPM1( std::mt19937& gen ){ return 2 * (dist(gen)%2) - 1; }
double uniform0to1( std::mt19937& gen ){
  double res = dist(gen);
  res /= dist.b();
  return res; }
inline int mod(int x, int n){ return  (n + (x % n))%n; }


int main( int argc, char *argv[] )
{
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);
  if(argc != 3) std::clog << "proceed with the default value." << std::endl;
  else{
    beta = std::atof(argv[1]);
    hExt = std::atof(argv[2]);
  }

  assert( dist.a()==0 );
  assert( dist.b()%2==1 );
  std::mt19937 gen[Ly][Lx];
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      gen[y][x].seed(Seed + Lx*y + x);
    }}

  std::clog <<
    "Lx: " << Lx << std::endl <<
    "Ly:  " << Ly << std::endl <<
    "Seed: " << Seed << std::endl <<
    "beta: " << beta << std::endl <<
    "hExt: " << hExt << std::endl <<
    "NConf: " << NConf << std::endl <<
    "Interval: " << Interval << std::endl
    ;

  // --------------------------------------


  /* === Site Variables === */
  int s[Ly][Lx];

  // initialize
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      s[y][x] = uniformPM1(gen[y][x]);
    }}
  printSpins( s );

  // // update & record
  // for(int k=0; k<NConf; k++) {
  //   const double r = heatBath(s);
  //   if(k%Interval==0) {
  //     const double b = getM(s);
  //     std::cout << std::setw (22) << b << " " << r << std::endl;
  //   }
  // }

  return 0;
}



double hEff(const int s[Ly][Lx],
            const int y, const int x){
  const int xp = mod(x+1,Lx);
  const int yp = mod(y+1,Ly);
  const int xm = mod(x-1,Lx);
  const int ym = mod(y-1,Ly);
  const int mEff = s[yp][x] + s[y][xp] + s[xp][yp] + s[ym][x] + s[y][xm] + s[xm][ym];
  return beta*mEff + hExt;
}



double heatBath(int s[Ly][Lx], std::mt19937 gen[Ly][Lx]){
  double rate = 0.0;
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      const double heff = hEff(s, y, x);
      const double p = std::exp(2.0*heff);
      const double r = uniform0to1(gen[y][x]);
      if( r<p/(1.0+p) ) {
        s[y][x] = +1;
        rate++;
      }
      else {
        s[y][x] = -1;
      }
    }}
  return rate/(Lx*Ly);
}


double getM(int s[Ly][Lx]){
  double res = 0;
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      res += s[y][x];
    }}
  return res / (Ly*Lx);
}



/********************************************************/


void printSpins(int s[Ly][Lx])
{
  for(int y=Ly-1; y>=0; y--){
    std::cout << "yg= " << std::setw(4) << y << ":   ";
    for(int x= 0 ; x<Lx; x++){
      std::cout << std::setw(5) << s[y][x];
    }
    std::cout << std::endl;
  }
}
