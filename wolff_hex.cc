#include <iostream>
#include <iomanip>

#include <omp.h>

#include "header.hpp"


int main( int argc, char *argv[] ){
  std::cout << std::scientific << std::setprecision(15);

#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads( nparallel );
#endif

  Lx = 3;
  Ly = 3;
  // Lx = 6 * 16;
  // Ly = 6 * 32;
  nu = 4; // PP, PA, AA, AP
  K = 2.0/3.0 * std::cos(M_PI/6.0);

  //
  const int Nbin = 102;
  const int binsize = 1e2;
  const int Nconf = Nbin*binsize;
  int ninit = 0;
  //

  // routine
  const int Nheatbath = 2;
  const int Nwolff = 4;
  const int Nrepeat = 50;
  //

  const int seed = 1;
  set_gen( seed );

  // init
  Spin s( Lx*Ly );
  s.random();


  Obs ssA( "ssA", binsize, [](const Spin& s0){ return s0.ss_even(0); } );
  Obs ssB( "ssB", binsize, [](const Spin& s0){ return s0.ss_even(1); } );
  Obs ssC( "ssC", binsize, [](const Spin& s0){ return s0.ss_even(2); } );


  // run
  for(int n=ninit; n<Nconf; n++){

    // update routine
    for(int jj=0; jj<Nrepeat; jj++){
      for(int ii=0; ii<Nheatbath; ii++) heatbath( s );
      for(int ii=0; ii<Nwolff; ii++) wolff( s );
    }

    // measurement
    ssA.meas( s );
    ssB.meas( s );
    ssC.meas( s );

    // write out
    if( (n+1)%binsize==0 ) {
      write( ssA, "./data/", (n+1)/binsize );
      write( ssB, "./data/", (n+1)/binsize );
      write( ssC, "./data/", (n+1)/binsize );
      std::clog << "iter: " << n+1 << std::endl;
    }

  } // end for n

  return 0;
}






// Idx i = 0, j;
// Idx xi, yi, xj, yj;
// get_xy(xi, yi, i);
// std::cout << "i = " << i << ", check = " << idx(xi,yi) << std::endl;

// Idx si = s[i];
// const int mu = 0;
// cshift( j, i, mu );
// get_xy(xj, yj, j);
// std::cout << "check = " << mod(xi-1,Lx) << std::endl;
// std::cout << "xi, yi = " << xi << ", " << yi << std::endl;
// std::cout << "xj, yj = " << xj << ", " << yj << std::endl;


// std::cout << "s = " << std::endl
//           << s.print() << std::endl;

// heatbath( s );
// wolff( s );

// std::cout << "s = " << std::endl
//           << s.print() << std::endl;


// {
//   double eps_tot = 0.0;
//   double h0_tot = 0.0;
//   Idx xx = 0, yy = 0;
//   int mu = 2;
//   int counter = 0;

//   for(int n=0; n<Niter; n++){
//     wolff( s );

//     if( (n+1)%interval==0 ) {
//       eps_tot += s.eps_even(mu);
//       h0_tot += s(xx, yy);
//       counter++;
//     }
//   }

//   std::cout << "K = " << K << std::endl;
//   std::cout << "eps = " << eps_tot/counter << std::endl;
//   std::cout << "h0 = " << h0_tot/counter << std::endl;
// }
