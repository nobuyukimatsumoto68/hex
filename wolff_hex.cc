#include <iostream>
#include <iomanip>
#include <filesystem>

#include <omp.h>

#include "header.hpp"


int main( int argc, char *argv[] ){
  std::cout << std::scientific << std::setprecision(15);

#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads( nparallel );
#endif

  // Lx = 6 * 8;
  // Ly = 6 * 16;

  // if (argc>1) nu = atoi(argv[1]);
  // if (argc>1) Lx = atoi(argv[2]);
  // if (argc>2) Ly = atoi(argv[3]);
  // if (argc>3) nparallel = atoi(argv[4]);

  const std::string description = "Lx"+std::to_string(Lx)+"Ly"+std::to_string(Ly)+"nu"+std::to_string(nu);
  const std::filesystem::path datadir = "./data_"+description+"/";
  std::filesystem::create_directories( datadir );

  //
  // const int Nbin = 1;
  // const int binsize = 1e1;
  //
  //

  // routine
  const int Nbin = 104;
  const int binsize = 1e3;
  // const int Nbin = 1;
  // const int binsize = 1e1;

  const int Nconf = Nbin*binsize;
  int ninit = 0;


  const int Nheatbath = 5;
  const int Nwolff = 5;
  const int Nrepeat = 20;
  //

  const int seed = 1;
  set_gen( seed );

  // init
  Spin s( Lx*Ly );
  s.random();

  // observables
  // Obs<Corr> ss_corr( "ss_corr", binsize, ss_corr_wrapper );
  // Obs<Scalar> eps_1pt( "eps_1pt", binsize, eps_1pt_wrapper );
  // Obs<Corr> epseps_corr( "epseps_corr", binsize, epseps_corr_wrapper );
  Obs<Scalar> ss_even0( "ss_even0", binsize, [](const Spin& s0){ return Scalar( s0.ss_even(0) ); } );
  Obs<Scalar> ss_even1( "ss_even1", binsize, [](const Spin& s0){ return Scalar( s0.ss_even(1) ); } );
  Obs<Scalar> ss_even2( "ss_even2", binsize, [](const Spin& s0){ return Scalar( s0.ss_even(2) ); } );
  Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );
  Obs<Scalar> eps_1pt( "eps_1pt", binsize, [](const Spin& s0){ return Scalar( s0.eps_1pt() ); } );
  Obs<Corr> epseps_corr( "epseps_corr", binsize, [](const Spin& s0){ return Corr( s0.epseps_corr() ); } );
  Obs<Scalar> TA( "TA", binsize, [](const Spin& s0){ return Scalar( s0.T_1pt(0) ); } );
  Obs<Scalar> TB( "TB", binsize, [](const Spin& s0){ return Scalar( s0.T_1pt(1) ); } );
  Obs<Scalar> TC( "TC", binsize, [](const Spin& s0){ return Scalar( s0.T_1pt(2) ); } );

  // run
  for(int n=ninit; n<Nconf; n++){

    // update routine
    for(int jj=0; jj<Nrepeat; jj++){
      for(int ii=0; ii<Nheatbath; ii++) heatbath( s );
      for(int ii=0; ii<Nwolff; ii++) wolff( s );
    }

    // measurement
    ss_even0.meas( s );
    ss_even1.meas( s );
    ss_even2.meas( s );
    ss_corr.meas( s );
    eps_1pt.meas( s );
    epseps_corr.meas( s );
    TA.meas( s );
    TB.meas( s );
    TC.meas( s );

    // write out and clear
    if( (n+1)%binsize==0 ) {
      ss_even0.write_and_clear( datadir, (n+1)/binsize );
      ss_even1.write_and_clear( datadir, (n+1)/binsize );
      ss_even2.write_and_clear( datadir, (n+1)/binsize );
      //
      ss_corr.write_and_clear( datadir, (n+1)/binsize );
      eps_1pt.write_and_clear( datadir, (n+1)/binsize );
      epseps_corr.write_and_clear( datadir, (n+1)/binsize );
      TA.write_and_clear( datadir, (n+1)/binsize );
      TB.write_and_clear( datadir, (n+1)/binsize );
      TC.write_and_clear( datadir, (n+1)/binsize );
      //
      std::clog << "iter: " << n+1 << std::endl;
    }

  } // end for n

  return 0;
}

