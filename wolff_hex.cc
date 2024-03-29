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

  Lx = 6 * 16;
  Ly = 6 * 32;

  // if (argc>1) nu = atoi(argv[1]);
  if (argc>1) Lx = atoi(argv[2]);
  if (argc>2) Ly = atoi(argv[3]);
  if (argc>3) nparallel = atoi(argv[4]);

  const std::string description = "Lx"+std::to_string(Lx)+"Ly"+std::to_string(Ly)+"nu"+std::to_string(nu);
  const std::filesystem::path datadir = "./data_"+description+"/";
  std::filesystem::create_directories( datadir );

  //
  const int Nbin = 42;
  const int binsize = 1e3;
  //
  const int Nconf = Nbin*binsize;
  int ninit = 0;
  //

  // routine
  const int Nheatbath = 2;
  const int Nwolff = 4;
  const int Nrepeat = 10;
  //

  const int seed = 1;
  set_gen( seed );

  // init
  Spin s( Lx*Ly );
  s.random();

  // observables
  Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );

  // run
  for(int n=ninit; n<Nconf; n++){

    // update routine
    for(int jj=0; jj<Nrepeat; jj++){
      for(int ii=0; ii<Nheatbath; ii++) heatbath( s );
      for(int ii=0; ii<Nwolff; ii++) wolff( s );
    }

    // measurement
    ss_corr.meas( s );

    // write out and clear
    if( (n+1)%binsize==0 ) {
      ss_corr.write_and_clear( datadir, (n+1)/binsize );
      //
      std::clog << "iter: " << n+1 << std::endl;
    }

  } // end for n

  return 0;
}

