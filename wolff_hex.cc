#include <iostream>
#include <iomanip>
#include <filesystem>

#include <chrono>
#include <limits>

#include <omp.h>

#include "header.hpp"


int main( int argc, char *argv[] ){

  std::cout << std::scientific << std::setprecision(15);
  std::cout << "short = " << std::numeric_limits<short>::digits10 << std::endl;
  std::cout << "int = " << std::numeric_limits<int>::digits10 << std::endl;
  std::cout << "Idx = " << std::numeric_limits<Idx>::digits10 << std::endl;
  // std::cout << "double = " << std::numeric_limits<double>::digits10 << std::endl;
  // std::cout << "long double = " << std::numeric_limits<long double>::digits10 << std::endl;

#ifdef _OPENMP
  // omp_set_dynamic(0);
  omp_set_num_threads( nparallel );
#endif

  // if (argc>1) nu = atoi(argv[1]);
  // if (argc>1) Lx = atoi(argv[2]);
  // if (argc>2) Ly = atoi(argv[3]);
  // if (argc>3) nparallel = atoi(argv[4]);

  const std::string description = "Lx"+std::to_string(Lx)+"Ly"+std::to_string(Ly)+"nu"+std::to_string(nu);
  const std::string datadir = "./data_"+description+"/";
  std::filesystem::create_directories( datadir );
  const std::string configdir = "./config_"+description+"/";
  std::filesystem::create_directories( configdir );

  // routine
  // const int Nbin = 1e4;
  // const int binsize = 1e4;
  const int Nbin = 1e1;
  const int binsize = 1e1;

  const bool if_read = false;
  const bool if_write = true;
  const int init_label = 0;

  const int Nheatbath = 4;
  const int Nwolff = 10;
  const int Nrepeat = 20;

  const int seed = 12;


  // init
  const int ninit = init_label * binsize;
  const int Nconf = Nbin*binsize;
  set_gen( seed );
  Spin s( Lx*Ly );

  if(if_read){
    const std::filesystem::path filepath = static_cast<std::string>(configdir)+"ckpoint"+std::to_string(init_label)+".dat";
    s.read( filepath );
  }
  else s.random();

  {
    std::ofstream of( description+".log", std::ios::out | std::ios::app );
    if(!of) assert(false);
    // of << std::numeric_limits<long double>::digits10 << std::endl;
    of << "Lx = " << Lx << std::endl
       << "Ly = " << Ly << std::endl
       << "nu = " << nu << std::endl
       << "beta_c = " << beta_c << std::endl
       << "Nbin = " << Nbin << std::endl
       << "binsize = " << binsize << std::endl
       << "if_read= " << if_read << std::endl
       << "init_label = " << init_label << std::endl
       << "Nheatbath = " << Nheatbath << std::endl
       << "Nwolff = " << Nwolff << std::endl
       << "Nrepeat = " << Nrepeat << std::endl
       << "seed = " << seed << std::endl;
    of.close();
  }

  // observables
  // Obs<Corr> ss_corr( "ss_corr", binsize v, ss_corr_wrapper );
  // Obs<Scalar> eps_1pt( "eps_1pt", binsize, eps_1pt_wrapper );
  // Obs<Corr> epseps_corr( "epseps_corr", binsize, epseps_corr_wrapper );
  // Obs<Scalar> ss_even0( "ss_even0", binsize, [](const Spin& s0){ return Scalar( s0.ss_even(0) ); } );
  // Obs<Scalar> ss_even1( "ss_even1", binsize, [](const Spin& s0){ return Scalar( s0.ss_even(1) ); } );
  // Obs<Scalar> ss_even2( "ss_even2", binsize, [](const Spin& s0){ return Scalar( s0.ss_even(2) ); } );
  //-----------
  Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );
  Obs<Scalar> eps_1pt( "eps_1pt", binsize, [](const Spin& s0){ return Scalar( s0.eps_1pt() ); } );
  Obs<Corr> epseps_corr( "epseps_corr", binsize, [](const Spin& s0){ return Corr( s0.epseps_corr() ); } );
  // Obs<Scalar> KA( "KA", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(0) ); } );
  // Obs<Scalar> KB( "KB", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(1) ); } );
  // Obs<Scalar> KC( "KC", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(2) ); } );
  Obs<Scalar> Txx( "Txx", binsize, [](const Spin& s0){ return Scalar( s0.Txx_1pt() ); } );
  Obs<Scalar> Txy( "Txy", binsize, [](const Spin& s0){ return Scalar( s0.Txy_1pt() ); } );
  // Obs<Scalar> TC( "TC", binsize, [](const Spin& s0){ return Scalar( s0.T_1pt(2) ); } );
  Obs<Corr> TxxTxx( "TxxTxx", binsize, [](const Spin& s0){ return Corr( s0.TxxTxx_corr() ); } );
  Obs<Corr> TxxTxy( "TxxTxy", binsize, [](const Spin& s0){ return Corr( s0.TxxTxy_corr() ); } );
  // Obs<Corr> TATC( "TATC", binsize, [](const Spin& s0){ return Corr( s0.TT_corr(0,2) ); } );
  // Obs<Corr> TBTA( "TBTA", binsize, [](const Spin& s0){ return Corr( s0.TT_corr(1,0) ); } );
  // Obs<Corr> TBTB( "TBTB", binsize, [](const Spin& s0){ return Corr( s0.TT_corr(1,1) ); } );
  // Obs<Corr> TBTC( "TBTC", binsize, [](const Spin& s0){ return Corr( s0.TT_corr(1,2) ); } );
  // Obs<Corr> TCTA( "TCTA", binsize, [](const Spin& s0){ return Corr( s0.TT_corr(2,0) ); } );
  // Obs<Corr> TCTB( "TCTB", binsize, [](const Spin& s0){ return Corr( s0.TT_corr(2,1) ); } );
  // Obs<Corr> TCTC( "TCTC", binsize, [](const Spin& s0){ return Corr( s0.TT_corr(2,2) ); } );
  //----------------------------


  // run
  auto start = std::chrono::steady_clock::now();
  for(int n=ninit; n<Nconf; n++){

    // update routine
    for(int jj=0; jj<Nrepeat; jj++){
      for(int ii=0; ii<Nheatbath; ii++) heatbath( s );
      for(int ii=0; ii<Nwolff; ii++) wolff( s );
    }

    // measurement
    // ss_even0.meas( s );
    // ss_even1.meas( s );
    // ss_even2.meas( s );
    // -----------
    ss_corr.meas( s );
    eps_1pt.meas( s );
    epseps_corr.meas( s );
    // KA.meas( s );
    // KB.meas( s );
    // KC.meas( s );
    Txx.meas( s );
    Txy.meas( s );
    // TC.meas( s );
    // //
    TxxTxx.meas( s );
    TxxTxy.meas( s );
    // TATC.meas( s );
    // TBTA.meas( s );
    // TBTB.meas( s );
    // TBTC.meas( s );
    // TCTA.meas( s );
    // TCTB.meas( s );
    // TCTC.meas( s );
    //------
    // TCTC.meas( s );

    // write out and clear
    if( (n+1)%binsize==0 && if_write ) {
      // ss_even0.write_and_clear( datadir, (n+1)/binsize );
      // ss_even1.write_and_clear( datadir, (n+1)/binsize );
      // ss_even2.write_and_clear( datadir, (n+1)/binsize );
      // -----------
      ss_corr.write_and_clear( datadir, (n+1)/binsize );
      eps_1pt.write_and_clear( datadir, (n+1)/binsize );
      epseps_corr.write_and_clear( datadir, (n+1)/binsize );
      // TA.write_and_clear( datadir, (n+1)/binsize );
      // TB.write_and_clear( datadir, (n+1)/binsize );
      // TC.write_and_clear( datadir, (n+1)/binsize );
      //
      // KA.write_and_clear( datadir, (n+1)/binsize );
      // KB.write_and_clear( datadir, (n+1)/binsize );
      // KC.write_and_clear( datadir, (n+1)/binsize );
      //
      Txx.write_and_clear( datadir, (n+1)/binsize );
      Txy.write_and_clear( datadir, (n+1)/binsize );
      //
      TxxTxx.write_and_clear( datadir, (n+1)/binsize );
      TxxTxy.write_and_clear( datadir, (n+1)/binsize );
      // TATC.write_and_clear( datadir, (n+1)/binsize );
      // TBTA.write_and_clear( datadir, (n+1)/binsize );
      // TBTB.write_and_clear( datadir, (n+1)/binsize );
      // TBTC.write_and_clear( datadir, (n+1)/binsize );
      // TCTA.write_and_clear( datadir, (n+1)/binsize );
      // TCTB.write_and_clear( datadir, (n+1)/binsize );
      // TCTC.write_and_clear( datadir, (n+1)/binsize );
      // -------
      // TCTC.write_and_clear( datadir, (n+1)/binsize );
      //
      std::clog << "iter: " << n+1 << std::endl;

      const int label = (n+1)/binsize;
      // const std::filesystem::path filepath = static_cast<std::string>(configdir)+"ckpoint"+std::to_string(label)+".dat";
      const std::string filepath = configdir+"ckpoint"+std::to_string(label)+".dat";
      s.ckpoint( filepath );

      const auto end = std::chrono::steady_clock::now();
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << elapsed_seconds.count() << " sec" << std::endl;
      start = std::chrono::steady_clock::now();
    }

  } // end for n


  return 0;
}

