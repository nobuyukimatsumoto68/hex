#include <iostream>
#include <iomanip>
#include <filesystem>

#include <chrono>
#include <limits>

#include <omp.h>

#include "header.hpp"


int main( int argc, char *argv[] ){

#ifdef _OPENMP
  // omp_set_dynamic(0);
  omp_set_num_threads( nparallel );
#endif

  set_all();

  Spin s( Lx*Ly );
  std::filesystem::path dir{"/projectnb/qfe/nmatsum/hex_gen/mult8/config_Lx24Ly24nu1tautil0.500000_0.866025/"};
  std::string file;
  int init_label=0;
  for(auto const& dir_entry : std::filesystem::directory_iterator{dir}){
    file = dir_entry.path();
    const std::size_t c1 = file.find("ckpoint");
    const std::size_t c2 = file.find(".dat");
    const int tmp = std::stoi( file.substr(c1+7, c2-c1-7) );
    if(init_label < tmp) init_label = tmp;
  }
  const std::filesystem::path filepath = static_cast<std::string>(dir)+"ckpoint"+std::to_string(init_label)+".dat";
  s.read( filepath );

  // observables
  int binsize=1;
  //-----------
  Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );

  Obs<Scalar> KA( "KA", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(0) ); } );
  Obs<Scalar> KB( "KB", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(1) ); } );
  Obs<Scalar> KC( "KC", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(2) ); } );

  Obs<Corr> KAKA( "KAKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,0) ); } );
  Obs<Corr> KAKB( "KAKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,1) ); } );
  Obs<Corr> KAKC( "KAKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,2) ); } );
  //
  Obs<Corr> KBKA( "KBKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,0) ); } );
  Obs<Corr> KBKB( "KBKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,1) ); } );
  Obs<Corr> KBKC( "KBKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,2) ); } );
  //
  Obs<Corr> KCKA( "KCKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,0) ); } );
  Obs<Corr> KCKB( "KCKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,1) ); } );
  Obs<Corr> KCKC( "KCKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,2) ); } );

  Obs<Corr> KA_ss( "KA_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 0) ); } );
  Obs<Corr> KB_ss( "KB_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 1) ); } );
  Obs<Corr> KC_ss( "KC_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 2) ); } );
  //----------------------------

  Idx x1 = 0, y1 = 0, x2 = Lx/2, y2 = 0;
  const double B = 1.2990381056766578;
  std::cout << "KA_ss = " << std::endl
            << ( Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 0) )*=B ).print() << std::endl;

  auto corrA = Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 0) )*=B;
  auto corrB = Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 1) )*=B;
  auto corrC = Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 2) )*=B;

  corrA += corrB;
  corrA += corrC;
  std::cout << "sum = " << std::endl
            << corrA.print() << std::endl;

  corrA *= 0.5 * 2.0/3.0;

  auto corr_ss = Corr( s.ss_corr() );
  std::cout << "corr_ss = " << std::endl
            << corr_ss.print() << std::endl;

  corrA += (corr_ss(x2-x1, y2-y1)*=(-1.0));
  std::cout << "subt = " << std::endl
            << corrA.print() << std::endl;


  return 0;
}

