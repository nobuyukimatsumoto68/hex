#include <iostream>
#include <iomanip>
#include <filesystem>

#include <chrono>
#include <limits>

#include <omp.h>

#include "header2.hpp"


int main( int argc, char *argv[] ){

#ifdef _OPENMP
  // omp_set_dynamic(0);
  omp_set_num_threads( 1 );
#endif

  mult = 8;
  Lx = 3*1*mult; // 12
  Ly = 3*1*mult;

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
  std::cout << "filepath = " << filepath << std::endl;

  // observables
  int binsize=1;
  //-----------
  // Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );

  // Obs<Scalar> KA( "KA", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(0) ); } );
  // Obs<Scalar> KB( "KB", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(1) ); } );
  // Obs<Scalar> KC( "KC", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(2) ); } );

  // Obs<Corr> KAKA( "KAKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,0) ); } );
  // Obs<Corr> KAKB( "KAKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,1) ); } );
  // Obs<Corr> KAKC( "KAKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,2) ); } );
  // //
  // Obs<Corr> KBKA( "KBKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,0) ); } );
  // Obs<Corr> KBKB( "KBKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,1) ); } );
  // Obs<Corr> KBKC( "KBKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,2) ); } );
  // //
  // Obs<Corr> KCKA( "KCKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,0) ); } );
  // Obs<Corr> KCKB( "KCKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,1) ); } );
  // Obs<Corr> KCKC( "KCKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,2) ); } );

  // Obs<Corr> KA_ss( "KA_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 0) ); } );
  // Obs<Corr> KB_ss( "KB_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 1) ); } );
  // Obs<Corr> KC_ss( "KC_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 2) ); } );

  Idx x1 = 0, y1 = 0, x2 = Lx/2, y2 = 0;
  Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );
  Obs<Scalar> eps_1pt( "eps_1pt", binsize, [](const Spin& s0){ return Scalar( s0.eps_1pt() ); } );
  Obs<Corr> epseps_corr( "epseps_corr", binsize, [](const Spin& s0){ return Corr( s0.epseps_corr() ); } );
  Obs<Scalar> TA( "TA", binsize, [](const Spin& s0){ return Scalar( s0.TM_1pt(0) ); } );
  Obs<Scalar> TB( "TB", binsize, [](const Spin& s0){ return Scalar( s0.TM_1pt(1) ); } );
  Obs<Scalar> TC( "TC", binsize, [](const Spin& s0){ return Scalar( s0.TM_1pt(2) ); } );

  Obs<Corr> TATA( "TATA", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(0,0) ); } );
  Obs<Corr> TATB( "TATB", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(0,1) ); } );
  Obs<Corr> TATC( "TATC", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(0,2) ); } );
  Obs<Corr> TBTA( "TBTA", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(1,0) ); } );
  Obs<Corr> TBTB( "TBTB", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(1,1) ); } );
  Obs<Corr> TBTC( "TBTC", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(1,2) ); } );
  Obs<Corr> TCTA( "TCTA", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(2,0) ); } );
  Obs<Corr> TCTB( "TCTB", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(2,1) ); } );
  Obs<Corr> TCTC( "TCTC", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(2,2) ); } );

  // Obs<Corr> eps_ss( "eps_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.eps_ss_corr(x1, y1, x2, y2) ); } );
  Obs<Corr> TA_ss( "TA_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.TM_ss_corr(0, x1, y1, x2, y2) ); } );
  Obs<Corr> TB_ss( "TB_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.TM_ss_corr(1, x1, y1, x2, y2) ); } );
  Obs<Corr> TC_ss( "TC_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.TM_ss_corr(2, x1, y1, x2, y2) ); } );

  // Obs<Corr> Txx_ss( "Txx_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.Txx_ss_corr(x1, y1, x2, y2) ); } );
  // Obs<Corr> Txy_ss( "Txy_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.Txy_ss_corr(x1, y1, x2, y2) ); } );
  // Obs<Corr> Txx_epseps( "Txx_epseps", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.Txx_epseps_corr(x1, y1, x2, y2) ); } );
  // Obs<Corr> Txy_epseps( "Txy_epseps", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.Txy_epseps_corr(x1, y1, x2, y2) ); } );

  //----------------------------

  // Idx x1 = 0, y1 = 0, x2 = Lx/2, y2 = 0;
  // const double B = 1.2990381056766578;
  // std::cout << "KA_ss = " << std::endl
  //           << ( Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 0) )*=B ).print() << std::endl;

  // std::cout << "spin = " << std::endl
  //           << s.print() << std::endl;

  {
    const int mu =0;
    std::cout << "B = " << DBetaDKappa[mu%3] << std::endl;
    std::cout << "div = " << -DBetaDKappa[mu%3] * std::tanh(Beta[mu%3]) - divs[mu%3] << std::endl;
  }

  // std::cout << "KA = " << std::endl
  //           << s.Knew_1pt(0) << std::endl;

  // std::cout << "KB = " << std::endl
  //           << s.Knew_1pt(1) << std::endl;

  // std::cout << "KC = " << std::endl
  //           << s.Knew_1pt(2) << std::endl;

  // std::cout << "eps = " << std::endl
  //           << s.eps_1pt() << std::endl;

  // std::cout << "TA = " << std::endl
  //           << s.TM_1pt(0) << std::endl;

  // std::cout << "TB = " << std::endl
  //           << s.TM_1pt(1) << std::endl;

  // std::cout << "TC = " << std::endl
  //           << s.TM_1pt(2) << std::endl;

  // std::cout << "TA_ss = " << std::endl
  //           << ( Corr( s.TM_ss_corr(0, 0.0, 0.0, Lx/2, 0.0) ) ).print() << std::endl;
  // std::cout << "TB_ss = " << std::endl
  //           << ( Corr( s.TM_ss_corr(1, 0.0, 0.0, Lx/2, 0.0) ) ).print() << std::endl;
  std::cout << "TC_ss = " << std::endl
            << ( Corr( s.TM_ss_corr(2, 0.0, 0.0, Lx/2, 0.0) ) ).print() << std::endl;

  // auto corrA = Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 0) )*=B;
  // auto corrB = Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 1) )*=B;
  // auto corrC = Corr( s.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 2) )*=B;

  // corrA += corrB;
  // corrA += corrC;
  // std::cout << "sum = " << std::endl
  //           << corrA.print() << std::endl;

  // corrA *= 0.5 * 2.0/3.0;

  // auto corr_ss = Corr( s.ss_corr() );
  // std::cout << "corr_ss = " << std::endl
  //           << corr_ss.print() << std::endl;

  // corrA += (corr_ss(x2-x1, y2-y1)*=(-1.0));
  // std::cout << "subt = " << std::endl
  //           << corrA.print() << std::endl;

  // std::cout << "KA_ss = " << std::endl
  //           << Corr( s.K_ss_corr(0, x1, y1, x2, y2) ).print() << std::endl;

  // auto corrA = Corr( s.K_ss_corr(0, x1, y1, x2, y2) );
  // auto corrB = Corr( s.K_ss_corr(1, x1, y1, x2, y2) );
  // auto corrC = Corr( s.K_ss_corr(2, x1, y1, x2, y2) );

  // corrA += corrB;
  // corrA += corrC;
  // std::cout << "sum = " << std::endl
  //           << corrA.print() << std::endl;

  // auto corr_ss = Corr( s.ss_corr() );
  // std::cout << "corr_ss = " << std::endl
  //           << corr_ss.print() << std::endl;

  // auto corr_ss = Corr( s.ss_corr() );
  // std::cout << "Ess = " << std::endl
  //           << corrE.print() << std::endl;

  // auto corrE = Corr( s.eps_ss_corr(x1, y1, x2, y2) );
  // std::cout << "Ess = " << std::endl
  //           << corrE.print() << std::endl;

  // auto TATA_corr = Corr( s.TMTN_corr(0, 0) );
  // std::cout << "TATA_corr = " << std::endl
  //           << TATA_corr.print() << std::endl;


  return 0;
}

