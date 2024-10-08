#include <iostream>
#include <iomanip>
#include <filesystem>

#include <chrono>
#include <limits>

#include <omp.h>

#include "header.hpp"


int main( int argc, char *argv[] ){

  int init_label = 0;
  bool if_read = false;
  if (argc>1){
    mult = atoi(argv[1]);
    Lx = 3*mult; // 12
    Ly = 3*mult;
  }
  if (argc>2){
    if_read=true;
    init_label = atoi(argv[2]);
  }

  std::cout << std::scientific << std::setprecision(15);
  std::cout << "int = " << std::numeric_limits<int>::digits10 << std::endl;
  std::cout << "Idx = " << std::numeric_limits<Idx>::digits10 << std::endl;
  std::cout << "double = " << std::numeric_limits<double>::digits10 << std::endl;
  std::cout << "Double = " << std::numeric_limits<Double>::digits10 << std::endl;

#ifdef _OPENMP
  // omp_set_dynamic(0);
  omp_set_num_threads( nparallel );
#endif

  set_all();
  const std::string description = "Lx"+std::to_string(Lx)+"Ly"+std::to_string(Ly)+"nu"+std::to_string(nu)+"tautil"+std::to_string(tautil1)+"_"+std::to_string(tautil2);
  const std::string datadir = "./data_"+description+"/";
  std::filesystem::create_directories( datadir );
  const std::string configdir = "./config_"+description+"/";
  std::filesystem::create_directories( configdir );

  // routine
  const int Nbin = 1e4;
  const int binsize = 1e4;
  // const int Nbin = 10;
  // const int binsize = 1;
  // const int Nbin = 1e2;
  // const int binsize = 1e2;

  const bool if_write = true;

  const int Nheatbath = 4;
  const int Nwolff = 10;
  // const int Nrepeat = 20;
  const int Nrepeat = 4;

  const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  // std::cout << seed << std::endl;

  // init
  const int Nconf = Nbin*binsize;
  set_gen( seed );
  Spin s( Lx*Ly );

  if(if_read){
    std::filesystem::path dir{configdir};
    std::string file;
    for(auto const& dir_entry : std::filesystem::directory_iterator{dir}){
      file = dir_entry.path();
      const std::size_t c1 = file.find("ckpoint");
      const std::size_t c2 = file.find(".dat");
      const int tmp = std::stoi( file.substr(c1+7, c2-c1-7) );
      if(init_label < tmp) init_label = tmp;
    }
    const std::filesystem::path filepath = static_cast<std::string>(configdir)+"ckpoint"+std::to_string(init_label)+".dat";
    s.read( filepath );
  }
  else s.random();

  int ninit = init_label * binsize;

  {
    std::ofstream of( description+".log", std::ios::out | std::ios::app );
    if(!of) assert(false);
    of << "Lx = " << Lx << std::endl
       << "Ly = " << Ly << std::endl
       << "nu = " << nu << std::endl
       // << "beta_c = " << beta_c << std::endl
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

  {
    std::ofstream of;
    of.open( description+"ell.dat", std::ios::out | std::ios::trunc);
    if(!of) assert(false);
    of << std::scientific << std::setprecision(15);

    of << ell0[0] << ", " << ell0[1] << std::endl
       << ell1[0] << ", " << ell1[1] << std::endl
       << ell2[0] << ", " << ell2[1] << std::endl;
  }

  {
    std::ofstream of;
    of.open( description+"kappa.dat", std::ios::out | std::ios::trunc);
    if(!of) assert(false);
    of << std::scientific << std::setprecision(15);

    of << kappa[0] << ", " << kappa[1] << ", " << kappa[2] << std::endl;
  }

  {
    std::ofstream of;
    of.open( description+"ellstar.dat", std::ios::out | std::ios::trunc);
    if(!of) assert(false);
    of << std::scientific << std::setprecision(15);

    of << ell_star0[0] << ", " << ell_star0[1] << std::endl
       << ell_star1[0] << ", " << ell_star1[1] << std::endl
       << ell_star2[0] << ", " << ell_star2[1] << std::endl;
  }

  {
    std::ofstream of;
    of.open( description+"e.dat", std::ios::out | std::ios::trunc);
    if(!of) assert(false);
    of << std::scientific << std::setprecision(15);

    of << e0[0] << ", " << e0[1] << std::endl
       << e1[0] << ", " << e1[1] << std::endl
       << e2[0] << ", " << e2[1] << std::endl;

  }

  {
    std::ofstream of;
    of.open( description+"cosH.dat", std::ios::out | std::ios::trunc);
    if(!of) assert(false);
    of << std::scientific << std::setprecision(15);

    of << cosH[0] << ", " << cosH[1] << ", " << cosH[2] << std::endl;
  }

  {
    std::ofstream of;
    of.open( description+"beta.dat", std::ios::out | std::ios::trunc);
    if(!of) assert(false);
    of << std::scientific << std::setprecision(15);

    of << Beta[0] << ", " << Beta[1] << ", " << Beta[2] << std::endl;
  }


  // observables
  //-----------
  Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );

  Obs<Scalar> KA( "KA", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(0) ); } );
  Obs<Scalar> KB( "KB", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(1) ); } );
  Obs<Scalar> KC( "KC", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(2) ); } );
  // Obs<Scalar> KD( "KD", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(3) ); } );
  // Obs<Scalar> KE( "KE", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(4) ); } );
  // Obs<Scalar> KF( "KF", binsize, [](const Spin& s0){ return Scalar( s0.K_1pt(5) ); } );

  Obs<Corr> KAKA( "KAKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,0) ); } );
  Obs<Corr> KAKB( "KAKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,1) ); } );
  Obs<Corr> KAKC( "KAKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,2) ); } );
  // Obs<Corr> KAKD( "KAKD", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,3) ); } );
  // Obs<Corr> KAKE( "KAKE", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,4) ); } );
  // Obs<Corr> KAKF( "KAKF", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(0,5) ); } );
  //
  Obs<Corr> KBKA( "KBKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,0) ); } );
  Obs<Corr> KBKB( "KBKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,1) ); } );
  Obs<Corr> KBKC( "KBKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,2) ); } );
  // Obs<Corr> KBKD( "KBKD", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,3) ); } );
  // Obs<Corr> KBKE( "KBKE", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,4) ); } );
  // Obs<Corr> KBKF( "KBKF", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(1,5) ); } );
  //
  Obs<Corr> KCKA( "KCKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,0) ); } );
  Obs<Corr> KCKB( "KCKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,1) ); } );
  Obs<Corr> KCKC( "KCKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,2) ); } );
  // Obs<Corr> KCKD( "KCKD", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,3) ); } );
  // Obs<Corr> KCKE( "KCKE", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,4) ); } );
  // Obs<Corr> KCKF( "KCKF", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(2,5) ); } );
  //
  // Obs<Corr> KDKA( "KDKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(3,0) ); } );
  // Obs<Corr> KDKB( "KDKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(3,1) ); } );
  // Obs<Corr> KDKC( "KDKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(3,2) ); } );
  // Obs<Corr> KDKD( "KDKD", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(3,3) ); } );
  // Obs<Corr> KDKE( "KDKE", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(3,4) ); } );
  // Obs<Corr> KDKF( "KDKF", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(3,5) ); } );
  // //
  // Obs<Corr> KEKA( "KEKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(4,0) ); } );
  // Obs<Corr> KEKB( "KEKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(4,1) ); } );
  // Obs<Corr> KEKC( "KEKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(4,2) ); } );
  // Obs<Corr> KEKD( "KEKD", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(4,3) ); } );
  // Obs<Corr> KEKE( "KEKE", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(4,4) ); } );
  // Obs<Corr> KEKF( "KEKF", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(4,5) ); } );
  // //
  // Obs<Corr> KFKA( "KFKA", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(5,0) ); } );
  // Obs<Corr> KFKB( "KFKB", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(5,1) ); } );
  // Obs<Corr> KFKC( "KFKC", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(5,2) ); } );
  // Obs<Corr> KFKD( "KFKD", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(5,3) ); } );
  // Obs<Corr> KFKE( "KFKE", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(5,4) ); } );
  // Obs<Corr> KFKF( "KFKF", binsize, [](const Spin& s0){ return Corr( s0.KK_corr(5,5) ); } );

  Obs<Corr> KA_ss( "KA_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 0) ); } );
  Obs<Corr> KB_ss( "KB_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 1) ); } );
  Obs<Corr> KC_ss( "KC_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(0.0, 0.0, Lx/2, 0.0, 2) ); } );
  // Obs<Corr> KA_ss( "KA_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(Lx/2, 0.0, 0) ); } );
  // Obs<Corr> KB_ss( "KB_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(Lx/2, 0.0, 1) ); } );
  // Obs<Corr> KC_ss( "KC_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(Lx/2, 0.0, 2) ); } );
  // Obs<Corr> KD_ss( "KD_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(Lx/2, 0.0, 3) ); } );
  // Obs<Corr> KE_ss( "KE_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(Lx/2, 0.0, 4) ); } );
  // Obs<Corr> KF_ss( "KF_ss", binsize, [](const Spin& s0){ return Corr( s0.K_ss_corr(Lx/2, 0.0, 5) ); } );
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
    ss_corr.meas( s );

    KA.meas( s );
    KB.meas( s );
    KC.meas( s );
    // KD.meas( s );
    // KE.meas( s );
    // KF.meas( s );

    KAKA.meas( s );
    KAKB.meas( s );
    KAKC.meas( s );
    // KAKD.meas( s );
    // KAKE.meas( s );
    // KAKF.meas( s );
    //
    KBKA.meas( s );
    KBKB.meas( s );
    KBKC.meas( s );
    // KBKD.meas( s );
    // KBKE.meas( s );
    // KBKF.meas( s );
    //
    KCKA.meas( s );
    KCKB.meas( s );
    KCKC.meas( s );
    // KCKD.meas( s );
    // KCKE.meas( s );
    // KCKF.meas( s );
    //
    // KDKA.meas( s );
    // KDKB.meas( s );
    // KDKC.meas( s );
    // KDKD.meas( s );
    // KDKE.meas( s );
    // KDKF.meas( s );
    // //
    // KEKA.meas( s );
    // KEKB.meas( s );
    // KEKC.meas( s );
    // KEKD.meas( s );
    // KEKE.meas( s );
    // KEKF.meas( s );
    // //
    // KFKA.meas( s );
    // KFKB.meas( s );
    // KFKC.meas( s );
    // KFKD.meas( s );
    // KFKE.meas( s );
    // KFKF.meas( s );

    KA_ss.meas( s );
    KB_ss.meas( s );
    KC_ss.meas( s );
    // KD_ss.meas( s );
    // KE_ss.meas( s );
    // KF_ss.meas( s );
    //------

    // write out and clear
    if( (n+1)%binsize==0 && if_write ) {
      // -----------
      ss_corr.write_and_clear( datadir, (n+1)/binsize );

      KA.write_and_clear( datadir, (n+1)/binsize );
      KB.write_and_clear( datadir, (n+1)/binsize );
      KC.write_and_clear( datadir, (n+1)/binsize );
      // KD.write_and_clear( datadir, (n+1)/binsize );
      // KE.write_and_clear( datadir, (n+1)/binsize );
      // KF.write_and_clear( datadir, (n+1)/binsize );

      KAKA.write_and_clear( datadir, (n+1)/binsize );
      KAKB.write_and_clear( datadir, (n+1)/binsize );
      KAKC.write_and_clear( datadir, (n+1)/binsize );
      // KAKD.write_and_clear( datadir, (n+1)/binsize );
      // KAKE.write_and_clear( datadir, (n+1)/binsize );
      // KAKF.write_and_clear( datadir, (n+1)/binsize );

      KBKA.write_and_clear( datadir, (n+1)/binsize );
      KBKB.write_and_clear( datadir, (n+1)/binsize );
      KBKC.write_and_clear( datadir, (n+1)/binsize );
      // KBKD.write_and_clear( datadir, (n+1)/binsize );
      // KBKE.write_and_clear( datadir, (n+1)/binsize );
      // KBKF.write_and_clear( datadir, (n+1)/binsize );

      KCKA.write_and_clear( datadir, (n+1)/binsize );
      KCKB.write_and_clear( datadir, (n+1)/binsize );
      KCKC.write_and_clear( datadir, (n+1)/binsize );
      // KCKD.write_and_clear( datadir, (n+1)/binsize );
      // KCKE.write_and_clear( datadir, (n+1)/binsize );
      // KCKF.write_and_clear( datadir, (n+1)/binsize );

      // KDKA.write_and_clear( datadir, (n+1)/binsize );
      // KDKB.write_and_clear( datadir, (n+1)/binsize );
      // KDKC.write_and_clear( datadir, (n+1)/binsize );
      // KDKD.write_and_clear( datadir, (n+1)/binsize );
      // KDKE.write_and_clear( datadir, (n+1)/binsize );
      // KDKF.write_and_clear( datadir, (n+1)/binsize );

      // KEKA.write_and_clear( datadir, (n+1)/binsize );
      // KEKB.write_and_clear( datadir, (n+1)/binsize );
      // KEKC.write_and_clear( datadir, (n+1)/binsize );
      // KEKD.write_and_clear( datadir, (n+1)/binsize );
      // KEKE.write_and_clear( datadir, (n+1)/binsize );
      // KEKF.write_and_clear( datadir, (n+1)/binsize );

      // KFKA.write_and_clear( datadir, (n+1)/binsize );
      // KFKB.write_and_clear( datadir, (n+1)/binsize );
      // KFKC.write_and_clear( datadir, (n+1)/binsize );
      // KFKD.write_and_clear( datadir, (n+1)/binsize );
      // KFKE.write_and_clear( datadir, (n+1)/binsize );
      // KFKF.write_and_clear( datadir, (n+1)/binsize );

      KA_ss.write_and_clear( datadir, (n+1)/binsize );
      KB_ss.write_and_clear( datadir, (n+1)/binsize );
      KC_ss.write_and_clear( datadir, (n+1)/binsize );
      // -------
      std::clog << "iter: " << n+1 << std::endl;

      const int label = (n+1)/binsize;
      const std::string filepath = configdir+"ckpoint"+std::to_string(label)+".dat";
      s.ckpoint( filepath );
      const int label2 = label-1;
      const std::string filepath2 = configdir+"ckpoint"+std::to_string(label2)+".dat";
      std::filesystem::remove( filepath2 );

      const auto end = std::chrono::steady_clock::now();
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << elapsed_seconds.count() << " sec" << std::endl;
      start = std::chrono::steady_clock::now();
    }

  } // end for n

  return 0;
}

