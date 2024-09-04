#include <iostream>
#include <iomanip>
#include <filesystem>

#include <chrono>
#include <limits>

#include <omp.h>

#include "header2.hpp"


int main( int argc, char *argv[] ){

  int init_label = 0;
  // bool if_read = false;
  bool if_read = true;
  // if (argc>1){
  //   mult = atoi(argv[1]);
  //   Lx = 3*mult; // 12
  //   Ly = 3*mult;
  // }
  // if (argc>2){
  //   if_read=true;
  //   init_label = atoi(argv[2]);
  // }

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
  std::vector<Obs<Scalar>*> scalars;
  std::vector<Obs<Corr>*> corrs;

  Idx x1 = 0, y1 = 0, x2 = Lx/2, y2 = 0;
  Obs<Corr> ss_corr( "ss_corr", binsize, [](const Spin& s0){ return Corr( s0.ss_corr() ); } );
  corrs.push_back( &ss_corr );
  Obs<Scalar> eps_1pt( "eps_1pt", binsize, [](const Spin& s0){ return Scalar( s0.eps_1pt() ); } );
  scalars.push_back( &eps_1pt );
  Obs<Corr> epseps_corr( "epseps_corr", binsize, [](const Spin& s0){ return Corr( s0.epseps_corr() ); } );
  corrs.push_back( &epseps_corr );
  Obs<Scalar> TA( "TA", binsize, [](const Spin& s0){ return Scalar( s0.TM_1pt(0) ); } );
  scalars.push_back( &TA );
  Obs<Scalar> TB( "TB", binsize, [](const Spin& s0){ return Scalar( s0.TM_1pt(1) ); } );
  scalars.push_back( &TB );
  Obs<Scalar> TC( "TC", binsize, [](const Spin& s0){ return Scalar( s0.TM_1pt(2) ); } );
  scalars.push_back( &TC );

  Obs<Corr> TATA( "TATA", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(0,0) ); } );
  corrs.push_back( &TATA );
  Obs<Corr> TATB( "TATB", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(0,1) ); } );
  corrs.push_back( &TATB );
  Obs<Corr> TATC( "TATC", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(0,2) ); } );
  corrs.push_back( &TATC );
  Obs<Corr> TBTA( "TBTA", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(1,0) ); } );
  corrs.push_back( &TBTA );
  Obs<Corr> TBTB( "TBTB", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(1,1) ); } );
  corrs.push_back( &TBTB );
  Obs<Corr> TBTC( "TBTC", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(1,2) ); } );
  corrs.push_back( &TBTC );
  Obs<Corr> TCTA( "TCTA", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(2,0) ); } );
  corrs.push_back( &TCTA );
  Obs<Corr> TCTB( "TCTB", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(2,1) ); } );
  corrs.push_back( &TCTB );
  Obs<Corr> TCTC( "TCTC", binsize, [](const Spin& s0){ return Corr( s0.TMTN_corr(2,2) ); } );
  corrs.push_back( &TCTC );

  Obs<Corr> TA_ss( "TA_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.TM_ss_corr(0, x1, y1, x2, y2) ); } );
  corrs.push_back( &TA_ss );
  Obs<Corr> TB_ss( "TB_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.TM_ss_corr(1, x1, y1, x2, y2) ); } );
  corrs.push_back( &TB_ss );
  Obs<Corr> TC_ss( "TC_ss", binsize, [x1,y1,x2,y2](const Spin& s0){ return Corr( s0.TM_ss_corr(2, x1, y1, x2, y2) ); } );
  corrs.push_back( &TC_ss );
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
    for(auto ptr : scalars) ptr->meas(s);
    for(auto ptr : corrs) ptr->meas(s);

    //------

    // write out and clear
    if( (n+1)%binsize==0 && if_write ) {
      // -----------
      for(auto ptr : scalars) ptr->write_and_clear( datadir, (n+1)/binsize );
      for(auto ptr : corrs) ptr->write_and_clear( datadir, (n+1)/binsize );
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

