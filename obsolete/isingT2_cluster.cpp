// mu =  x,y,z, -x, -y, -z = 0,1,2,3,4,5

/* ---- LATTICE ---- */
#define Lx 4
#define Ly 5
#define Dim 2
#define Three 3
#define N Lx*Ly

/* ---- MONTE CARLO ---- */
#define NConf 10000
#define Interval 2
#define Seed 1

#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>
#include <stack>

double Beta = 0.274653;

void printSpins(const int s[Ly][Lx]);

double heatBath(int s[Ly][Lx]);
double getM(const int s[Ly][Lx]);

// random number
std::mt19937 gen;
std::uniform_int_distribution<> dist01; // [a,b] = [0,2147483647]; we check a%2=0, b%2=1.
// std::uniform_int_distribution<int> dist(0,9);
// std::uniform_real_distribution<double> distribution(0.0,1.0);
int uniformPM1(){ return 2 * (dist(gen)%2) - 1; }

int uniformInt0N(){ return dist(gen)*N/; }


double uniform0to1(){ return (double)dist(gen)/dist.b(); }

int neighbor(int n, int mu);
int FlipCluster( int* s );


int main( int argc, char *argv[] ){
  std::cout << std::scientific << std::setprecision(15);

  if(argc == 2){
    Beta = std::atof(argv[1]);
  }
  else{
    std::cout << "# proceed with the default value." << std::endl; // default nRank is 1
  }


  assert( dist.a()==0 );
  assert( dist.b()%2==1 );
  gen.seed(Seed);

  std::cout <<
    "# Lx: " << Lx << std::endl <<
    "# Ly:  " << Ly << std::endl <<
    "# Seed: " << Seed << std::endl <<
    "# Beta: " << Beta << std::endl <<
    "# NConf: " << NConf << std::endl <<
    "# Interval: " << Interval << std::endl
    ;

  // --------------------------------------

  /* === Site Variables === */
  int s[Ly][Lx];

  // initialize
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      s[y][x] = uniformPM1();
    }}

  printSpins(s);

  // // update & record
  // std::cout << "#" << std::setw(21) << "mag per site"
  //           << std::setw(23) << "acceptance" << std::endl;
  // for(int k=0; k<NConf; k++) {
  //   const double r = heatBath(s, gen);
  //   // std::cout << "# k = " << k << std::endl;
  //   // printSpins(s);

  //   if((k+1)%Interval==0) {
  //     const double b = getM(s);
  //     std::cout << std::setw(22) << b << " "
  //               << std::setw(22) << r << std::endl;
  //   }
  // }

  return 0;
}




int neighbor(int n, int mu)
{
  const int sign = -2*(mu/Three) + 1;
  const int dir = mu % Three;

  int x = n%Lx;
  x += sign*(dir+1)%2;
  x = (x+Lx)%Lx;

  int y = n/Lx;
  y += sign*(dir+1)/2;
  y = (y+Ly)%Ly;

  return Lx*y + x;
}




int FlipCluster(int* s){
  bool isCluster[N];
  for(int n=0; n<N; n++) isCluster[n] = false;
  std::stack<int> newlyAdded;

  // init (a)
  const int init = uniformInt0N(gen);
  int orig_sign = s[init];
  is_cluster[init] = true;
  newlyAdded.push(init);

  while( newlyAdded.size() != 0 ){
    int n = newlyAdded.top();
    newlyAdded.pop();
    wolff_cluster_size++;

    s[n] = -s[n]; // flip when visited

    // try to add neighbors
    for(int mu = 0; mu < Dim*Three; mu++){
      const int nn = neighbor(n,mu);
      if( s[nn] == s[n] || is_cluster[nn] ) continue; // s[x]*sR[y]<0 or y in c
      const double r = uniform01(gen);
      if( r < exp(-2.0 * Beta) ) continue; // reject
      is_cluster[nn] = true;
      newlyAdded.push(nn);
    }
  }

  const int signed_size =  orig_sign * wolff_cluster_size;
  return signed_size;
}







// test code

// int n = Lx*2+0;
// std::cout << "n = " << n << std::endl;
// for(int mu=0; mu<Dim*Three; mu++){
//   std::cout << mu << ", " << neighbor(n,mu) << std::endl;
// }





// int main()
// { 
//   int s[Ly][Lx];
//   double beta; //  Triagular critical point β_crit = ln(3)/4 =  0.274653
//   beta = log(3.0)/4.0;
//   ofstream out_Heat;
//   out_Heat.open("Heat.dat");
//     ofstream out_Wolff;
//   out_Wolff.open("Wolff.dat");
  
//    cout<<" Ising on Lx by Ly  triagular lattice: "<< Lx <<"  " << Ly <<" critical point " << log(3.0)/4.0 << endl;
 
//   srand(137);
//   // Hot start
//     for(int y = 0; y< Ly; y++)
//       for(int x = 0; x < Lx; x++)
// 	s[y][x] = 2 * ( rand() % 2) - 1;

//  #if Debug
//     cout << " Magnetization is " << getMag(s) << endl;
//    printSpins(s);
   
//     testMap(s);
//     cout << "After TestMap Magnetization is " << getMag(s) << endl;
//     printSpins(s);
   
// #endif
   
//     beta = beta + 0.02;

//     double AveMag =0.0;
//     int heatMag;

//    for(int iter = 0;iter < 100000; iter++)
//     {
//        HeatBath(s,beta);
//         heatMag =  getMag(s);
// 	AveMag += (double) abs(heatMag);
       
//     if(iter%1000  == 0 )
// 	{
// 	  cout << iter;
// 	  cout <<" In Heat Bath   " << heatMag <<"   " << endl;
// 	  out_Heat << AveMag/1000.0 << endl;
// 	  AveMag = 0.0;
// 	}     
//     }

//    int WolffMag  =  getMag(s);

//    double AveAbsMag = 0.0;
   
//    for(int iter = 0;iter < 100000 ; iter++)
//      {
//        int  signed_flip =  FlipCluster(s,beta);
//        WolffMag += - 2*signed_flip;
   
//       AveMag += (double) WolffMag;
//        AveAbsMag += abs( (double) WolffMag);
  
//        if(iter%1000  == 0 )
//        {
// 	 cout << iter;
// 	 cout <<" In Wolff signed_flip =   " << signed_flip << "Mag and  Update mag "<< getMag(s) << " " << WolffMag << endl;
// 	 out_Wolff << AveAbsMag/1000.0 <<"  " << AveMag/1000.0 <<endl;
// 	   AveMag = 0.0;
// 	   AveAbsMag = 0.0;   
//        }
//      }
//    return 0;
// }
  
// /*******************************************************
// Ising exp[beta sum_<i,j> s_i s_j]  and h_i = sum_{j nn of i} s_j
// Defining  prob =  exp[2*beta*h_i]
// prob to aligning with h_i
// is  exp[beta*h]/(exp[beta*h] + (exp[-beta*h])   = pup/(1 + pup)
// *******************************************************/


// void  HeatBath(int s[Ly][Lx],double beta)
// {
//   int h = 0;
//   int xp,xm,yp, ym;
//   double xran, prob;
  
//       for(int y = 0; y<Ly; y++)
// 	for(int x= 0 ; x<Lx; x++)
// 	  { xp = mod(x+1,Lx); xm =  mod(x-1,Lx); yp  = mod(y+1,Ly); ym =  mod(y-1,Ly);
// 	    h = s[y][xp] +  s[yp][x] +  s[yp][xp] +  s[y][xm] +  s[ym][x] +  s[ym][xm];
// 	    xran = (double)rand()/(double)RAND_MAX;
// 	    prob = exp(2.0*beta* h); 
// 	    s[y][x] = xran < prob/(1.0 + prob) ? +1 : -1;
// 	  } 
// }


// int  FlipCluster( int s[Ly][Lx], double beta){
//   int *spt;
//   spt = *s;  // external field
//   bool is_cluster[N];
//   int  wolff_cluster_size = 0;
//   int mag_sign;
//   double K[2 *Three];
//   for(int mu =0;mu < 2 *Three; mu++) K[mu] = beta;

//    // index to site n =  x * Lx*y so that value is spt[n] == s[y][x];
//   for(int n= 0; n<N; n++) is_cluster[n] = false;  // clear is_cluster
  
//   // create the stack
//   std::stack<int> stack;
//   // choose a random site and add it to the cluster
//   int site =    rand() % (N); // careful with defines
//    mag_sign = spt[site];
//    is_cluster[site] = true;
//    stack.push(site);
//    //  cout << " pushed in stack site = " << site << endl;
   
  
//     //     cout << "Seed Site = " << site << endl;

 
//   while (stack.size() != 0) {
//     int  n = stack.top();
//     stack.pop();
//       wolff_cluster_size += 1;

//     // flip the spin
//       spt[n] = - spt[n];

//     // try to add neighbors
//     for (int mu = 0; mu < 2*Three; mu++) {
//       int nn = neighbor(n,mu);
//       // skip if the site not allinged or  already pushed into clustered
//        if ( spt[nn] == spt[n] || is_cluster[nn]  )    continue;
//        double rate = -2.0 * K[mu]; // link coupling
//        double xran = (double)rand()/(double)RAND_MAX;	    
//       if (rate >= 0.0 || xran < exp(rate)) continue;
//       is_cluster[nn] = true;
//       stack.push(nn);
//       //     cout << " pushed in stack nn = " << nn << endl;
//     }
//   }

//    int  signed_flip =  mag_sign*wolff_cluster_size;
    
//   return signed_flip; // signed cluster
// }


// int  getMag( int s[Ly][Lx])
// {
//  int  mag = 0;
//   for(int y = 0; y<Ly; y++)
//     for(int x= 0 ; x<Lx; x++)
//       mag += s[y][x];
//   return mag;
// }

// void printSpins(int s[Ly][Lx])
// {
// 	cout<<"\n-------------------------------------------- \n";

// 	for(int y = 0; y<Ly; y++)
// 	  {
// 	    for(int x= 0 ; x<Lx; x++)
// 	      cout<< s[y][x] <<"   ";
// 	    cout << endl; 
// 	  }					 
// } 
  
//  void testMap(int s[Ly][Lx])
// {
//   int *spt;
//   spt = *s;
//   cout<<"\n Insidet testMap \n";
// 	for(int y = 0; y<Ly; y++)
// 	  {
// 	    for(int x= 0 ; x<Lx; x++)
// 	      cout<< spt[ x + y*Lx] <<"   ";
// 	    cout << endl; 
// 	  }

// 	cout<<"\n Remap insdie estMap \n";
//       	for(int y = 0; y<Ly; y++)
// 	  {
// 	    for(int x= 0 ; x<Lx; x++){
// 	      spt[x + y*Lx] = x + y*Lx;
// 	       cout << spt[x + y*Lx] << "   ";
// 	    }
// 	     cout << endl;
// 	  }
		
// }
 
 
