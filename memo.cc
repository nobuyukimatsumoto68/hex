
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

//   std::cout << "beta = " << beta << std::endl;
//   std::cout << "eps = " << eps_tot/counter << std::endl;
//   std::cout << "h0 = " << h0_tot/counter << std::endl;
// }



