#include <iostream>
#include <string>
#include <vector>
// #include <stdio.h>

#include "appl_grid/appl_grid.h"

#include "TH1D.h"


// lhapdf routines
extern "C" 
{  
  void initmypdf_(const char* name, const int& iset);
  //  void initpdfset_( const char* );
  //  void initpdf_( const int& );

  void evolvepdf_(const double& , const double& , double* );
  double alphaspdf_(const double& Q);
} 


// wrapper to get the x*pdf from gavins evolution code.
// in fact, don't actually need this wrapper any longer - we could
// just pass the routine directly

// void GetPdf(const double& x, const double& Q, double* xf) { 
//   evolvepdf_( x, Q, xf);    
//   return; 
// }



int main(int argc, char** argv) { 

  if ( argc<2 ) return -1;

  std::cout << "reading grid" << std::endl; 

  // get name of grid from user and create from grid file
  appl::grid g(argv[1]);
  g.trim(); // trim away uneeded memory

  std::cout << "setting up lhapdf" << std::endl; 

  // initialise lhapdf

  const string _pdfname = "PDFsets/cteq6mE.LHgrid";  
  int iset = 0;

  // need a fortran wrapper to do the fortran string passing 
  // properly
  initmypdf_(_pdfname.c_str(), iset);
  // initpdfset_(_pdfname.c_str());
  // initpdf_( iset );

  // do the convolution into a vector
  std::cout << "doing standalone convolution..." << std::endl; 
  std::vector<double>  xsec = g.vconvolute( evolvepdf_, alphaspdf_ ); 

  // or get into a histogram
  //  TH1D*               hxsec = g.convolute( GetPdf, alphaspdf_ ); 
  //  hxsec->SetName("xsec");

  return 0;
}