#include <iostream>
#include <string>
#include <vector>
// #include <stdio.h>

#include "appl_grid/appl_grid.h"

#include "TH1D.h"


// lhapdf routines
#include "LHAPDF/LHAPDF.h"
extern "C" void evolvepdf_(const double& , const double& , double* ); 
extern "C" double alphaspdf_(const double& Q);



int main(int argc, char** argv) { 

  if ( argc<2 ) return -1;

  std::cout << "reading grid" << std::endl; 

  

  // get name of grid from user and create from grid file
  appl::grid g(argv[1]);

  for ( int i=2 ; i<argc ; i++ ) {
    //  std::cout << "adding grid " << argv[i] << std::endl; 
    g += appl::grid(argv[i]);
  }


  g.trim(); // trim away uneeded memory

  std::cout << "setting up lhapdf" << std::endl; 

  // initialise lhapdf

  const string _pdfname = "cteq6mE.LHgrid";  
  int iset = 0;

  LHAPDF::initPDFSet( _pdfname, iset );
  // initpdfset_(_pdfname.c_str());
  // initpdf_( iset );

  // do the convolution into a vector
  std::cout << "doing standalone convolution..." << std::endl; 
  std::vector<double>  xsec = g.vconvolute( evolvepdf_, alphaspdf_ ); 

  //  or get into a histogram
  //  TH1D* hxsec = g.convolute( evolvepdf_, alphaspdf_ ); 
  //  hxsec->SetName("xsec");

  return 0;
}
