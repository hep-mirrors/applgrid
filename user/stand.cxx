#include <iostream>
#include <string>
#include <vector>

#include "appl_grid/appl_grid.h"

#include "TH1D.h"


// lhapdf routines
#include "LHAPDF/LHAPDF.h"
extern "C" void evolvepdf_(const double& , const double& , double* ); 
extern "C" double alphaspdf_(const double& Q);



int main(int argc, char** argv) { 

  // use a default atlas inclusive grid
  std::string gridname = "atlas-incljets04-eta1.root";

  if ( argc>1 ) gridname = argv[1];

  std::cout << "reading grid " << gridname << std::endl; 

  // get name of grid from user and create from grid file
  appl::grid g(gridname);

  g.trim(); // trim away uneeded memory

  std::cout << "setting up lhapdf" << std::endl; 

  // initialise lhapdf

  const string _pdfname = "cteq6mE.LHgrid";  
  int iset = 0;

  LHAPDF::initPDFSet( _pdfname, iset );
  // initpdfset_(_pdfname.c_str());
  // initpdf_( iset );

  // do the convolution into a vector
  std::cout << "doing standalone convolution" << std::endl; 
  std::cout << g.getDocumentation() << std::endl; 
  std::vector<double>  xsec = g.vconvolute( evolvepdf_, alphaspdf_ ); 

  // print results
  for ( int i=0 ; i<xsec.size() ; i++ ) { 
    std::cout << "xsec(" << i << ")=" << xsec[i] << std::endl;
  }

  //  or get into a histogram
  //  TH1D* hxsec = g.convolute( evolvepdf_, alphaspdf_ ); 
  //  hxsec->SetName("xsec");

  return 0;
}
