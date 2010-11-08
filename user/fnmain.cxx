

#include <iostream>
#include <stdio.h>

#include "appl_grid/fastnlo.h"

#include "TFile.h"
#include "TH1D.h"


// lhapdf routines
#include "LHAPDF/LHAPDF.h"
extern "C" void evolvepdf_(const double& , const double& , double* );
extern "C" double alphaspdf_(const double& Q);



int main(int argc, char** argv) { 

  if ( argc<2 ) { 
    std::cerr << "usage: fnlo fastnlogrid.tab" << std::endl;
    return -1;
  }

  std::cout << "set up lhapdf..." << std::endl;

  const std::string _pdfname = "cteq6mE.LHgrid";  
  int iset = 0;
  LHAPDF::initPDFSet( _pdfname, iset );
  
  // don't need to hard code, can use a runtime parameter...
  std::string gridname = argv[1];
  // std::string gridname = "fnt1007midp.tab";
  // std::string gridname = "fnt1008midp.tab";
  // std::string gridname = "fnt2004.tab";
  // std::string gridname = "fnl2004.tab";

  fastnlo f( gridname );

  std::vector<appl::grid*> g = f.grids(); 

  //  g.push_back( new appl::grid("atlas-incljets04-eta1.root") );

  /// histograms and file for seeing the results
  std::vector<TH1D*> hc(g.size());
  std::string foutname = "appl.root";
  TFile fout( foutname.c_str(),"recreate");

  for ( int i=0 ; i<g.size() ; i++ ) { 

    // trim the grids (not actually needed, done in 
    // the fastnlo constructor
    g[i]->trim();

    char hname[64];
    sprintf(hname, "hist%02d", i);

    /// optionally print out the grid documentation  
    std::cout << "\n" << g[i]->getDocumentation() << std::endl;

    /// perform the convolution
    hc[i] = g[i]->convolute( evolvepdf_, alphaspdf_ );
    hc[i]->SetName(hname);
    hc[i]->SetDirectory(&fout);
    hc[i]->Write();
    
    /// print out the results
    for ( int j=1 ; j<=hc[i]->GetNbinsX() ; j++ ) { 
      std::cout << "xsec(" << j-1 << ")=" << hc[i]->GetBinContent(j) << std::endl;
    }

  }

  std::cout << "writing file " << foutname << std::endl; 
  fout.Close();
  
  return 0;
}

