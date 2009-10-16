

#include "appl_grid/fastnlo.h"
#include "TFile.h"
#include "TH1D.h"

#include <iostream>
#include <stdio.h>


// lhapdf routines
extern "C" 
{  
  void initmypdf_(const char* name, const int& iset);
  //  void initpdfset_( const char* );
  //  void initpdf_( const int& );

  void evolvepdf_(const double& , const double& , double* );
  double alphaspdf_(const double& Q);
} 




int main(int argc, char** argv) { 

  if ( argc<2 ) { 
    std::cerr << "usage: fnlo fastnlogrid.tab" << std::endl;
    return -1;
  }

  std::cout << "set up lhapdf..." << std::endl;

  const std::string _pdfname = "PDFsets/cteq61.LHgrid";  
  int iset = 0;
  initmypdf_(_pdfname.c_str(), iset);
  
  // don't need to hard code, can use a runtime parameter...
  std::string gridname = argv[1];
  // std::string gridname = "fnt1007midp.tab";
  // std::string gridname = "fnt1008midp.tab";
  // std::string gridname = "fnt2004.tab";
  // std::string gridname = "fnl2004.tab";

  fastnlo f( gridname );

  std::vector<appl::grid*> g = f.grids(); 

  std::vector<TH1D*> hc(g.size());
  
  std::string foutname = "appl.root";

  TFile fout( foutname.c_str(),"recreate");
  for ( int i=0 ; i<g.size() ; i++ ) { 

    // trim the grids (not actually needed, done in 
    // the fastnlo constructor
    g[i]->trim();

    char hname[64];
    sprintf(hname, "hist%02d", i);

    hc[i] = g[i]->convolute( evolvepdf_, alphaspdf_ );
    hc[i]->SetName(hname);
    hc[i]->SetDirectory(&fout);
    hc[i]->Write();

    std::cout << "central" << std::endl;
    for ( int j=0 ; j<hc[i]->GetNbinsX() ; j++ ) { 
      std::cout << "xsec " << j << "\t" << hc[i]->GetBinContent(j+1) << std::endl; 
    } 

    std::cout << "scale " << 0.5 << std::endl;
    
    std::string hsname = hname;
    hsname += "-0.5";

    hc[i] = g[i]->convolute( evolvepdf_, alphaspdf_, g[i]->nloops(), 0.5, 0.5 );
    hc[i]->SetName(hsname.c_str());
    hc[i]->SetDirectory(&fout);
    hc[i]->Write();

    for ( int j=0 ; j<hc[i]->GetNbinsX() ; j++ ) { 
      std::cout << "xsec " << j << "\t" << hc[i]->GetBinContent(j+1) << std::endl; 
    } 

    std::cout << "scale " << 0.5 << std::endl;
    
    hsname = hname;
    hsname += "-2.0";

    hc[i] = g[i]->convolute( evolvepdf_, alphaspdf_, g[i]->nloops(), 2.0, 2.0 );
    hc[i]->SetName(hsname.c_str());
    hc[i]->SetDirectory(&fout);
    hc[i]->Write();

    for ( int j=0 ; j<hc[i]->GetNbinsX() ; j++ ) { 
      std::cout << "xsec " << j << "\t" << hc[i]->GetBinContent(j+1) << std::endl; 
    } 

  }

  std::cout << "writing file " << foutname << std::endl; 
  fout.Close();
  
  return 0;
}
