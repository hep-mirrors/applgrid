#include <iostream>
#include <string.h>

using std::cerr;
using std::cout;
using std::endl;

#include "appl_grid/appl_grid.h"
using appl::grid;

#include "appl_grid/appl_pdf.h"
using appl::appl_pdf;

#include "appl_grid/appl_timer.h"

#include <TH1D.h>
#include <TFile.h>
//
//
//

extern "C" 
{
  
  void evolvepdf_(const double& , const double& , double* );
  void initmypdf_(const char* name, const int& set);
  
  double alphaspdf_(const double& Q);
} 

//
//
//
//static const double pb_fac = 3.89379656e8 ;    // conversion GeV^2 -> pb  

static const int nScales = 5;

static const double mur[nScales] = {1.0, 0.5, 2.0, 1.0, 0.5};
static const double muf[nScales] = {1.0, 0.5, 2.0, 0.5, 1.0};
//static const double mur[nScales] = {1.00, 0.50, 1.00};
//static const double muf[nScales] = {1.00, 1.00, 0.50};
//static const double mur[nScales] = {1.00, 0.50, 1.00, 2.00, 1.00, 2.00};
//static const double muf[nScales] = {1.00, 1.00, 0.50, 1.00, 2.00, 2.00};
static const int nLoops = 1;
static const int nFlavours = 5;


TH1D* reference;
TH1D* gridObservable;
TH1D* ratio;

TH1D* soft_nlo;
TH1D* grid_nlo;
TH1D* ratio_nlo;

TH1D* soft_lo;
TH1D* grid_lo;
TH1D* ratio_lo;

//
//
void removeErrors(TH1D* hh)
{
  int nBins = hh->GetNbinsX();
  for (int iBin = 0; iBin <= nBins + 1; iBin++)
    {
      hh->SetBinError(iBin,0.);
    }
}



// wrapper to get the basic pdf rather than x*pdf
// gavins evolution code


void GetPdf(const double& x, const double& Q, double* xf) { 
  evolvepdf_( x, Q, xf);    
  return; 
}




int main(int argc, char** argv) 
{ 
  
  cout << "main()" << endl;
  string gridFilePrefix = "./output/";
  //string gridFilePrefix = "./parameters/fig6/";
  string gridFileName = gridFilePrefix; 
  string  outFileName = gridFilePrefix + "_Standalone_"; 
  
  if ( argc<2 ) 
    {
      cerr << "weight file not specified" << endl; 
      return -1;
    }
  else
    {
      gridFileName += argv[1];
      outFileName += argv[1];
    }
  
  gridFileName += ".root";
  outFileName += ".root";
  if (FILE* testfile = fopen(gridFileName.c_str(),"r"))
    {
      fclose(testfile);
    }
  else
    {
      cout <<"No such ROOT-file "<<gridFileName<<endl;
      return -1;
    }
  
//  const string _pdfname = "PDFsets/cteq66.LHgrid";  
  const string _pdfname = "PDFsets/cteq6mE.LHgrid";  
  int Npdf = 0;
  // setup gavins code
  initmypdf_(_pdfname.c_str(), Npdf);
  cout << "lhapdf..." << endl;
  //
  //
  // read in the grid
  //
  //
  TFile outfile(outFileName.c_str(), "recreate");
  TFile infile(gridFileName.c_str(), "read");
  
  grid* g = new grid(gridFileName);
  
  //g->print();
  
  g->trim();
  
  reference = (TH1D*)infile.Get("grid/reference");
  //
  // observable
  //  
  cout << "convoluting observable" << endl;

  gridObservable = (TH1D*)g->convolute(GetPdf, alphaspdf_);

  gridObservable->SetName("gridObservable");
  gridObservable->SetTitle("Observable calculated via APPLgrid");
  
  ratio = (TH1D*)reference->Clone("ratio");
  ratio->Divide(gridObservable);
  ratio->SetTitle("reference/grid ratio");

  outfile.cd();
  reference->Write();
  gridObservable->Write();
  removeErrors(ratio);
  ratio->Write();
  //
  //   writing results
  //
  //  outfile.Write();
  outfile.Close();
  infile.Close();
  
  delete g;
  return 0;
}








