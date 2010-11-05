#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"

// #include "appl_grid/appl_pdf.h"


#include "LHAPDF/LHAPDF.h"
extern "C" void evolvepdf_(const double& , const double& , double* ); 
extern "C" double alphaspdf_(const double& Q);

#include "appl_grid/appl_timer.h"
#include "hoppet_v1.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TPad.h>


extern "C" 
{
  
  //void dglapeval_(const double& _x, const double& _Q, double* f);
  //void dglapevalsplit_(const double& _x, const double& _Q, const int&, const int&, double* f);
  void initmypdf_(const char* name, const int& set);
  
  double alphaspdf_(const double& Q);
} 




// wrapper to get the x*pdf from gavins evolution code.
// in fact, don't actually need this wrapper any longer - we could
// just pass the routine directly

void GetPdf(const double& x, const double& Q, double* f) { 

  hoppeteval_( x, Q, f);    

  return; 
}





TH1D* divide( const TH1D* h1, const TH1D* h2 ) {

  if ( h1==NULL || h2==NULL ) return NULL;
 
  TH1D* h = (TH1D*)h1->Clone();

  std::cout << "histograms " << h1->GetTitle() << " " << h2->GetTitle() << std::endl;

  

  for ( int i=1 ; i<=h1->GetNbinsX() ; i++ ) { 
    double b  = h2->GetBinContent(i);
    double be = h2->GetBinError(i);
    double t  = h1->GetBinContent(i);
    double te = h1->GetBinError(i);

    double r  = ( b!=0 ? t/b : 0 );
    //    double re = ( b!=0 ? sqrt((r+1)*r/b) : 0 );
    double re = 0;
    
    h->SetBinContent( i, r );
    h->SetBinError( i, re ) ;

    //    if ( DBG ) std::cout << "\tx=" << h->GetBinCenter(i) << "\tratio=" << r << std::endl;
  } 

  double hmin = h->GetBinContent(1);
  double hmax = h->GetBinContent(1);
  
  for ( int i=2 ; i<=h->GetNbinsX() ; i++ ) { 
    double d = h->GetBinContent(i);
    if ( hmin>d ) hmin=d;
    if ( hmax<d ) hmax=d; 
  }

  std::cout << "\tmin ratio = " << hmin << "\tmax ratio = " << hmax << std::endl;
  
  if ( h->GetMaximum()<1.01 ) h->SetMaximum(1.01);
  if ( h->GetMinimum()>0.99 ) h->SetMinimum(0.99);

  return h;
}




int main(int argc, char** argv) { 

  if ( argc<2 ) return -1;

  appl::grid g(argv[1]);
  g.trim();

  // get grid cms energy

  double gridEnergy = g.getCMSScale();

  std::cout << "grid CMS energy " << gridEnergy << std::endl; 


  // get all the reference histograms

  TFile* f;
  if ( argc>2 ) f = new TFile(argv[2]);
  else          f = new TFile(argv[1]);

  TFile* fout = new TFile("xsec.root", "recreate");
  Directory ref("reference");
  ref.push();



  TH1D* reference = (TH1D*)f->Get("grid/reference"); reference->Write();


  ref.pop();


  
  // now calculate all the cross sections
  
  // setup lhapdf etc
  // const string _pdfname = "lhapdf/PDFsets/cteq6mE.LHgrid";  
  const string _pdfname = "PDFsets/cteq6mE.LHgrid";  
  int Npdf = 0;
  // setup gavins code
  initmypdf_(_pdfname.c_str(), Npdf);

  // const string _pdfname = "cteq6mE.LHgrid";  
  // LHAPDF::initPDFSet( _pdfname, Npdf );

  bool first = true;

  TCanvas* ratioc = new TCanvas("ratio",     "ratio",     500, 500);
  TCanvas* refc   = new TCanvas("reference", "reference", 500, 500);

  Directory xsecdir("xsec");
  xsecdir.push();


 
  //  TH1D* xsec = g.convolute( GetPdf, alphaspdf_ , nLoops ); xsec->SetName("xsec");
  TH1D* xsec = g.convolute( GetPdf, alphaspdf_ ); 
  xsec->SetName("xsec");
  xsec->SetTitle(reference->GetTitle());
  
  xsecdir.pop();

  // now take all the ratios etc

  Directory ratiodir("ratio");
  ratiodir.push();

  TH1D* ratio = divide( xsec, reference ); if ( ratio ) ratio->SetName("ratio");
  

  // some test output 
 
  ratioc->cd();
  if ( first ) ratio->DrawCopy();
  else         ratio->DrawCopy("same");

  refc->cd();
  if ( first ) reference->DrawCopy();
  else         reference->DrawCopy("same");

  first = false;

  // done
  
  ratiodir.pop();


  ratioc->cd();
  gPad->Print("ratio.gif");

  refc->cd();
  gPad->SetLogy();
  gPad->Print("reference.gif");

 
  fout->Write(); 

  return 0;
}
