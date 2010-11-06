#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>

#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"

// #include "appl_grid/appl_pdf.h"

#include "LHAPDF/LHAPDF.h"
#include "appl_grid/appl_timer.h"
// #include "hoppet_v1.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TPad.h>


extern "C" 
{
  
  //void dglapeval_(const double& _x, const double& _Q, double* f);
  //void dglapevalsplit_(const double& _x, const double& _Q, const int&, const int&, double* f);
  void evolvepdf_(const double& , const double& , double* );
  //  void initmypdf_(const char* name, const int& set);
  
  double alphaspdf_(const double& Q);
} 


#define DBG true


// static const double pb_fac = 3.89379656e8 ;    // conversion GeV^2 -> pb  

#include "scales.h"

static const int nLoops    = 1;
static const int nFlavours = 5;



// wrapper to get the x*pdf from gavins evolution code.
// in fact, don't actually need this wrapper any longer - we could
// just pass the routine directly

void GetPdf(const double& x, const double& Q, double* xf) { 
  //  double xf[13];
  //  hoppeteval_( x, Q, xf);    
  //  hoppeteval_( x, Q, xf);    
  evolvepdf_( x, Q, xf);    
  //if (debug) cout << "\t evo=" << xf[6];
  //if (debug) cout << " x= "<<" Q= "<<Q<<"\tdgl=" << xf[6] << endl;
  //  double invx=0.;
  //  if (x!=0.) invx=1./x;
  //  for ( int i=0; i<13 ; i++ ) f[i] = xf[i]*invx;
  return; 
}




void increment( TH1D* h1, const TH1D* h2 ) {

  if ( h1==NULL || h2==NULL ) return;
 
  if ( DBG ) std::cout << "increment histograms " << h1->GetTitle() << " " << h2->GetTitle() << std::endl;

  for ( int i=1 ; i<=h1->GetNbinsX() ; i++ ) { 
    double b1 = h1->GetBinContent(i);
    double b2 = h2->GetBinContent(i);

    h1->SetBinContent(i, b1+b2); 
  }

}



TH1D* divide( const TH1D* h1, const TH1D* h2 ) {

  if ( h1==NULL || h2==NULL ) return NULL;
 
  TH1D* h = (TH1D*)h1->Clone();

  if ( DBG ) std::cout << "histograms " << h1->GetTitle() << " " << h2->GetTitle() << std::endl;

  

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

  if ( DBG ) std::cout << "\tmin ratio = " << hmin << "\tmax ratio = " << hmax << std::endl;
  
  if ( h->GetMaximum()<1.01 ) h->SetMaximum(1.01);
  if ( h->GetMinimum()>0.99 ) h->SetMinimum(0.99);

  return h;
}



// extern double Escale;

int main(int argc, char** argv) { 

  if ( argc<2 ) return -1;

  appl::grid g(argv[1]);
  g.trim();
  //  g *= pb_fac;

  // get grid cms energy

  double gridEnergy = g.getCMSScale();

  std::cout << "grid CMS energy " << gridEnergy << std::endl; 

  // get energy scale argument if needed

  double Escale = 1;
  if ( argc>3 ) Escale = std::atof(argv[3]);


  // get all the reference histograms

  TFile* f;
  if ( argc>2 ) f = new TFile(argv[2]);
  else          f = new TFile(argv[1]);
  TFile* fout = new TFile("xsec.root", "recreate");
  Directory ref("reference");
  ref.push();



  TH1D* reference = (TH1D*)f->Get("grid/reference"); reference->Write();

  TH1D* soft_scale[Nscales];

  std::vector<TH1D*> soft_sub;
  std::vector<TH1D*> soft_subscale[Nscales];

  // get number of sub processes

  int Nsub = g.subProcesses();

  // get all the subprocess histos for the overall cross section
  
  for ( int i=0 ; i<Nsub ; i++  ) {
    char hname[64];
    sprintf( hname, "addReference/soft_sub_%d", i); 
    soft_sub.push_back( (TH1D*)f->Get(hname) );
    soft_sub.back()->Write();
  }


 
  // now get histo for scale in the scale variation.
  // NB: should really write this sort of information into the file
  //     so we can decide, at run time, how many scales there are, 
  //     and what values they take. Some TVectors would be ideal, 
  //     and all these reference files should really be in their
  //     own file, and the grid kept by itself in it's own special
  //     grid file

  for ( int i=Nscales ; i-- ; ) { 

    // get histo for this scale 
    
    char hname[64];
    sprintf( hname, "addReference/soft_scale_%d", i); 

    soft_scale[i] = (TH1D*)f->Get(hname);
    soft_scale[i]->Write();

  }

  ref.pop();


  
  // now calculate all the cross sections
  
  // setup lhapdf etc
  // const string _pdfname = "lhapdf/PDFsets/cteq6mE.LHgrid";  
  //  const string _pdfname = "PDFsets/cteq6mE.LHgrid";  
  const string _pdfname = "cteq6mE.LHgrid";  
  int Npdf = 0;
  // setup gavins code
  LHAPDF::initPDFSet(_pdfname.c_str(), Npdf);


  bool first = true;


  TCanvas* ratioc = new TCanvas("ratio",     "ratio",     500, 500);
  TCanvas* refc   = new TCanvas("reference", "reference", 500, 500);

  Directory xsecdir("xsec");
  xsecdir.push();


  

  //  TH1D* xsec = g.convolute( GetPdf, alphaspdf_ , nLoops ); xsec->SetName("xsec");
  TH1D* xsec = g.convolute( Escale, GetPdf, alphaspdf_ ); 
  xsec->SetName("xsec");
  xsec->SetTitle(reference->GetTitle());

  TH1D* xsec_scale[Nscales];

  std::vector<TH1D*> xsec_sub;
  std::vector<TH1D*> xsec_subscale[Nscales];

  // subprocess histos
  
  for ( int i=0 ; i<Nsub ; i++ ) {
    char hname[64];
    sprintf( hname, "xsec_sub_%d", i); 
    //    xsec_sub.push_back(  g.convolute_subproc( i, GetPdf, alphaspdf_ , nLoops) );
    xsec_sub.push_back(  g.convolute_subproc( i, Escale, GetPdf, alphaspdf_ ) );
    xsec_sub.back()->SetName(hname);
    xsec_sub.back()->SetTitle(soft_sub[i]->GetTitle());
  }


  for ( int i=Nscales ; i-- ; ) { 
    
    // get histo for this scale 
    
    char hname[64];
    sprintf( hname, "xsec_scale_%d", i); 

    xsec_scale[i] = g.convolute( Escale, GetPdf, alphaspdf_ , nLoops, mur[i], muf[i] );
    xsec_scale[i]->SetName(hname);
    xsec_scale[i]->SetTitle(soft_scale[i]->GetTitle());

  }
  
  xsecdir.pop();




  // now take all the ratios etc

  Directory ratiodir("ratio");
  ratiodir.push();

  TH1D* ratio = divide( xsec, reference ); if ( ratio ) ratio->SetName("ratio");
  
  TH1D* ratio_scale[Nscales];

  std::vector<TH1D*> ratio_sub;
  std::vector<TH1D*> ratio_subscale[Nscales];
  
  for ( int i=0 ; i<Nsub ; i++ ) {
    char hname[64];
    sprintf( hname, "ratio_sub_%d", i); 
    ratio_sub.push_back(  divide( xsec_sub[i], soft_sub[i]) );
    if ( ratio_sub.back() ) ratio_sub.back()->SetName(hname);
  }

  for ( int i=Nscales ; i-- ; ) { 
    
    // get histo for this scale 
    
    char hname[64];
    sprintf( hname, "ratio_scale_%d", i); 

    ratio_scale[i] = divide( xsec_scale[i], soft_scale[i] );
    if ( ratio_scale[i] ) ratio_scale[i]->SetName(hname);

  }

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
