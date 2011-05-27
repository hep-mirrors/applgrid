#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>

#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"


#include "appl_grid/appl_timer.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TPad.h>

#include "LHAPDF/LHAPDF.h"
#include "hoppet_v1.h"

#define DBG true

#include "scales.h"

extern "C" void evolvepdf_(const double& , const double& , double* ); 
extern "C" double alphaspdf_(const double& Q);


static const int nLoops    = 1;
static const int nFlavours = 5;


void GetPdf(const double& x, const double& Q, double* xf) { 
  //  evolvePDF( x, Q, xf);    
  evolvepdf_( x, Q, xf);    
  //hoppeteval_( x, Q, xf);    
  return; 
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

  TString gridName(argv[1]);
  std::cout << "------------------------------------------------"<< std::endl;
  std::cout << "------------------------------------------------"<< std::endl;
  std::cout << "------------------------------------------------"<< std::endl;
          std::cout << "readingfile = " << gridName.Data()<< std::endl;
  std::cout << "------------------------------------------------"<< std::endl;
  std::cout << "------------------------------------------------"<< std::endl;
  std::cout << "------------------------------------------------"<< std::endl;
  
  appl::grid g(gridName.Data());
  g.trim();

  double gridEnergy = g.getCMSScale();
  std::cout << "grid CMS energy " << gridEnergy << std::endl; 


  // get all the reference histograms

  TFile* f;
  f = new TFile(argv[1]);

  TString outFileName = gridName(0,gridName.Index(".root"));
  outFileName += "_standSimple.root";

  TFile* fout = new TFile(outFileName, "recreate");

  Directory ref("reference");
  ref.push();

  TH1D* reference = (TH1D*)f->Get("grid/reference");
  
  TH1D* referenceSubProc[50];

  for ( int i=0 ; i< g.subProcesses() ; i++  ) {
    char hname[64];
    sprintf( hname, "procReference_%d", i); 

    referenceSubProc[i] = (TH1D*)f->Get(hname);
    referenceSubProc[i]->Write();
  }

  ref.pop();


  // now calculate all the cross sections
  //  const string _pdfname = "PDFsets/cteq66.LHgrid";  
  const string _pdfname = "cteq66.LHgrid"; // with new LHAPDF, don't need to specify or link directory   
  int Npdf = 0;
  //  setup gavins code
  //  initPDFset(_pdfname.c_str());
  //  initPDF(Npdf);


  LHAPDF::initPDFSet( _pdfname, Npdf );


  Directory xsDir("xSection");
  xsDir.push();
  
  TH1D* xsec = g.convolute(GetPdf, alphaspdf_, nLoops); 
  xsec->SetName("xsec");
  xsec->SetTitle(reference->GetTitle());

  TH1D* xsec_subproc[50];
  for ( int i=0; i < g.subProcesses() ; i++ ) 
    { 
      char hname[64];
      sprintf( hname, "xsec_subproc_%d", i); 
      
      xsec_subproc[i] = g.convolute_subproc(i,GetPdf, alphaspdf_ , nLoops, mur[i], muf[i]);
      xsec_subproc[i]->SetName(hname);
      xsec_subproc[i]->SetTitle(referenceSubProc[i]->GetTitle());
      xsec_subproc[i]->Write();
    }
  xsec->Write();  
  xsDir.pop();

  // now take all the ratios etc

  Directory ratiodir("ratio");
  ratiodir.push();

  TH1D* ratio = divide( xsec, reference ); 
  if ( ratio ) {
    ratio->SetName("ratio");
    ratio->Write();
  }
  
  TH1D* ratio_subproc[Nscales];

  for ( int i=0; i < Nscales ; i++ ) 
    { 
      
      char hname[64];
      sprintf( hname, "ratio_subproc_%d", i); 
      
      ratio_subproc[i] = divide( xsec_subproc[i], referenceSubProc[i] );
      if ( ratio_subproc[i] ) 
	{
	  ratio_subproc[i]->SetName(hname);
	  ratio_subproc[i]->Write();
	}
    }

  ratiodir.pop();

  fout->Write();
  fout->Close();

  return 0;
}
