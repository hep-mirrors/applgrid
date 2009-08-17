#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"

// #include "appl_grid/appl_pdf.h"

#include "appl_grid/appl_timer.h"

#include <TH1D.h>
#include <TFile.h>

extern "C" 
{
  
  void dglapeval_(const double& _x, const double& _Q, double* f);
  void dglapevalsplit_(const double& _x, const double& _Q, const int&, const int&, double* f);
  void initmypdf_(const char* name, const int& set);
  
  double alphaspdf_(const double& Q);
} 


#define DBG true


static const double pb_fac = 3.89379656e8 ;    // conversion GeV^2 -> pb  
static const int nScales = 5;
// static const double mur[nScales] = {1.0, 0.5, 2.00, 1.0, 0.5};
// static const double muf[nScales] = {1.0, 0.5, 2.00, 0.5, 1.0};

// static const double mur[nScales] = { 1.0, 0.5, 2.0, 1.0, 1.0 };
// static const double muf[nScales] = { 1.0, 1.0, 1.0, 0.5, 2.0 };

static const double mur[nScales] = { 1.0, 0.5, 2.0, 1.0, 0.5 };
static const double muf[nScales] = { 1.0, 0.5, 2.0, 0.5, 1.0 };

static const int nLoops    = 1;
static const int nFlavours = 5;



// wrapper to get the basic pdf rather than x*pdf
// gavins evolution code

void GetPdf(const double& x, const double& Q, double* f) { 
  
  double xf[13];

  dglapeval_( x, Q, xf);    
  //if (debug) cout << "\t evo=" << xf[6];
  //if (debug) cout << " x= "<<" Q= "<<Q<<"\tdgl=" << xf[6] << endl;
 double invx=0.;
 if (x!=0.) invx=1./x;
 for ( int i=0; i<13 ; i++ ) f[i] = xf[i]*invx;
 return; 
}

void GetPdfSplit(const double& x, const double& Q, double* f) { 
  
 const bool debug=false;
 if (debug)cout<<" x= "<<x<<" Q= "<<Q<<endl;
 double xf[13];
 dglapevalsplit_( x, Q, nLoops, nFlavours, xf); 
 double invx=0.;
 if (x!=0.) invx=1./x;
 for (int i=0; i<13 ; i++ ) f[i] = xf[i]*invx;
 if (debug){
   for (int i=0; i<13 ; i++ )
   cout<<i<<" f= "<<f[i]<<endl;
 }
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

    if ( DBG ) std::cout << "\tx=" << h->GetBinCenter(i) << "\tratio=" << r << std::endl;
  } 
  
  if ( h->GetMaximum()<1.1 ) h->SetMaximum(1.1);
  if ( h->GetMinimum()>0.9 ) h->SetMinimum(0.9);

  return h;
}


int main(int argc, char** argv) { 

  if ( argc<2 ) return -1;

  appl::grid g(argv[1]);
  g.trim();

  // get all the reference histograms

  TFile* f = new TFile(argv[1]);
  TFile* fout = new TFile("xsec.root", "recreate");
  Directory ref("reference");
  ref.push();

  TH1D* reference = (TH1D*)f->Get("grid/reference"); reference->Write();

  TH1D* soft_scale[nScales];

  std::vector<TH1D*> soft_sub;
  std::vector<TH1D*> soft_subscale[nScales];

  // get number of sub proc

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

  for ( int i=nScales ; i-- ; ) { 

    // get histo for this scale 
    
    char hname[64];
    sprintf( hname, "addReference/soft_scale_%d", i); 

    soft_scale[i] = (TH1D*)f->Get(hname);
    soft_scale[i]->Write();

    // get all the subprocess histos for this scale

    for ( int j=0 ; j<Nsub ; j++ ) {
      char hsname[64];
      sprintf( hsname, "addReference/soft_subscale_%d_%d", j, i); 
      soft_subscale[i].push_back( (TH1D*)f->Get(hsname) );
      soft_subscale[i].back()->Write();
    }
    
  }

  ref.pop();


  
  // now calculate all the cross sections
  
  // setup lhapdf etc
  const string _pdfname = "lhapdf/PDFsets/cteq6mE.LHgrid";  
  int Npdf = 0;
  // setup gavins code
  initmypdf_(_pdfname.c_str(), Npdf);


  Directory xsecdir("xsec");
  xsecdir.push();

  g *= pb_fac;

  //  TH1D* xsec = g.convolute( GetPdf, alphaspdf_ , nLoops ); xsec->SetName("xsec");
  TH1D* xsec = g.convolute( GetPdf, alphaspdf_ ); 
  xsec->SetName("xsec");
  xsec->SetTitle(reference->GetTitle());

  TH1D* xsec_scale[nScales];

  std::vector<TH1D*> xsec_sub;
  std::vector<TH1D*> xsec_subscale[nScales];
  
  for ( int i=0 ; i<Nsub ; i++ ) {
    char hname[64];
    sprintf( hname, "xsec_sub_%d", i); 
    //    xsec_sub.push_back(  g.convolute_subproc( i, GetPdf, alphaspdf_ , nLoops) );
    xsec_sub.push_back(  g.convolute_subproc( i, GetPdf, alphaspdf_ ) );
    xsec_sub.back()->SetName(hname);
    xsec_sub.back()->SetTitle(soft_sub[i]->GetTitle());
  }

  for ( int i=nScales ; i-- ; ) { 
    
    // get histo for this scale 
    
    char hname[64];
    sprintf( hname, "xsec_scale_%d", i); 

    xsec_scale[i] = g.convolute( GetPdf, alphaspdf_ , nLoops, mur[i], muf[i], GetPdfSplit );
    xsec_scale[i]->SetName(hname);
    xsec_scale[i]->SetTitle(soft_scale[i]->GetTitle());

    // get all the subprocess histos for this scale

    for ( int j=0 ; j<Nsub ; j++ ) {
      char hsname[64];
      sprintf( hsname, "xsec_subscale_%d_%d", j, i); 
      xsec_subscale[i].push_back( g.convolute_subproc( i, GetPdf, alphaspdf_ , nLoops, mur[i], muf[i], GetPdfSplit  ) );
      xsec_subscale[i].back()->SetName(hsname);
      xsec_subscale[i].back()->SetTitle(soft_subscale[i][j]->GetTitle());
    }
  }


  // and specially for the sums over the subprocesses

  TH1D* xsecsum_scale[nScales];

  TH1D* xsecsum = (TH1D*)xsec_sub[0]->Clone(); 
  xsecsum->SetName("xsecsum"); xsecsum->SetTitle(xsec->GetTitle());

  for ( int i=0 ; i<nScales ; i++ ) { 
    char hname[64];
    sprintf( hname, "xsecsum_scale_%d", i); 
    //    xsec_sub.push_back(  g.convolute_subproc( i, GetPdf, alphaspdf_ , nLoops) );
    xsecsum_scale[i] = (TH1D*)xsec_subscale[i][0]->Clone();
    xsecsum_scale[i]->SetName(hname);
    xsecsum_scale[i]->SetTitle(xsec_scale[i]->GetTitle());
  }
  
  for ( int j=1 ; j<Nsub ; j++ ) { 
    increment( xsecsum, xsec_sub[j] );

    for ( int i=0 ; i<nScales ; i++ ) increment( xsecsum_scale[i],  xsec_subscale[i][j] ); 
  }

  
  xsecdir.pop();



  // now take all the ratios etc

  Directory ratiodir("ratio");
  ratiodir.push();

  TH1D* ratio = divide( xsec, reference ); if ( ratio ) ratio->SetName("ratio");
  
  TH1D* ratio_scale[nScales];

  std::vector<TH1D*> ratio_sub;
  std::vector<TH1D*> ratio_subscale[nScales];
  
  for ( int i=0 ; i<Nsub ; i++ ) {
    char hname[64];
    sprintf( hname, "ratio_sub_%d", i); 
    ratio_sub.push_back(  divide( xsec_sub[i], soft_sub[i]) );
    if ( ratio_sub.back() ) ratio_sub.back()->SetName(hname);
  }

  for ( int i=nScales ; i-- ; ) { 
    
    // get histo for this scale 
    
    char hname[64];
    sprintf( hname, "ratio_scale_%d", i); 

    ratio_scale[i] = divide( xsec_scale[i], soft_scale[i] );
    if ( ratio_scale[i] ) ratio_scale[i]->SetName(hname);

    // get all the subprocess histos for this scale

    for ( int j=0 ; j<Nsub ; j++ ) {
      char hsname[64];
      sprintf( hsname, "ratio_subscale_%d_%d", j, i); 
      ratio_subscale[i].push_back( divide( xsec_subscale[i][j], soft_subscale[i][j] ) );
      if ( ratio_subscale[i].back() ) ratio_subscale[i].back()->SetName(hsname);
    }
  }



  TH1D* subratio = divide( xsec, xsecsum ); if ( subratio ) subratio->SetName("subratio");

  TH1D* subratio_scale[nScales];
  
  for ( int i=nScales ; i-- ; ) { 
    
    // get histo for this scale 
    
    char hname[64];
    sprintf( hname, "subratio_scale_%d", i); 
    
    subratio_scale[i] = divide( xsec_scale[i], xsecsum_scale[i] );
    if ( subratio_scale[i] ) subratio_scale[i]->SetName(hname);
  }

  
  ratiodir.pop();

  fout->Write(); 

  return 0;
}
