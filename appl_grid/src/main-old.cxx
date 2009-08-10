//   main.cxx        
//
//   test grid program                    
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: main.cxx, v   Tue Feb 20 20:18:03 GMT 2007 sutt


#include "appl_grid/appl_grid.h"
using appl::grid;

#include "appl_grid/appl_pdf.h"

#include <vector>
using std::vector;

#include <cmath>
using std::pow;
using std::fabs;
using std::log;

// #include <lcg.h>


#include "appl_grid/appl_timer.h"

using namespace appl;

#include "TFile.h"
#include "TH1D.h"


double _fy(double x) { return log(1/x-1); } 
double _fx(double y) { return 1/(1+exp(y)); } 

double _fy1(double x) { return 1-x; } 
double _fx1(double y) { return 1-y; } 

// double _fun(double x) { return pow(x,0.15)*pow((1-x),-3.3); } 
// double _fun(double x) { return 1; }

// extern "C" void pdflo_(const double& x, const double& Q, double* f);

extern "C" void   dglapeval_(const double& _x, const double& _Q, double* f);

void pdf(const double& x, const double& Q, double* f) { 
  //  for ( int i=0 ; i<13 ; i++ ) {  f[i] = pow(x,-0.5)*pow((1-x),3); }
  for ( int i=0 ; i<13 ; i++ )  f[i] = pow(1-x, 4.2)*pow(x, -1.4)*sqrt(Q);
}


double _alphas(const double& Q) {  return 1/log(Q); }

// double _alphas(const double& Q) {  return 1; }



void  _genpdf(const double* f1, const double* f2, double* F) {  
 
  double Q1=0, Q1bar=0, Q2=0, Q2bar=0, D=0, Dbar=0; 

  for ( int i=0 ; i<6 ; i++ )  {  Q1    += f1[i];   Q2    += f2[i];  }
  for ( int i=6 ; i<12 ; i++ ) {  Q1bar += f1[i];   Q2bar += f2[i];  }

  for ( int i=0 ; i<12 ; i++ ) D += f1[i]*f2[i];

  for ( int i=0 ; i<6 ; i++ ) { 
    Dbar += f1[i]*f2[i+6];  
    Dbar += f1[i+6]*f2[i]; 
  }
  
  // spin terms, qq->1/6*1/6  qg->1/6*1/16 gg->1/16*1/16
  // static double s[3] = { 1./36, 1./96, 1./256 };
  static double s[3] = { 1, 1, 1 };
  
  F[0] = s[0] *   D;
  F[1] = s[0] * ( Q1*Q2+Q1bar*Q2bar-D );
  F[2] = s[0] * ( Dbar );
  F[3] = s[0] * ( Q1*Q2bar+Q1bar*Q2-Dbar );
  F[4] = s[1] * ( (Q1+Q1bar)*f2[12] + f1[12]*(Q2+Q2bar) );
  F[5] = s[2] *  f1[12]*f2[12];

}







void sig(double x1, double x2, double Q, double s[13][13], double e) { 
 
  double y1 = _fy(x1);
  double y2 = _fy(x2);

  if ( y1<y2 ) { 
    double t = y1;
    y1 = y2;
    y2 = t;
  } 

  double y = 20*(y1-2.85)*(y1-2.85) + 0.75*(y2-2.85)*(y2-2.85) - 9;

  y = 0;
  if ( x2>0.001 && x2<0.8 && x1>0.001 && x1<0.01 ) y=1/log(Q);
  if ( x1>0.001 && x1<0.8 && x2>0.001 && x2<0.01 ) y=1/log(Q);

  for ( int i1=0 ; i1<13 ; i1++ ) { 
    for ( int i2=0 ; i2<13 ; i2++ ) s[i1][i2] = y;
  }

}





void binwidth(TH1D* h) {
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) { 
    double delta = h->GetBinLowEdge(i+1)-h->GetBinLowEdge(i);
    h->SetBinContent(i, h->GetBinContent(i)/delta);    
    h->SetBinError(i, h->GetBinError(i)/delta);    
    cout << "bin[" << i << "]=" << h->GetBinContent(i) << endl;
  }
} 


void zeroTH1D(TH1D* h) {
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) { 
    h->SetBinContent(i, 0);
    h->SetBinError(i, 0); 
  }
} 




// small, fast 32 bit linear congruential random 
// number generator 
// NB: should only be used for testing, 
//     not for serious work, 

#define _LCG_SEED 1
#define _LCG_M 4294967296.0

static int seed = _LCG_SEED;

double lcg() {
  static unsigned  x=seed;
  x=69069*x+1;
  return ((double)x)/_LCG_M;
}





double fys(double x) { return x*x; }  
double fxs(double y) { return sqrt(y); }  

int main() {


  cout << "running" << endl;

  const int Nobs = 25;

  vector<double> etbins;
  double         aetbins[Nobs+1];
  
  for ( int i=0 ; i<Nobs+1 ; i++ ) { 

    double d=i*sqrt(7000-10)/Nobs;
    etbins.push_back(d*d+10);

    // etbins.push_back( exp(i*(log(7000.)-log(10.))/Nobs + log(10.)) );
    // etbins.push_back( i*(7000.-10.)/Nobs + 10. );
    aetbins[i] = etbins[i];
  }
  

  cout << "creating grid" << endl;



  grid* gp;

  //  igrid::symmetrise();

  //  igrid::reweight(false);

  bool oldgrid=true;

  // can add new user transform pairs 
  igrid::add_transform("s1", _fx, _fy);
  igrid::add_transform("s2", _fx1, _fy1);


  double Q2min = 2500;
  double Q2max = 12250000;

  double xmin = 1e-5;
  double xmax = 0.99;


  double DelQ2 = Q2max-Q2min;
  double Delx  = xmax-xmin; 

  const int Np = 6;

  
  FILE* fptr=fopen("duff.root", "r+");
  if ( fptr ) {
    gp = new grid("duff.root");
    //    gp->optimise(); 

    zeroTH1D(gp->getReference());

    cout << *gp << endl;
  }
  else {
    
    oldgrid=false;

    gp = new grid(etbins, 
		  11, Q2min, Q2max, 4,
		  21,  xmin,  xmax, 4, 
		  "nlojet",
		  0, 0,
		  "f2");

    gp->run()=0;

    // symmetrise
    gp->symmetrise(false);
  }

  grid& g = (*gp);





  TFile f("ftest.root", "recreate");

  TH1D* h  = new TH1D( "h",  "h", g.Nobs(), aetbins );

  // now generate N randon events over x1, x2 and Q2
  // and fill the grid, and make a histogram.

  // Then, loop over the grid, doing the convolution
  // and write another histogram

  // keep it simple here, to allow debugging, then 
  // move on  

  
  double f1[13], f2[13];

  double sigma[13][13];

  
  cout << "filling grid" << endl;
  
  for ( int it=0 ; it<1 ; it++ ) { 
    
    int N = 150000;
    
    if ( g.isOptimised() ) N=20000;
    
    cout << "it=" << it << "\tN=" << N << endl;
    
    struct timeval ittimer = appl_timer_start();  
    
    double fweight = 1.0/N;

    cout << "fweight=" << fweight << endl;

    double percen = 0;
  
    int counted = 0;
    
    for ( int i=0 ; i<N ; i++ ) { 
      
      double x1 = lcg()*Delx+xmin;
      double x2 = lcg()*Delx+xmin;
      double Q2 = lcg()*DelQ2+Q2min;
      
      double Q = sqrt(Q2);
      double e = 6990*lcg()+10; 
      
      if ( e<100 || Q2<2500 ) continue;
      
      pdf(x1, Q, f1);
      pdf(x2, Q, f2);
      
      int nip=Np;
      double gen[Np];
      //    _genpdf(f1,f2,gen);
      //      appl::_genpdf_jetrad(f1,f2,gen);
      //      appl::_genpdf_nlojet(f1,f2,gen);

      //      for ( int ip=0 ; ip<Np ; ip++ ) cout << "gen[" << ip << "]=" << gen[ip] << endl;
      
      sig(x1, x2, Q, sigma, e);
      
      if ( sigma[0][0]==0 ) continue;
      
      double _alpha = _alphas(Q)/(2*M_PI);
      
      counted++;
      
      //      if ( it>0 && counted>1 ) break;
      
      double w = 0; 
      for ( int ip=0 ; ip<Np ; ip++ ) w += _alpha*(gen[ip]*sigma[ip][ip])*fweight;

      if ( it>0 ) 
      { 
	// fill external reference histo
	h->Fill(e, w);
	// fill the internal gtrid reference histogram
	g.getReference()->Fill(e, w);
      }

      double F[Np];
      for ( int ip=0 ; ip<Np ; ip++ ) F[ip]=fweight*sigma[ip][ip];

      if ( g.isOptimised() ) g.fill( x1, x2, Q*Q, e, F, 0);
      else                   g.fill_phasespace( x1, x2, Q*Q, e, F, 0);

    }
    
    double ittime = appl_timer_stop(ittimer);
    cout << "iteration " << it << "\ttime=" << ittime << " ms" << endl;

    cout << "it=" << it << "\toptimised=" << g.isOptimised() << endl; 
   
    if ( g.isOptimised() ) g.run()++;
    if ( it<1 && !g.isOptimised() )  g.optimise();
  }
  
  binwidth(h);
  binwidth(g.getReference());


  cout << "filled grid" << endl;
  
  // Now we can try to build the cross section 
  // from the grid we've just created
  
  cout << "trim" << endl;
  g.trim();

  //  cout << "untrim" << endl;
  //  g.untrim();

  TH1D* hgg = NULL;

  cout << "convoluting " << endl;

  struct timeval conv_timer = appl_timer_start(); 
 
  //  hgg = g.convolute(pdf, appl::_genpdf_jetrad, _alphas, 1, 0);
  //  hgg = g.convolute(pdf, appl::_genpdf_jetrad, _alphas, 1, 0);
  hgg = g.convolute(pdf, _alphas, 1, 0);

  //  hgg = g.convolute(pdf, _genpdf, _alphas, 1, 0);
  
  double conv_time = appl_timer_stop(conv_timer);

  cout << "convolution time " << conv_time << " ms" << endl;
    
  hgg->SetName("hgg");
  

  // scale convoluted cross section by the number 
  // of runs - could be built into the grid itself
  cout << "runs=" << g.run() << endl;
  hgg->Scale(g.run());

  //  could have done it this way
  //  g *= 1.0/g.run();

  
  cout << "trim" << endl;

  TH1D* hcopy = new TH1D(*hgg);
  hcopy->SetName("diff");

  double mean = 0;  
  double dN   = 0;
  
  for ( int i=1 ; i<=hcopy->GetNbinsX() ; i++ ) { 
    double dh   = h->GetBinContent(i);
    double dhgg = hgg->GetBinContent(i);
  
    if ( dh ) { 
      double r = dhgg/dh;
      mean += r;
      dN++;
      if ( dh ) hcopy->SetBinContent(i, r);      
    }
  }
  
  cout << "mean=" << mean/dN << endl;

  hcopy->SetMinimum(0.999);
  hcopy->SetMaximum(1.001);
  
  f.Write();
  f.Close();
  
 
  //  if ( !oldgrid ) g.optimise(10,21,21); 

  cout << "writing" << endl;
  g.Write("duff.root");

  return 0;
} 
