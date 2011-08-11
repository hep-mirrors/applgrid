//   main.cxx        
//
//   test grid program                    
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: tst1.cxx, v0   Tue Feb 20 20:18:03 GMT 2007 sutt


#include "appl_grid/appl_grid.h"
using appl::grid;

#include "appl_grid/appl_pdf.h"

#include <vector>
using std::vector;

#include <cmath>
// using std::pow;
// using std::fabs;
// using std::log;

#include "lcg.h"


#include "gtttimer.h"

#include "TFile.h"
#include "TH1D.h"


void pdf(const double& x, const double& Q, double* f) { 
  //  for ( int i=0 ; i<13 ; i++ ) {  f[i] = pow(x,-0.5)*pow((1-x),3); }
  for ( int i=0 ; i<13 ; i++ )  f[i] = pow(1-x, 4.2)*pow(x, -1.4)*sqrt(Q);
}


double _alphas(const double& Q) {  return 1/log(Q); }
// double _alphas(const double& Q) {  return 1; }



double _fy(double x) { return log(1/x-1); } 
double _fx(double y) { return 1/(1+exp(y)); } 


void sig(double x1, double x2, double Q, double s[13][13], double e) { 
 
  double y1 = _fy(x1);
  double y2 = _fy(x2);

  //  if ( y1<y2 ) { 
  //    double t = y1;
  //    y1 = y2;
  //    y2 = t;
  //  } 

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
    cout << "bin[" << i << "]=" << h->GetBinContent(i) << endl;
  }
} 



int main() {
  cout << "running" << endl;

  grid g("duff.root");
  g.trim();

  TFile f("ctest.root", "recreate");

  // Now we can try to build the cross section 
  // from the grid we've just read

  cout << "convoluting " << endl;

  struct timeval conv_timer = gtttimer_start();  
  TH1D* hgg = g.convolute(pdf, _alphas, 1, 0);
  double conv_time = gtttimer_stop(conv_timer);

  hgg->SetName("hgg");
  
  cout << "convolution time " << conv_time << " ms" << endl;
    
  cout << "runs=" << g.run() << endl;
  hgg->Scale(g.run());
  
  f.Write();
  f.Close();
  
  cout << "writing" << endl;
 
  return 0;
} 
