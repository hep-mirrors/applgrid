//
//   standalone.cxx        
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: main.cxx, v   Sat Mar 15 15:51:49 GMT 2008 sutt


#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "appl_grid/appl_grid.h"
using appl::grid;

#include "appl_grid/appl_pdf.h"
using appl::appl_pdf;

#include "appl_grid/appl_timer.h"
// #include <gtttimer.h>

#include <TH1D.h>
#include <TFile.h>

#include "parden.h"

extern "C" void   dglapeval_(const double& _x, const double& _Q, double* f);
extern "C" void   initmypdf_(const char* name, const int& set);
extern "C" double alphaspdf_(const double& Q);




double rmsdev(TH1D* h, TH1D* href, int cheat=0) { 
  
  double sum=0;

  for ( int i=1 ; i<=h->GetNbinsX()-cheat ; i++ ) { 
    double hc = h->GetBinContent(i);
    double hr = href->GetBinContent(i);

    double diff = 0;
    if ( hr!=0 ) diff = (hc/hr-1);

    if ( hr!=0 ) sum += diff*diff;
  }
 
  return sqrt(sum/h->GetNbinsX());
}





double maxdev(TH1D* h, TH1D* href, int cheat=0) { 
  
  double sum=0;

  for ( int i=1 ; i<=h->GetNbinsX()-cheat ; i++ ) { 
    double hc = h->GetBinContent(i);
    double hr = href->GetBinContent(i);

    double diff = 0;
    if ( hr!=0 ) diff = (hc/hr-1);

    if ( fabs(diff)>sum ) sum = fabs(diff);
  }
 
  return sum;
}



// wrapper to get the basic pdf rather than x*pdf
// gavins evolution code
 
void lhapdf(const double& x, const double& Q, double* f) { 

  double xf[13];

  //  evolvepdf_(x, Q, xf); 
  //  cout << "\tevo=" << xf[6];

  dglapeval_( x, Q, xf); 
  //  cout << "\tdgl=" << xf[6] << endl;

  //  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i]/x;
  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i]/x;

  return;
}



// wrapper to get the basic pdf rather than x*pdf
// lhapdf evolution 

void lhapdf2(const double& _x, const double& _Q, double* f) { 
  double x = _x; 
  double Q = _Q;
  
  double xf[13];

  evolvepdf_(x, Q, xf); 
  //  cout << "\tevo=" << xf[6];

  //  dglapeval_(_x, _Q, xf); 
  //  cout << "\tdgl=" << xf[6] << endl;

  //  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i]/x;
  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i]/x;

  return;
}



int main(int argc, char** argv) { 

  cout << "main()" << endl;

  if ( argc<2 ) return -1;

  string filename(argv[1]);

  cout << "argv[1]=" << filename << "<" << endl;

  
  // set up lhapdf :(

  //  const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  //  const string _pdfname = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  //  const string _pdfname = "/home/sutt/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  const string _pdfname = "cteq6mE.LHgrid";  
  //  const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6m.LHpdf";  
  //  initpdfset_(_pdfname);

  // const string _pdfname = "cteq6mE.LHgrid";  


  int Npdf = 0;

  // setup gavins code
  initmypdf_(_pdfname.c_str(), Npdf);

  cout << "lhapdf..." << endl;

  // initialise lhapdf if needed
  // initpdfset_(_pdfname.c_str());

  //  cout << "set up lhapdf" << endl;

  //  initpdfset_(_pdfname.c_str());

  //  double R = 0;
  //  getrenfac_(R);
  
  //  int iset = 1;
  //  initpdf_(iset);


  // read in the grid
  

  grid g0(filename);

  TFile outfile("cross-duff.root", "recreate");

  g0.trim();
  
  TH1D* htot3   = g0.convolute(lhapdf2, alphaspdf_ );   

  string labels[16] = { "0", "1", "2", "3", 
			"4", "5", "6", "7", 
			"8", "9", "a", "b", 
			"c", "d", "e", "f"  };

  for ( int i=0 ; i<appl_pdf::getpdf(g0.getGenpdf())->Nproc() ; i++ )  {

    TH1D* h = g0.convolute_subproc(i, lhapdf2, alphaspdf_ );  
    h->SetName(labels[i].c_str());

  }

  outfile.Write();
  outfile.Close();
  
  return 0;
}








