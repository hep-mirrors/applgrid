
//   main.cxx        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: main.cxx, v   Sat Mar 15 15:51:49 GMT 2008 sutt


#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include <appl_grid.h>
using appl::grid;

#include <appl_pdf.h>
using appl::appl_pdf;

// #include <genpdf.h>
#include "parden.h"

#include <gtttimer.h>

#include <TH1D.h>
#include <TFile.h>


extern "C" void   dglapeval_(const double& _x, const double& _Q, double* f);
extern "C" void   initmypdf_(const char* name, const int& set);
extern "C" double alphaspdf_(const double& Q);

void lhapdf(const double& _x, const double& _Q, double* _f) { 
  double x = _x; 
  double Q = _Q;
  
  double f[13];

  //  evolvepdf_(x, Q, f); 
  dglapeval_(x, Q, f); 

  for ( int i=0 ;i<13 ; i++ ) _f[i] = f[i]/x;

  return;
}


// jetrad remapping mapping
void jetradpdf(const double& _x, const double& _Q, double* _f) { 

  double f[13];

  lhapdf(_x, _Q, f); 
  
  _f[0] = f[8]; // u
  _f[1] = f[7]; // d
  _f[2] = f[9]; // s
  _f[3] = f[10]; // c
  _f[4] = f[11]; // b
  _f[5] = f[12]; // t

  _f[6]  = f[4];  // ubar 
  _f[7]  = f[5];  // dbar
  _f[8]  = f[3];  // sbar
  _f[9]  = f[2];  // cbar
  _f[10] = f[1];  // bbar
  _f[11] = f[0];  // tbar

  _f[12] = f[6]; // gluon

} 



int main() { 

  cout << "main()" << endl;

  // set up lhapdf :(

  //  const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  // const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6m.LHpdf";  
  //  initpdfset_(_pdfname);


  int Npdf = 0;
  //  numberpdf_(Npdf);

  initmypdf_(_pdfname, Npdf);

  //  double R = 0;
  //  getrenfac_(R);
  
  int iset = 1;
  initpdf_(iset);



#if 0
  //  appl_pdf* mcfm = new appl_pdf("mcfm", 6, mcfm_pdf__);
  appl_pdf* mcfm = appl_pdf::getpdf("mcfm_w");

  cout << "pdf: " << mcfm->name() << endl;

  double f0[13];
  double f1[13];

  double x0 = 0.1211427;
  double x1 = 0.0005462219;
  double Q0 = sqrt(6400.01);
 
  lhapdf( x0, Q0, f0 );
  lhapdf( x1, Q0, f1 );


  double H[6];
  mcfm->evaluate(f0, f1, H);
  for ( int i=0 ; i<mcfm->Nproc() ; i++ ) cout << "\tH[" << i << "]=" << H[i] << endl;

  mcfm->evaluate(f1, f0, H);
  for ( int i=0 ; i<mcfm->Nproc() ; i++ ) cout << "\tH[" << i << "]=" << H[i] << endl;

  delete mcfm;

  return 0;

#endif

  // read in the grid

  // appl::addtopdfmap("mcfm_pdf", mcfm_pdf__);
  // appl::addtopdfmap("mcfm_w", mcfm_pdf__);
 
  // mcfm W gridw
  grid g0("weight_eta4.root");
  // grid g0("weight.root");
  //  grid g0("weight_pt4c.root");
  //  grid g0("weight_pt4f.root");
  
  // jetrad grid
  // grid g0("pdf-grid.root");

  unsigned trim_size = g0.size(); 

  cout << "trim size=" << g0.size() << endl;

  g0.untrim();

  unsigned untrim_size = g0.size(); 

  cout << "trim size=" << g0.size() << endl;
  
  cout << "ratio=" << trim_size*1.0/untrim_size << endl;
  

  //  g0.untrim();
  //  g1.untrim();
  //  g2.untrim();

  
  cout << "doing reference" << endl;
   
  TFile ftmp("weight.root");
  //  TH1D* ref_tmp = (TH1D*)ftmp.Get("grid/obs_bins");
  //  TH1D* reference = new TH1D(*ref_tmp);
  TH1D* reference = (TH1D*)ftmp.Get("grid/obs_bins");
  reference->SetName("reference");
  reference->SetTitle("reference histogram");

  cout << "done reference" << endl;


  TFile outfile("cross.root", "recreate");

  reference->Write();

  cout << "convoluting" << endl;

  // mcfm
  int lowest_order = 0;
  int nloops       = 1;

  cout << "untrim" << endl;
  g0.untrim();
  TH1D* htot   = g0.convolute(lhapdf, alphaspdf_, lowest_order, nloops);   
  htot->SetName("eta");

  cout << "trim" << endl;
  g0.trim();
  TH1D* htot_untrim = g0.convolute(lhapdf, alphaspdf_, lowest_order, nloops);   
  htot_untrim->SetName("eta_untrim");
  

  // jetrad
  //  int lowest_order = 0;
  //  int nloops       = 0;
  //  TH1D* htot = g0.convolute(jetradpdf, alphaspdf_, lowest_order, nloops);   
  //  htot->SetName("pt");


  string num[6] = { "0", "1", "2", "3", "4", "5" }; 

  // leading order and nlo parts by subprocess
  for ( int i=0 ; i<6 ; i++ ) { 
    string name("eta-");
    name += num[i];

    string name_lo("lo-eta-");
    name_lo += num[i];

    string name_nlo("nlo-eta-");
    name_nlo += num[i];

    TH1D* h = g0.convolute_subproc(i, lhapdf, alphaspdf_, lowest_order, nloops);   
    h->SetName(name.c_str());

    TH1D* hlo = g0.convolute_subproc(i, lhapdf, alphaspdf_, lowest_order, 0);   
    hlo->SetName(name_lo.c_str());

    //    TH1D* hnlo = g0.convolute_subproc(i, lhapdf, alphaspdf_, lowest_order, -1);   
    //    hnlo->SetName(name_nlo.c_str());
  }

  outfile.Write();
  outfile.Close();

  return 0;
}








