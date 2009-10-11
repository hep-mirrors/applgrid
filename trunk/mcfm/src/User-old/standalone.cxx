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

#include <appl_grid.h>
using appl::grid;

#include <appl_pdf.h>
using appl::appl_pdf;

// #include <sppl_timer.h>

#include <TH1D.h>
#include <TFile.h>

#include "parden.h"

extern "C" void   dglapeval_(const double& _x, const double& _Q, double* f);
extern "C" void   initmypdf_(const char* name, const int& set);
extern "C" double alphaspdf_(const double& Q);




// wrapper to get the basic pdf rather than x*pdf
// gavins evolution code
 
void lhapdf(const double& x, const double& Q, double* f) { 

  double xf[13];

  //  evolvepdf_(x, Q, xf); 
  //  cout << "\tevo=" << xf[6];

  dglapeval_( x, Q, xf); 
  //  cout << "\tdgl=" << xf[6] << endl;

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

  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i]/x;

  return;
}



int main() { 

  cout << "main()" << endl;

  // set up lhapdf :(

  //  const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  const string _pdfname = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  // const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6m.LHpdf";  
  //  initpdfset_(_pdfname);


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

  string _filename[3] = { "weight_eta4.root", 
			  "weight_pt4c.root", 
			  "weight_pt4f.root" };

  string _label[3]    = { "eta4", "pt4c", "pt4f" };


  TFile outfile("cross.root", "recreate");

  // do all three grids, eta, pt central and forward
  for ( int ig=0 ; ig<3 ; ig++ ) { 

    string& label    = _label[ig];

    string  filename = "weight_" + label + ".root";
    
    cout << label << " -----------------" << endl;

    // mcfm W grids
    grid g0(filename);

    cout << g0 << endl;

    outfile.cd();
    Directory d(label);
    d.push();

    
    g0.trim();
    unsigned trim_size = g0.size(); 
    
    g0.untrim();
    unsigned untrim_size = g0.size(); 

    cout << "trimmed size="     << trim_size 
	 << "\tuntrimmed size=" << untrim_size
	 << "\tratio=" << trim_size*1.0/untrim_size << endl;
  

    //    cout << "doing reference" << endl;
   
    TFile ftmp(filename.c_str());

    outfile.cd();
    d.push();

    //  TH1D* ref_tmp = (TH1D*)ftmp.Get("grid/obs_bins");
    //  TH1D* reference = new TH1D(*ref_tmp);
    //  WHY!!!!! are there root "directories" that get changed into 
    //  and out of etc. why not just have histo->Write(file) to write
    //  to a file, or even better, file->Save(histo)
    TH1D* reference = (TH1D*)ftmp.Get("grid/reference");
    reference->SetName((label+"-reference").c_str());
    reference->SetTitle((label+" reference histogram").c_str());
    
    //    cout << "done reference" << endl;
    
    reference->Write();

    cout << "convoluting" << endl;
    
    // lhapdf
    TH1D* htot   = g0.convolute(lhapdf, alphaspdf_ );   
    htot->SetName((label+"-lha").c_str());

    // gspdf
    TH1D* htot2   = g0.convolute(lhapdf2, alphaspdf_ );   
    htot2->SetName((label+"-sal").c_str());
    
    //    htot->Write();

    // jetrad
    //  int lowest_order = 0;
    //  int nloops       = 0;
    //  TH1D* htot = g0.convolute(jetradpdf, alphaspdf_, lowest_order, nloops);   
    //  htot->SetName("pt");
    
    
    d.push();

    // do all the subprocesses also for the leading order 
    // and nlo parts by subprocess

    string num[6] = { "0", "1", "2", "3", "4", "5" }; 

    for ( int i=0 ; i<6 ; i++ ) { 

      Directory dsub("sub-"+num[i]);
      dsub.push();

      string name(label+"-");
      name += num[i];
      
      string name_lo("lo-"+label+"-");
      // string name_lo(label+"-");
      name_lo += num[i];
      
      string name_nlo("nlo-"+label+"-");
      // string name_nlo(label+"-");
      name_nlo += num[i];

      // subprocess
      
      TH1D* h = g0.convolute_subproc(i, lhapdf, alphaspdf_ );   
      h->SetName(name.c_str());
      // h->Write();
      
      // leading order contribution 
      // dlo.push();
      TH1D* hlo = g0.convolute_subproc(i, lhapdf, alphaspdf_, 0 );   
      hlo->SetName(name_lo.c_str());
      //  hlo->Write();
      //  dlo.pop();

      // nlo contribution 
      // dnlo.push();
      TH1D* hnlo = g0.convolute_subproc(i, lhapdf, alphaspdf_, -1 );   
      hnlo->SetName(name_nlo.c_str());
      //  hnlo->Write();
      //  dnlo.pop();

      dsub.pop();
    }
    d.pop();
  }
  
  outfile.Write();
  outfile.Close();
  
  return 0;
}








