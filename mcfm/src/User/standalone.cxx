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




TH1D* diffs(TH1D* h0, TH1D* h1) { 
  TH1D* h = new TH1D(*h0);
  
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) { 
    double hc = h0->GetBinContent(i);
    double hr = h1->GetBinContent(i);
    
    double diff = 0;
    if ( hr!=0 ) diff = (hc/hr-1);
    
    h->SetBinContent(i, diff);
  }
  
  return h;
}



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
// lhapdf evolution 
 
void lhapdf(const double& x, const double& Q, double* f) { 
  double _x = x; 
  double _Q = Q;

  double xf[13];

  //  evolvepdf_(x, Q, xf); 
  //  cout << "\tevo=" << xf[6];

  evolvepdf_(_x, _Q, xf); 
  //  cout << "\tdgl=" << xf[6] << endl;

  //  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i]/x;
  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i];

  return;
}



// wrapper to get the basic pdf rather than x*pdf
// gavins evolution code

void gavin(const double& x, const double& Q, double* f) { 
  
  double xf[13];

  dglapeval_( x, Q, xf); 
  //  cout << "\tevo=" << xf[6];

  //  dglapeval_(_x, _Q, xf); 
  //  cout << "\tdgl=" << xf[6] << endl;

  //  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i]/x;
  for ( int i=0 ;i<13 ; i++ ) f[i] = xf[i];

  return;
}



int main(int argc, char** argv) { 

  cout << "main()" << endl;

  if ( argc<2 ) { 
    cerr << "weight file not specified" << endl; 
    return -1;
  }

  // set up lhapdf :(

  //  const char _pdfname[256] = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  const string _pdfname = "/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid";  
  // const string _pdfname = "cteq6mE.LHgrid";  
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


  string _label[3]    = { "eta4", "pt4c", "pt4f" };

  string outname = "cross";


  TFile outfile("cross.root", "recreate");

  // do all three grids, eta, pt central and forward
  for ( int ig=0 ; ig<3 ; ig++ ) { 

    string& label    = _label[ig];

    string  filename = "weight_" + label + ".root";
    
    filename = argv[1] + filename;

    cout << label << " -----------------" << endl;

    // 261

    // mcfm W grids
    grid g0(filename);

    // 263

    struct timeval trimtimer = appl::appl_timer_start();
    g0.untrim();
    double trimtime = appl::appl_timer_stop(trimtimer);

    struct timeval ttrimtimer = appl::appl_timer_start();
    g0.trim();
    double ttrimtime = appl::appl_timer_stop(ttrimtimer);

    cout << "trimtime=" << trimtime << "\tttrimtime=" << ttrimtime << endl;


    // 269

    //    for ( ; ; ) ; 

    cout << g0 << endl;
    g0.trim();


    //    g0.print();


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

    g0.trim();
    
    // lhapdf
    cout << "dglapeval" << endl;
    TH1D* htot   = g0.convolute(gavin, alphaspdf_ );   
    htot->SetName((label+"-lha").c_str());

    // gspdf
    cout << "lhapdf" << endl;
    TH1D* htot2   = g0.convolute(lhapdf, alphaspdf_ );   
    htot2->SetName((label+"-sal").c_str());
    

    TH1D* dff = diffs(htot, reference);
    dff->SetName((label+"-difference").c_str());
    dff->SetMaximum( 0.001);
    dff->SetMinimum(-0.001);

    //    htot->Write();

    //    g0.untrim();

    //    TH1D* htot3   = g0.convolute(lhapdf, alphaspdf_ );   
    //    delete htot3;

    //    g0.trim();


    if ( ig<2 ) { 
    cout << "size = " << untrim_size << "\t " << trim_size << "\t" 
	 << ( untrim_size>0 ? 1.0*trim_size/untrim_size : 0 ) 
	 << "\tdeviation rms = " << rmsdev(htot, reference) << "\t max = " << maxdev(htot,reference) << endl; 
    }    
    else { 
    cout << "size = " << untrim_size << "\t " << trim_size << "\t" 
	 << ( untrim_size>0 ? 1.0*trim_size/untrim_size : 0 ) 
	 << "\tdeviation rms = " << rmsdev(htot, reference, 1) << "\t max = " << maxdev(htot,reference, 1) << endl; 
    }

    // jetrad
    //  int lowest_order = 0;
    //  int nloops       = 0;
    //  TH1D* htot = g0.convolute(jetradpdf, alphaspdf_, lowest_order, nloops);   
    //  htot->SetName("pt");
    
    
    d.push();

    // do all the subprocesses also for the leading order 
    // and nlo parts by subprocess

#if 0

    for ( int i=0 ; i<g0.subProcesses() ; i++ ) { 

      char _sub[8];
      sprintf(_sub, "%d", i);
      string num(_sub);

      cout << "sub process " << num << endl;

      Directory dsub("sub-"+num);
      dsub.push();

      string name(label+"-");
      name += num;
      
      string name_lo("lo-"+label+"-");
      // string name_lo(label+"-");
      name_lo += num;
      
      string name_nlo("nlo-"+label+"-");
      // string name_nlo(label+"-");
      name_nlo += num;

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

#endif

    d.pop();
  }
  
  outfile.Write();
  outfile.Close();
  
  return 0;
}








