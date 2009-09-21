// #include "/Users/sutt/rootmacros/brun.h"
#include <stdio.h>


void plot() { 

  //  gStyle = setstyle(0);  
  //  gStyle->cd();

  TCanvas* tc = new TCanvas("plots", "plots", 1000, 300);  
  tc->Draw();
  tc->Divide(3);
  tc->cd();
  

  TFile f("xsec.root");

  int offset=0;
  for ( int i=0 ; i<3 ; i++ ) {
    char hname[64];
    //    if ( i>0 ) offset = 2;
    sprintf( hname, "ratio/ratio_scale_%d", i+offset);
    TH1D* h=(TH1D*)f.Get(hname);
    tc->cd(i+1);
    h->DrawCopy();
    tc->Update();
  }

}
