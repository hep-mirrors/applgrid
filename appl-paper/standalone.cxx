#include <iostream>
#include <string.h>

using std::cerr; 
using std::cout;
using std::endl;  

#include "appl_grid/appl_grid.h"
using appl::grid;

#include "appl_grid/appl_pdf.h"
using appl::appl_pdf;

#include "appl_grid/appl_timer.h"

#include <TH1D.h>
#include <TFile.h>
//
//
bool debug=true;
//
extern "C" 
{
  
  void dglapeval_(const double& _x, const double& _Q, double* f);
  void dglapevalsplit_(const double& _x, const double& _Q, const int&, const int&, double* f);
  void initmypdf_(const char* name, const int& set);
  
  double alphaspdf_(const double& Q);
} 
// 
//
//
static const double pb_fac = 3.89379656e8 ;    // conversion GeV^2 -> pb  
static const int nScales = 5;
static const double mur[nScales] = {1.0, 0.5, 2.00, 1.0, 0.5};
static const double muf[nScales] = {1.0, 0.5, 2.00, 0.5, 1.0};

static const int nLoops = 1;
static const int nFlavours = 5;

// reference histograms for 7 subprocesses
TH1D* soft_sub  [7];
TH1D* grid_sub  [7];
TH1D* ratio_sub [7];
// reference histogram for 7 subprocesses and  3 scales
TH1D* soft_subscale  [7][nScales];            
TH1D* grid_subscale  [7][nScales];            
TH1D* ratio_subscale [7][nScales];            
// reference histograms to test renormalisation & factorisation scales deps
TH1D* soft_scale  [nScales];
TH1D* grid_scale  [nScales];
TH1D* ratio_scale [nScales];

TH1D* reference;
TH1D* gridObservable;
TH1D* ratio;

TH1D* soft_nlo;
TH1D* grid_nlo;
TH1D* ratio_nlo;

TH1D* soft_lo;
TH1D* grid_lo;
TH1D* ratio_lo;

void readReferenceHistograms(TFile* myfile)
{

  char histname[3000];  

  cout << "readReferenceHistograms() \t\t pwd=" << gDirectory->GetName() << endl;
  string refDirName = "addReference";
  Directory refDir(refDirName);
  refDir.push();
  cout << "\t readReferenceHistograms() \t pwd=" << gDirectory->GetName() << endl;  

//TC #if 0
  for(int is=0; is<7; is++) 
    {
      sprintf(histname,"%s/soft_sub_%i",refDirName.c_str(),is);
      soft_sub[is] = (TH1D*)(myfile->Get(histname));
      for(int ir=0; ir < nScales; ir++)
	{
	  sprintf(histname,"%s/soft_subscale_%d_%d", refDirName.c_str(), is, ir);
	  soft_subscale[is][ir] = (TH1D*)myfile->Get(histname);
          if (soft_subscale[is][ir]) cout<<" read "<<histname<<endl;
	  if(is == 0)
	    {
	      sprintf(histname,"%s/soft_scale_%d",refDirName.c_str(),ir);
	      soft_scale[ir] = (TH1D*)myfile->Get(histname);
              if ( soft_scale[ir]) cout<<" read "<<histname<<endl;
	    }
        
	}
    }
//TC #endif

  sprintf(histname,"%s/soft_%s",refDirName.c_str(), "lo");
  soft_lo = (TH1D*)myfile->Get(histname);
  if (soft_lo) cout<<" read "<<histname<<endl;

  sprintf(histname,"%s/soft_%s",refDirName.c_str(), "nlo");
  soft_nlo = (TH1D*)myfile->Get(histname);
  if (soft_nlo) cout<<" read "<<histname<<endl;
   
  refDir.pop();
  
  reference = (TH1D*)myfile->Get("grid/reference");

}

//
//
void removeErrors(TH1D* hh)
{
  int nBins = hh->GetNbinsX();
  for (int iBin = 0; iBin <= nBins + 1; iBin++)
    {
      hh->SetBinError(iBin,0.);
    }
}



// wrapper to get the basic pdf rather than x*pdf
// gavins evolution code

void gavin(const double& x, const double& Q, double* f) { 
  
  double xf[13];

  dglapeval_( x, Q, xf);    
  //if (debug) cout << "\t evo=" << xf[6];
  //if (debug) cout << " x= "<<" Q= "<<Q<<"\tdgl=" << xf[6] << endl;
 double invx=0.;
 if (x!=0.) invx=1./x;
 for ( int i=0; i<13 ; i++ ) f[i] = xf[i]*invx;
 return; 
}

void gavinSplit(const double& x, const double& Q, double* f) { 
  
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

  void Dump_ratio(TH1D *ratio){

  if (!ratio) {cout<<" Histo not found "<<endl; return;}
  ratio->Print();

  Int_t nbin=ratio->GetNbinsX();
  for (Int_t i=0; i<nbin; i++) {
    cout<<ratio->GetName()<<" x= "<<ratio->GetBinCenter(i)<<" y= "<<ratio->GetBinContent(i)<<endl;
  }
  return;
 }


int main(int argc, char** argv) 
{ 
  
  cout << "main()" << endl;
  //  string gridFilePrefix = "./output/";
  //  string gridFilePrefix = "./parameters/fig6/";
  //  string gridFileName = gridFilePrefix + "weight_"; 
  string gridFilePrefix = "./";
  string gridFileName;
  // string  outFileName = 
  string  outFileName = "out.root";
  
  if ( argc<2 ) 
    {
      cerr << "weight file not specified" << endl; 
      return -1;
    }
  else
    {
      gridFileName += argv[1];
      // outFileName += argv[1];
    }
  
  //  gridFileName += ".root";
  //  outFileName += ".root";
  if (FILE* testfile = fopen(gridFileName.c_str(),"r"))
    {
      fclose(testfile);
    }
  else
    {
      cout <<"No such ROOT-file "<<gridFileName<<endl;
      return -1;
    }
  
  const string _pdfname = "PDFsets/cteq6mE.LHgrid";  
  //  const string _pdfname = "PDFsets/cteq6mE.LHgrid";  
  int Npdf = 0;
  // setup gavins code
  initmypdf_(_pdfname.c_str(), Npdf);
  cout << "lhapdf..." << endl;
  //
  //
  // read in the grid
  //
  //
  TFile outfile(outFileName.c_str(), "recreate");
  TFile infile(gridFileName.c_str(), "read");
  
  grid* g = new grid(gridFileName);
  //g.print();
  
  g->trim();
  
  readReferenceHistograms(&infile);
  
  //
  // observable
  //  
  cout << "convoluting observable" << endl;

  gridObservable = (TH1D*)g->convolute(gavin, alphaspdf_);
  gridObservable->Scale(pb_fac);
  gridObservable->SetName("gridObservable");
  gridObservable->SetTitle("Observable calculated via APPLgrid");
  
  ratio = (TH1D*)reference->Clone("ratio");
  ratio->Divide(gridObservable);
  ratio->SetTitle("reference/grid ratio");

  cout<<" Ratio reference/grid "<<endl;
  Dump_ratio(ratio);

  outfile.cd();
  reference->Write();
  gridObservable->Write();
  removeErrors(ratio);
  ratio->Write();
  //
  //  convoluting per order
  //

//   cout << "convoluting per order ...";

//   grid_lo = g->convolute(gavin, alphaspdf_, 0);
//   grid_lo->Scale(pb_fac);
//   grid_lo->SetName("grid_lo");
//   grid_lo->SetTitle("Observable @LO calculated via APPLgrid");

//   if (!soft_lo) {cout<<" soft_lo not found ! "<<endl; return 0;}
//   ratio_lo = (TH1D*)soft_lo->Clone("ratio_lo");
//   ratio_lo->Divide(grid_lo);
//   ratio_lo->SetTitle("reference/grid ratio @LO");

//   cout <<"LO done, working with NLO ...";

//   grid_nlo = g->convolute(gavin, alphaspdf_, -1);
//   grid_nlo->Scale(pb_fac);
//   grid_nlo->SetName("grid_nlo");
//   grid_nlo->SetTitle("Observable @NLO calculated via APPLgrid");

//   if (!soft_nlo) {cout<<" soft_nlo not found ! "<<endl; return 0;}
//   ratio_nlo = (TH1D*)soft_nlo->Clone("ratio_nlo");

//   ratio_nlo->Divide(grid_nlo);
//   ratio_nlo->SetTitle("reference/grid ratio @NLO");

//   cout<<" done"<<endl;

//   outfile.cd();
   Directory order("order");
//   order.push();


 

//   soft_lo->Write();
//   grid_lo->Write();
//   removeErrors(ratio_lo);
//   ratio_lo->Write();

//   soft_nlo->Write();
//   grid_nlo->Write();
//   removeErrors(ratio_nlo);
//   ratio_nlo->Write();



  order.pop();
  //
  // subprocesses   
  //
  cout << "convoluting observable per subprocess" << endl;
  char histtitle[1000];
  for (int iSub = 0; iSub < 7; iSub++)
    {
      cout<<"\tworking with subprocess "<<iSub<<" ..."<<endl;
      //
      grid_sub[iSub] = g->convolute_subproc(iSub, gavin, alphaspdf_);
      grid_sub[iSub]->Scale(pb_fac);

      sprintf(histtitle,"grid_sub_%d",iSub);
      grid_sub[iSub]->SetName(histtitle);
      sprintf(histtitle,"Observable calculated via APPLgrid for subProcess: %d", iSub);
      grid_sub[iSub]->SetTitle(histtitle);
      //

      sprintf(histtitle,"ratio_sub_%d",iSub);
      if (!soft_sub[iSub]) {cout<<" soft_sub["<<iSub<<"])  not found ! "<<endl; return 0;}
      ratio_sub[iSub] = (TH1D*)soft_sub[iSub]->Clone(histtitle);
      ratio_sub[iSub]->Divide(grid_sub[iSub]);
      sprintf(histtitle,"reference/grid ratio for subProcess: %d", iSub);
      ratio_sub[iSub]->SetTitle(histtitle);

    }



  outfile.cd();
  Directory dirSubProcess("subProcess");
  dirSubProcess.push();
  for (int iSub = 0; iSub < 7; iSub++)
    {
      grid_sub[iSub]->Write();
      //soft_sub[iSub]->Write();

      removeErrors(ratio_sub[iSub]);
      ratio_sub[iSub]->Write();
    }
  dirSubProcess.pop();


  //
  // renormalisation and factorisation scales
  //
  cout << "convoluting observable per scale and subprocess" << endl;

  for (int iScale = 0; iScale < nScales; iScale++)
    {
      cout<<"\t  working with scale "<<iScale<<" ..."<<endl;
      grid_scale[iScale] = g->convolute(gavin, alphaspdf_, nLoops, 
					mur[iScale], muf[iScale], gavinSplit);
      grid_scale[iScale]->Scale(pb_fac);

      sprintf(histtitle,"grid_scale_%d",iScale);
      grid_scale[iScale]->SetName(histtitle);
      sprintf(histtitle, "Observable calculated via APPLgrid for scales: #mu_{R} = %3.2f,  #mu_{F} = %3.2f", mur[iScale], muf[iScale]);
      grid_scale[iScale]->SetTitle(histtitle);
      //      grid_scale->Print();
      //
      sprintf(histtitle,"ratio_scale_%d",iScale);
      if (!soft_scale[iScale]) {cout<<" soft_scale["<<iScale<<"])  not found ! "<<endl; return 0;}
      ratio_scale[iScale] = (TH1D*)soft_scale[iScale]->Clone(histtitle);

      ratio_scale[iScale]->Divide(grid_scale[iScale]);

      sprintf(histtitle,"reference/grid ratio for  scales: #mu_{R} = %3.2f,  #mu_{F} = %3.2f", mur[iScale], muf[iScale]);
      ratio_scale[iScale]->SetTitle(histtitle);

      cout<<" Ratio reference/grid "<<" mur= "<<mur[iScale]<<" muf= "<<muf[iScale]
          <<endl;
      Dump_ratio(ratio_scale[iScale]); 

      for (int iSub = 0; iSub < 7; iSub++)
	{
	  cout<<"\t\tworking with subprocess "<<iSub<<" ..."<<endl;
	  //
	  grid_subscale[iSub][iScale] = 
//
//TC>>>>>>>>
// 	    g->convolute_subproc(iSub, gavin, alphaspdf_,
// 				 nLoops, 
// 				 std::sqrt(mur[iScale]), 
// 				 std::sqrt(muf[iScale]), gavinSplit);
 	    g->convolute_subproc(iSub, gavin, alphaspdf_,
 				 nLoops, 
 				 mur[iScale], 
 				 muf[iScale], gavinSplit);
	  grid_subscale[iSub][iScale]->Scale(pb_fac);
//TC<<<<<<<<
	  sprintf(histtitle, "grid_subscale_%d_%d", iSub, iScale);
	  grid_subscale[iSub][iScale]->SetName(histtitle);
	  sprintf(histtitle,"Observable calculated via APPLgrid for scales: #mu_{R} = %3.2f,  #mu_{F} = %3.2f in subProcess = %d", mur[iScale], muf[iScale], iSub);
	  grid_subscale[iSub][iScale]->SetTitle(histtitle);
	  
	  //
	  sprintf(histtitle,"ratio_subscale_%d_%d", iSub, iScale);
          if (!soft_subscale[iSub][iScale]) {cout<<" soft_subscale["<<iSub<<"]["<<iScale<<"])  not found ! "<<endl; return 0;}
          ratio_subscale[iSub][iScale] = (TH1D*)soft_subscale[iSub][iScale]->Clone(histtitle);
	  ratio_subscale[iSub][iScale]->Divide(grid_subscale[iSub][iScale]);
	  sprintf(histtitle,"reference/grid ratio for  scales: #mu_{R} = %3.2f,  #mu_{F} = %3.2f in subProcess = %d", mur[iScale], muf[iScale], iSub);
	  ratio_subscale[iSub][iScale]->SetTitle(histtitle);
	}
    }


  outfile.cd();
  Directory dirScale("scale");
  dirScale.push();
  for (int iScale = 0; iScale < nScales; iScale++)
    {
      grid_scale [iScale]->Write();
      soft_scale [iScale]->Write();
      removeErrors(ratio_scale[iScale]);
      ratio_scale[iScale]->Write();
    }
  gDirectory->ls();
  dirScale.pop();

  outfile.cd();
  Directory dirSubScale("subScale");
  dirSubScale.push(); 
  for (int iScale = 0; iScale < nScales; iScale++)
    {
      for (int iSub = 0; iSub < 7; iSub++)
	{
	  grid_subscale [iSub][iScale]->Write();
          soft_subscale [iSub][iScale]->Write();
	  removeErrors(ratio_subscale[iSub][iScale]);
	  ratio_subscale[iSub][iScale]->Write();
	}
    }
  dirSubScale.pop();
  //
  //   writing results
  //
  //  outfile.Write();
  outfile.Close();
  infile.Close();
  
  delete g;
  return 0;
}








