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
#include <TGraph.h>
#include <TFile.h>
//
//
//
extern "C" 
{
  
  void evolvepdf_(const double& , const double& , double* );
  void initmypdf_(const char* name, const int& set);
  void initmypdfset_(const char* name);
  void initmypdfmember_(const int& set);
  
  double alphaspdf_(const double& Q);
} 
//
//
//
static const int debug = 0;
static const int nScales = 9;
static const double muR[nScales] = {1.0, 0.5, 1.0, 2.0, 0.5, 2.0, 0.5, 1.0, 2.0};
static const double muF[nScales] = {1.0, 0.5, 0.5, 0.5, 1.0, 1.0, 2.0, 2.0, 2.0};
//
static const int nLoops = 1;
static const int nFlavours = 5;
//
static const int nPDFs = 5;
static const string PDFName[nPDFs] = {"cteq6mE.LHgrid", "cteq66.LHgrid", "HERAPDF01.LHpdf", "NNPDF10_100.LHgrid", "MSTW2008nlo90cl.LHgrid"};
static const int nPDFMembers[nPDFs]={41, 45, 23, 101, 41};
//
TH1D* gridObservable, *reference, *ratioGridObservable;
TH1D* sigma[nPDFs][120] = {0}, *ratio[nPDFs][120] = {0};
TH1D* sigmaScale[nPDFs][nScales] = {0}, *ratioScale[nPDFs][nScales] = {0};
//
void removeErrors(TH1D* hh)
{
  int nBins = hh->GetNbinsX();
  for (int iBin = 0; iBin <= nBins + 1; iBin++)
    {
      hh->SetBinError(iBin,0.);
    }
}
#include <sstream>

template <class T>
inline std::string to_String (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

// wrapper to get the basic pdf rather than x*pdf
// gavins evolution code


void GetPdf(const double& x, const double& Q, double* xf) { 

  evolvepdf_( x, Q, xf);    
  return; 
}

//
//--------------------------------------------------------------------------
//
int main(int argc, char** argv) 
{ 
  cout << " pdf errors and scale variation standalone main()..." << endl;
  //
  // preparing root file
  //
  //
  //
  //
  string gridFilePrefix = "./output/";
  string gridFileName = gridFilePrefix + "trijet_"; 
  string  outFileName = gridFilePrefix + "spVar_trijet_"; 
  string observableName = "";
  if ( argc<2 ) 
    {
      cerr << "weight file not specified" << endl; 
      return -1;
    }
  else
    {
      observableName += argv[1];
      gridFileName += observableName + "_3.root";
      outFileName  += observableName + "_3.root";
    }
  //
  //
  //  
  if (FILE* testfile = fopen(gridFileName.c_str(),"r"))
    {
      fclose(testfile);
    }
  else
    {
      cout <<"No such ROOT-file "<<gridFileName<<endl;
      return -1;
    }
  //
  // read in the grid and create output file
  //
  TFile infile(gridFileName.c_str(), "read");
  grid* g = new grid(gridFileName);
  g->trim();

  TFile outfile(outFileName.c_str(), "recreate");

  std::cout<<"Readinf grid...done"<<std::endl;

  //
  //--------------------------------------------------------------------------------
  //
  //
  //--------------------------------------------------------------------------------
  initmypdfset_( ("PDFsets/"+PDFName[0]).c_str());
  initmypdfmember_(0);
  gridObservable = (TH1D*)g->convolute(GetPdf, alphaspdf_);
  reference = (TH1D*)infile.Get("grid/reference");
  ratioGridObservable = (TH1D*)gridObservable->Clone("ratioGridObservable");
  ratioGridObservable->Divide(reference);
  removeErrors(ratioGridObservable);
  gridObservable->SetName("applgrid");
  reference->SetName("nlojet");
  ratioGridObservable->SetName("ratio");
  //
  //
  //--------------------------------------------------------------------------------
  int minBin = 1;

  for (int iPDF = 1; iPDF < nPDFs ;iPDF++)
    {
      string PDFset = "PDFsets/";
      PDFset += PDFName[iPDF];
      initmypdfset_( PDFset.c_str());

      for (int iPDFMember = 0; iPDFMember < nPDFMembers[iPDF] ;iPDFMember++)
	{
	  //  initialize pdfs
	  //	  initmypdf_( PDFset.c_str(), iPDFMember);
	  initmypdfmember_(iPDFMember);
	  // observable
	  sigma[iPDF][iPDFMember] = (TH1D*)g->convolute(GetPdf, alphaspdf_);
      
	  string sigmaName = "sigma_" + PDFName[iPDF] + "_" + to_String(iPDFMember);
	  sigma[iPDF][iPDFMember]->SetName(sigmaName.c_str());

	  string ratioName =  "ratio_" + PDFName[iPDF] + "_" + to_String(iPDFMember);
	  ratio[iPDF][iPDFMember] = (TH1D*)sigma[iPDF][iPDFMember]->Clone(ratioName.c_str());
	  ratio[iPDF][iPDFMember]->Divide(sigma[iPDF][0]);
	  removeErrors(ratio[iPDF][iPDFMember]);

	  string sigmaTitle = "#sigma ( " + to_String(PDFName[iPDF]) + "( member = " + to_String(iPDFMember) + "))";
	  sigma[iPDF][iPDFMember]->SetTitle(sigmaTitle.c_str());
	  string ratioTitle = "#frac{#sigma( " + to_String(PDFName[iPDF]) + "( member = " + to_String(iPDFMember) + "))}"+ "{#sigma_{0}( " + to_String(PDFName[0]) + " member = " + to_String(0) + ")}";
	  ratio[iPDF][iPDFMember]->SetTitle(ratioTitle.c_str());
  	}

      initmypdfmember_(0);
      for (int iScale = 0; iScale < nScales; iScale++)
	{
	  sigmaScale[iPDF][iScale] = (TH1D*)g->convolute(GetPdf, alphaspdf_, nLoops, muR[iScale], muF[iScale]);
	  string sigmaScaleName = "sigmaScale_" + to_String(iPDF) + "_" + to_String(iScale);
	  sigmaScale[iPDF][iScale]->SetName(sigmaScaleName.c_str());
	  
	  string ratioScaleName =  "ratioScale_" + to_String(iPDF) + "_" + to_String(iScale);
	  ratioScale[iPDF][iScale] = (TH1D*)sigmaScale[iPDF][iScale]->Clone(ratioScaleName.c_str());
	  ratioScale[iPDF][iScale]->Divide(sigma[iPDF][0]);
	}
    }
  //
  // writing out results
  //
  outfile.cd();

  gridObservable->Write();
  reference->Write();
  ratioGridObservable->Write();

  for (int iPDF = 1; iPDF < nPDFs ;iPDF++)
    {
      Directory currentPDF(PDFName[iPDF]);
      currentPDF.push();

      for (int iPDFMember = 0; iPDFMember < nPDFMembers[iPDF] ;iPDFMember++)
	{
	  sigma[iPDF][iPDFMember]->Write();
	  ratio[iPDF][iPDFMember]->Write();
	}
      for (int iScale = 0; iScale < nScales; iScale++)
	{
	  sigmaScale[iPDF][iScale]->Write();
	  ratioScale[iPDF][iScale]->Write();
	}
      currentPDF.pop();
    }

  outfile.Close();
  infile.Close();
  //
  //
  //  
  delete g;
  return 0;
}








