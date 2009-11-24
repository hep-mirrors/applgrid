#include <string>
using std::string;

#include "nlogrid.h"
#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"
using appl::grid;
using appl::igrid;

#include <bits/hhc-process.h>

#ifdef USELHAPDF 
#include "pdf-genlha.h"
#else
#include "pdf-cteq6.h"
#endif

#include "TFile.h"

static const int nloops = 1;
static const string pdf_function = "nlojet";
//
//grid parameters
//
static const int NQ2bins = 10, iOrderQ2 = 3;
static const double Q2low = 1.0e4, Q2up = 25e6;

static const int  Nxbins = 30,  iOrderx = 3;
static const double xlow = 1e-5, xup = 1.0;

static const double apramval = 5.;
static const bool pdfWeight = false;
//
//
//
const int nlogrid::numberOfJets = 3;           // MIN number of jets in event
//    
//observable
//
static const int nobsMass = 13;                 //no of bins
double xbinsMass[nobsMass + 1] =
  {                                      //bin-information, lower edges...
    300., 349., 405., 469., 542., 
    625., 720., 827., 949.,1087.,
    1243.,1420.,1620.,2500.
  };

static const int nobsPT = 24;                 //no of bins
double xbinsPT[nobsPT + 1] =
  {                                      //bin-information, lower edges...
    10., 12., 14., 16., 18., 20., 24., 28., 32., 36., 40., 45., 50., 55., 60., 65., 70., 75.,
    100., 120., 140., 160., 180.,200., 500.
  };

static const int nobsDPhi = 30;                 //no of bins
double dp = 0.05, sp = M_PI - dp;
double xbinsDPhiFine[nobsDPhi + 1] =
  {                                      //bin-information, lower edges...
    sp, 
    sp + 1./nobsDPhi*2.*dp, sp + 2./nobsDPhi*2.*dp,sp + 3./nobsDPhi*2.*dp,sp + 4./nobsDPhi*2.*dp,
    sp + 5./nobsDPhi*2.*dp, sp + 6./nobsDPhi*2.*dp,sp + 7./nobsDPhi*2.*dp,sp + 8./nobsDPhi*2.*dp,
    sp + 9./nobsDPhi*2.*dp, sp + 10./nobsDPhi*2.*dp,sp + 11./nobsDPhi*2.*dp,sp + 12./nobsDPhi*2.*dp,
    sp + 13./nobsDPhi*2.*dp, sp + 14./nobsDPhi*2.*dp,sp + 15./nobsDPhi*2.*dp,sp + 16./nobsDPhi*2.*dp,
    sp + 17./nobsDPhi*2.*dp, sp + 18./nobsDPhi*2.*dp,sp + 19./nobsDPhi*2.*dp,sp + 20./nobsDPhi*2.*dp,
    sp + 21./nobsDPhi*2.*dp, sp + 22./nobsDPhi*2.*dp,sp + 23./nobsDPhi*2.*dp,sp + 24./nobsDPhi*2.*dp,
    sp + 25./nobsDPhi*2.*dp, sp + 26./nobsDPhi*2.*dp,sp + 27./nobsDPhi*2.*dp,sp + 28./nobsDPhi*2.*dp,
    sp + 29./nobsDPhi*2.*dp, 3.14159 + dp
  };

double SP2 = M_PI/2;
double EP2 = M_PI;
double binw = (EP2 - SP2)/nobsDPhi;
double xbinsDPhi[nobsDPhi + 1] =
  {                                      //bin-information, lower edges...
    SP2, 
    SP2 + 1.*binw, SP2 + 2.*binw, SP2 + 3.*binw, SP2 + 4.*binw,
    SP2 + 5.*binw, SP2 + 6.*binw, SP2 + 7.*binw, SP2 + 8.*binw,
    SP2 + 9.*binw, SP2 + 10.*binw, SP2 + 11.*binw, SP2 + 12.*binw,
    SP2 + 13.*binw, SP2 + 14.*binw, SP2 + 15.*binw, SP2 + 16.*binw,
    SP2 + 17.*binw, SP2 + 18.*binw, SP2 + 19.*binw, SP2 + 20.*binw,
    SP2 + 21.*binw, SP2 + 22.*binw, SP2 + 23.*binw, SP2 + 24.*binw,
    SP2 + 25.*binw, SP2 + 26.*binw, SP2 + 27.*binw, SP2 + 28.*binw,
    SP2 + 29.*binw, EP2
  };
// clone a histogram, scale it and write it, then delete it
void ScaleWrite(TH1D* h, double d ) { 
  std::string name = std::string("soft_") + std::string(h->GetName());
  TH1D* _h = (TH1D*)h->Clone(name.c_str());
  _h->Scale(d);
  _h->Write("",TObject::kOverwrite);
  delete _h;
}
//
//
//
//  Constructor & destructor
//
//
//
nlogrid::nlogrid(std::string inputName)
{
  //lowest order
  int lowest_order = 2;
  if (numberOfJets > 2) lowest_order = numberOfJets;
  // number of events
  numOfEvents = 0;
  // file to store grid
  observableName = inputName;
  fullFileName = "./output/" +  observableName + ".root";

  refDirName = "addReference";

  Directory nlogridDir(observableName);
  nlogridDir.push();

  FILE * testfile; 
  testfile = fopen (fullFileName.c_str(),"r");
  
  if (testfile == NULL) 
    {
      cout<<"Creating new grid... "<<endl;
      
      mode=0;
      

      appl::igrid::transformvar(apramval);
      appl::igrid::reweight(pdfWeight);
      
      if (inputName.find("_mass_") != string::npos)
	{
	  gridObject = 
	    new appl::grid(
			   nobsMass,xbinsMass,
			   NQ2bins, Q2low, Q2up, iOrderQ2,    
			   Nxbins, xlow, xup, iOrderx,       
			   pdf_function, lowest_order, nloops
			   );
	}
      else if (inputName.find("_pt_") != string::npos)
	{
	  gridObject = 
	    new appl::grid(
			   nobsPT,xbinsPT,
			   NQ2bins, Q2low, Q2up, iOrderQ2,    
			   Nxbins, xlow, xup, iOrderx,       
			   pdf_function, lowest_order, nloops
			   );
	}
      else if (inputName.find("_phi_") != string::npos)
	{
	  gridObject = 
	    new appl::grid(
			   nobsDPhi,xbinsDPhi,
			   NQ2bins, Q2low, Q2up, iOrderQ2,    
			   Nxbins, xlow, xup, iOrderx,       
			   pdf_function, lowest_order, nloops
			   );
	}
      else if (inputName.find("_phiFine_") != string::npos)
	{
	  gridObject = 
	    new appl::grid(
			   nobsDPhi,xbinsDPhiFine,
			   NQ2bins, Q2low, Q2up, iOrderQ2,    
			   Nxbins, xlow, xup, iOrderx,       
			   pdf_function, lowest_order, nloops
			   );
	}
      else
	{
	  std::cout <<"------------------------------------------"<< std::endl;
	  std::cout <<"\t\t\t WRONG OBSERVABLE ( " << inputName<<" )"<<std::endl;
	  std::cout <<"------------------------------------------"<< std::endl;
	  exit(0);
	}

      cout<<"DONE Creating new grid"<<endl;
    }
  else
    {
      cout<<"Creating optimized grid... "<<endl;
      
      mode=1;

      fclose(testfile);

      gridObject = new appl::grid(fullFileName);
      nlogridDir.push();
      if (gridObject->isOptimised())
	{
	  delete gridObject;
	  std::cout <<"Grid is aready optimised. Quitting ..."<<std::endl;
	  exit(0);
	}
      // reseting reference histgram
      TH1D* htemp = gridObject->getReference();
      for (int i = 0; i <= htemp->GetNbinsX() + 1; i++) htemp->SetBinContent(i,0);
      htemp->SetEntries(0);

      gridObject->optimise();

      cout<<"DONE Creating optimized grid"<<endl;
    }
  
  //Write out grid information to screen
  //cout<<*gridObject<<endl;
  nlogridDir.pop();
}
//
// destructor
//
nlogrid::~nlogrid()
{
  delete gridObject;
}
//
//    Grid persistency
//
void nlogrid::writeGrid(long int nRuns)
{
  gridObject->run() = nRuns;

  Directory obs(observableName);
  obs.push();
  gridObject->Write(fullFileName);
  obs.pop();

  string goMode = (mode == 0 ? "Non-o": "O");
  goMode += "ptimised";

  std::cout<<"\tGridObject ( "<<goMode<<" ) saved after "<<nRuns
	   <<" number of events."
	   <<"\n Weightgrid File = "<<fullFileName<<" ."
	   <<std::endl;
}
void nlogrid::writeHistToGridFile(TH1D* h)
{
    TFile myfile(fullFileName.c_str(),"UPDATE");
    h->Write();
    myfile.Close();
}
//
//
//
void nlogrid::fillPhaseSpace(const double &x1,
			     const double &x2,
			     const double &Q2,
			     const double &obs,
			     const double *weight,
			     const int    &iorder)
{
  gridObject->fill_phasespace(x1, x2, Q2, obs, weight, iorder);
}
//
//  fill weights
//
//
void nlogrid::fillWeights(const double &x1,
			     const double &x2,
			     const double &Q2,
			     const double &obs,
			     const double *weight,
			     const int    &iorder)
{
  gridObject->fill(x1, x2, Q2, obs, weight, iorder);
}
//
// fill reference
//
void nlogrid::fillReferenceHistograms(const int &iorder,
				      const double &SCALE2,
				      const double &obs, 
				      const nlo::amplitude_hhc& amp
				      )
{
  
  double binwidth = gridObject->deltaobs(gridObject->obsbin(obs));
  
  
  for(int is = 0; is < 7; is++)
    {
      double weightval = amp(SCALE2,SCALE2)[is]/binwidth;
      gridObject->getReference()->Fill(obs, weightval);
    }

}




