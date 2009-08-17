#include <sys/stat.h>

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
#include "TVectorT.h"

//
//     renormalisation and factorisation scales
//
static const double mur[Nscales] = {1.0,0.5,2.0,1.0,0.5};
static const double muf[Nscales] = {1.0,0.5,2.0,0.5,1.0};
// static const double mur[Nscales] = {1.0,0.5,2.0,1.0,1.0};
// static const double muf[Nscales] = {1.0,1.0,1.0,0.5,2.0};
//
//grid parameters
//
static const int NQ2bins = 10, iOrderQ2 = 5;
static const double Q2low = 1.0e4, Q2up = 4.9e7;

static const int  Nxbins = 30,  iOrderx = 5;
static const double xlow = 1e-5, xup = 1.0;

static const double apramval=5.;
static const bool pdfWeight = false;

static const string pdf_function = "nlojet";

static const int lowest_order = 2, nloops = 1;
//    
//observable
//
const int nobs=24;                 //no of bins
double xbins[nobs+1]=
  {                                //bin-information, lower edges...
    100., 200., 300., 400., 500., 
    600., 700., 800., 900.,1000.,
    1100.,1200.,1300.,1400.,1500.,
    1600.,1700.,1800.,1900.,2000.,
    2500.,3000.,3500.,4000.,4500.
  };
//
//
//
//  Constructor & destructor
//
//
//
nlogrid::nlogrid(std::string inputName)
{
  // number of events
  numOfEvents = 0;
  // file to store grid
  observableName = inputName;
  fullFileName = "./output/" +  observableName + ".root";

  refDirName = "addReference";

  Directory nlogridDir(observableName);
  nlogridDir.push();

  //  FILE * testfile; 
  //  testfile = fopen (fullFileName.c_str(),"r");
  //  if (testfile == NULL) 
  
  // test if file exists, don't need to try to open it, 
  struct stat stFileInfo;
  if ( stat(fullFileName.c_str(),&stFileInfo) )   
    {
      cout<<"Creating new grid... "<<endl;
      
      mode=0;
      

      appl::igrid::transformvar(apramval);
      appl::igrid::reweight(pdfWeight);
      

      gridObject = 
	new appl::grid(
		       nobs,xbins,
		       NQ2bins, Q2low, Q2up, iOrderQ2,    
		       Nxbins, xlow, xup, iOrderx,       
		       pdf_function, lowest_order, nloops
		       );
      cout<<"DONE Creating new grid"<<endl;
    }
  else
    {
      cout<<"Creating optimized grid... "<<endl;
      
      mode=1;

      gridObject = new appl::grid(fullFileName);
      nlogridDir.push();
      if (gridObject->isOptimised())
	{
	  delete gridObject;
	  std::cout <<"Grid is aready optimised. Quitting ..."<<std::endl;
	  exit(0); // should really throw an exception
	}
      // reseting reference histgram
      TH1D* htemp = gridObject->getReference();
      for (int i = 0; i <= htemp->GetNbinsX() + 1; i++) htemp->SetBinContent(i,0);

      gridObject->optimise();

      cout<<"DONE Creating optimized grid"<<endl;
    }
  
#ifdef REN_REFERENCE
  // booking reference histos
  bookReferenceHistograms();  
#endif
  //Write out grid information to screen
  //cout<<*gridObject<<endl;
  nlogridDir.pop();
}
//
// destructor
//
nlogrid::~nlogrid()
{
  cout<<"\t\t\t\t nlogrid destructor"<<endl;

#ifdef REN_REFERENCE
  deleteReferenceHistograms();
#endif

  delete gridObject;

  cout<<" \n\n\n\n\n\n \t\t\t Calculation finished!!!! \n\n\n\n\n\n"<<endl;
}
//
//    Grid persistency
//
void nlogrid::writeGrid(long int& nRuns)
{
  gridObject->run() = nRuns;

  Directory obs(observableName);
  obs.push();
  gridObject->Write(fullFileName);

  obs.pop();

#ifdef REN_REFERENCE
  //  if (mode)
    { 
      writeReferenceHistograms();
    }
#endif

  string goMode = (mode == 0 ? "Non-o": "O");
  goMode += "ptimised";

  std::cout<<"\tGridObject ( "<<goMode<<" ) saved after "<<nRuns
	   <<" number of events."
	   <<"\n Weightgrid File = "<<fullFileName<<" ."
	   <<std::endl;
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
//
//
//
//  Reference histograms manipulations
//
//
//
//
#ifdef REN_REFERENCE

void nlogrid::bookReferenceHistograms()
{
  Directory refDir(refDirName);
  refDir.push();
  
  char histname[3000]; char htit[3000];
  for(int is = 0; is < 7; is++) 
    {// loop over subprocesses
      sprintf(histname,"soft_sub_%d",is);
      sprintf(htit,"addRefHist. subProcess = %d",is);
      soft_sub[is] = new TH1D(histname, htit, nobs, xbins);
      soft_sub[is]->Sumw2();
      for(int ir = 0; ir <  Nscales; ir++)
	{ // loop ren scal variations
	      sprintf(histname, "soft_subscale_%d_%d", is, ir);
	      sprintf(htit,"addRefHist. scales: #mu_{R} = %3.2f,  #mu_{F} = %3.2f, subProcess = %d ", mur[ir], muf[ir], is);
	      soft_subscale[is][ir] = new TH1D(histname, htit, nobs, xbins);
	      soft_subscale[is][ir]->Sumw2();
	      if (is == 0)
		{
		  sprintf(histname,"soft_scale_%d",ir);
		  sprintf(htit,"addRefHist. scales  #mu_{R} = %3.2f, #mu_{F} = %3.2f", mur[ir], muf[ir]);
		  soft_scale[ir] = new TH1D(histname, htit, nobs, xbins);
		  soft_scale[ir]->Sumw2();
		}
	}           // ir
    }               // is

  refDir.pop();

}                   // bookReferenceHistograms

void nlogrid::deleteReferenceHistograms()
{
  
  for(int is = 0; is < 7; is++) 
    {  
      delete soft_sub[is];
      soft_sub[is] = NULL;
      for(int ir = 0; ir < Nscales; ir++) 
	{ 
	  delete soft_subscale[is][ir];
	  soft_subscale[is][ir] = NULL;
	  
	  if(is == 0)
	    {
	      delete soft_scale[ir];
	      soft_scale[ir] = NULL;
	    }
	}
    }
}

void nlogrid::writeReferenceHistograms()
{
  Directory obs(observableName);
  obs.push();

  //cout << "nlogrid::writeReferenceHistograms() \t\t pwd=" << gDirectory->GetName() << endl;

  TFile myfile(fullFileName.c_str(),"UPDATE");

  // write the scales to the histogram;

  TVectorT<double> mur_v(Nscales);
  TVectorT<double> muf_v(Nscales);
  
  for ( int i=0 ; i<Nscales ; i++ ) { 
    mur_v(i) = mur[i];
    muf_v(i) = muf[i];
  }
  
  //  mur_v.Write("mu_r");
  //  muf_v.Write("mu_f");
  
  mur_v.Write("mu_r");
  muf_v.Write("mu_f");

  Directory refDir(refDirName);
  refDir.push();
  //  cout << "nlogrid::writeReferenceHistograms() \t\t pwd=" << gDirectory->GetName() << endl;  
  
  for(int is=0; is<7; is++) 
    {
      soft_sub[is]->Write("",TObject::kOverwrite);
      for(int ir=0; ir< Nscales; ir++)
	{
	  soft_subscale[is][ir]->Write("",TObject::kOverwrite);
	  if(is == 0)
	    {
	      soft_scale[ir]->Write("",TObject::kOverwrite);   
	    }
	}
    }
  refDir.pop();
  myfile.Close();

  obs.pop();
}

#endif 

void nlogrid::fillReferenceHistograms(const int &iorder,
				      const double &SCALE2,
				      const double &obs, 
				      const nlo::amplitude_hhc& amp
				      )
{
  
  bool debug=false;
  double binwidth = gridObject->deltaobs(gridObject->obsbin(obs));
  
  
  amp(SCALE2,SCALE2);
  for(int is = 0; is < 7; is++) {
   double weightval = amp(SCALE2,SCALE2)[is]/binwidth;
   gridObject->getReference()->Fill(obs, weightval);
   soft_sub[is]->Fill(obs, weightval);
  }
      
#ifdef REN_REFERENCE
  for(int ir = 0; ir < Nscales; ir++){
    //amp(mur[ir]*mur[ir]*SCALE2, muf[ir]*muf[ir]*SCALE2);// needed ?
    if (debug) cout<<amp.contrib()<<" obs= "<<obs
                   <<" ...filling reference "
                   <<" mur["<<ir<<"]= "<<mur[ir]
                   <<" muf["<<ir<<"]= "<<muf[ir]<<endl;

   for(int is = 0; is < 7; is++){
    double weightval = amp(mur[ir]*mur[ir]*SCALE2, muf[ir]*muf[ir]*SCALE2)[is]/binwidth;
    if (debug) cout<<" is= "<<is
                   <<" weight= "<< weightval<<endl;

    soft_subscale[is][ir]->Fill(obs, weightval);
    soft_scale       [ir]->Fill(obs, weightval);
   }
  }
#endif 
}




