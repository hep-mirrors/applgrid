//
//
//
//------ DON'T TOUCH THIS PART! ------//------ DON'T TOUCH THIS PART! ------
//------ DON'T TOUCH THIS PART! ------//------ DON'T TOUCH THIS PART! ------
//------ DON'T TOUCH THIS PART! ------//------ DON'T TOUCH THIS PART! ------
//
//
//
#include <bits/hhc-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>
//#include <nlo++-module_add.h>
//----- used namespaces -----
using namespace nlo;
using namespace std;
//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_hhc *, double&);
user_base_hhc * userfunc();

//typedef unsigned long int (*module_add_type)(bool, const list<basic_string<char> >&, const basic_string<char>&);
//extern  module_add_type module_add;


//----- array of the symbols symbols -----
extern "C"{
struct { 
  const char *name;
  void *address;
} user_defined_functions[] = 
  {
    //   process index: hhc for hadron-hadron --> jets
    {"procindex", (void *) "hhc"},
    
    //   input function 
    {"inputfunc", (void *) inputfunc},
    
    //   phase space input function 
    {"psinput", (void *) psinput},
    
    //   user defined functions
    {"userfunc",  (void *) userfunc},
    
    //   module to generate the readable result
    //    {"main_module_add", (void *) module_add},
    
    //  end of the list
    {0, 0}
  };
};
//
//
//------ USER DEFINED PART STARTS HERE ------
//------ USER DEFINED PART STARTS HERE ------
//------ USER DEFINED PART STARTS HERE ------
//
//
#include <algorithm>
#include <kT_clus.h>

#include "pdf-genlha.h"

#include <string.h>
#include <stdio.h>
#include <time.h>
#include <bits/phys-cone_seedless.h>

#include "nlogrid.h"
#include "fjClustering.h"
#include <stdlib.h>

#include <TH1D.h>
//
//
//
//const long int saveAfterEvents = 100001;
//const long int maxNumberEvents = 10000;
const  long int saveAfterEvents = 10001;
const long int maxNumberEvents =  10000;

const int debug = 0;
const double pb_fac = 3.89379656e8 ;                          // conversion GeV^2 -> pb  

static const int nGrids = 4;
const std::string obsName[nGrids] = {"mass","pt","phi","phiFine"};
const std::string gridNamePrefix[5] = {"zero","single","di","tri","quatro"};

//static const int sqrts =  1800;                               // initial energy
//  static const int sqrts =  5000;                             // initial energy
//  static const int sqrts =  7000;                             // initial energy
static const int sqrts = 10000;                               // initial energy
//  static const int sqrts = 14000;                             // initial energy
//  static const int sqrts = 16000;                             // initial energy
  //  static const int sqrts = 18000;                           // initial energy

#include <sstream>

template <class T>
inline std::string to_String (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}


//
//   user class
//
class UserHHC : public user1d_hhc
{
public:
  UserHHC();
  ~UserHHC();

  long int eventNb;                    // number of event

  void initfunc(unsigned int);
  void userfunc(const event_hhc&, const amplitude_hhc&);
  

  double x1;                        // initial parton1 momentum fraction
  double x2;                        // initial parton2 momentum fraction
  double SCALE2[nGrids];            // highest PT^2 of jets in event
  // cut for events
  double pT_cut[20][nGrids], eta_cut[20][nGrids];  // min PT and eta range of resolved jets
  //
  long int selectedEventsNb[nGrids];
  bool eventSelected[nGrids];       // event has passed the cuts
  nlogrid* mygrid[nGrids];          // grid for coeficients
  
private:

  //  cone_seedless *jetclus;           // jet clustering algorithm
  fjClustering *jetclus;           // fastjet jet clustering algorithms
  void clusterJets(const event_hhc&);
  

  pdf_genlha *pdf;                    
  
  std::string gridName[nGrids];

  
  typedef lorentzvector<double> _Lv;
  bounded_vector<_Lv> cj, pj[nGrids];                         // the jet structure in lab. frame
  bounded_vector<unsigned int> jetjet;
  //////  void fill_jets(int, double, double , const amplitude_hhc&);
  struct pT_sort {
    bool operator()(const _Lv& p1, const _Lv& p2) const {
      return p1.perp2() > p2.perp2();
    }
  };
  struct E_sort {
    bool operator()(const _Lv& p1, const _Lv& p2) const {
      return p1.T() > p2.T();
    }
  };
  
};
//
//  destructor
//
UserHHC::~UserHHC() 
{
  delete jetclus;
  delete pdf;
  delete [] mygrid;
  cout<<"Calculation finished."<<endl;
}
//
//  constructor
//
UserHHC::UserHHC() : pdf(0)  
{
  cout<<"NLOJET++ FillGrid started..."<<endl;
  // create jet algorythm
  //  jetclus = new cone_seedless(0,1);
//  jetclus = new fjClustering();
  jetclus = new fjClustering(fastjet::antikt_algorithm, 0.6, fastjet::E_scheme, fastjet::Best);
  //---------------------------------------------------------------------------------------
  /* pt */
  pT_cut[0][0] = 100.0;                        //
  pT_cut[1][0] = 100.0;                        //  interesting jet parameters for mass distr
  if (nlogrid::numberOfJets == 3)              //
    pT_cut[2][0] = 20.0;                       // 
  if (nlogrid::numberOfJets == 4)              //
    pT_cut[3][0] = 20.0;                       // 
  /* eta */
  eta_cut[0][0] = 0.50;                        //
  eta_cut[1][0] = 0.50;                        //
  if (nlogrid::numberOfJets == 3)              //
    eta_cut[2][0] = 0.50;                      //
  if (nlogrid::numberOfJets == 4)              //
    eta_cut[3][0] = 0.50;                      //
  //---------------------------------------------------------------------------------------
  /* pt */
  pT_cut[0][1] = 100.0;                        //
  pT_cut[1][1] = 100.0;                        //  interesting jet parameters for pt_{njet} distr
  if (nlogrid::numberOfJets == 3)              // 
    pT_cut[2][1] = 20.0;                       //
  if (nlogrid::numberOfJets == 4)              // 
    pT_cut[3][1] = 20.0;                       //
  /* eta */
  eta_cut[0][1] = 0.50;                        //
  eta_cut[1][1] = 0.50;                        //
  if (nlogrid::numberOfJets == 3)              //
    eta_cut[2][1] = 0.50;                      //
  if (nlogrid::numberOfJets == 4)              //
    eta_cut[3][1] = 0.50;                      //
  //---------------------------------------------------------------------------------------
  /* pt */
  pT_cut[0][2] = 100.0;                        //
  pT_cut[1][2] = 100.0;                        //  interesting jet parameters for deltaphi distr
  if (nlogrid::numberOfJets == 3)              // 
    pT_cut[2][2] = 20.0;                       //
  if (nlogrid::numberOfJets == 4)              // 
    pT_cut[3][2] = 20.0;                       //
  /* eta */
  eta_cut[0][2] = 0.50;                        //
  eta_cut[1][2] = 0.50;                        //
  if (nlogrid::numberOfJets == 3)              //
    eta_cut[2][2] = 0.50;                      //
  if (nlogrid::numberOfJets == 4)              //
    eta_cut[3][2] = 0.50;                      //
  //---------------------------------------------------------------------------------------
  /* pt */
  pT_cut[0][3] = 100.0;                        //
  pT_cut[1][3] = 100.0;                        //  interesting jet parameters for deltaphiFine distr
  if (nlogrid::numberOfJets == 3)              // 
    pT_cut[2][3] = 20.0;                       //
  if (nlogrid::numberOfJets == 4)              // 
    pT_cut[3][3] = 20.0;                       //
  /* eta */
  eta_cut[0][3] = 0.50;                        //
  eta_cut[1][3] = 0.50;                        //
  if (nlogrid::numberOfJets == 3)              //
    eta_cut[2][3] = 0.50;                      //
  if (nlogrid::numberOfJets == 4)              //
    eta_cut[3][3] = 0.50;                      //

  


 
  // create grid structure
  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      gridName[iGrid] = gridNamePrefix[nlogrid::numberOfJets] + "jet_" + obsName[iGrid] + "_" + to_String(nlogrid::numberOfJets);
      mygrid[iGrid] = new nlogrid(gridName[iGrid]);
      selectedEventsNb[iGrid] = 0;
    }
  //Zero the counter for number of events
  eventNb = 0;
  if (debug) 
    {
      cout<<"UserHHC::UserHHC() \t\t gDirectory = ";
      gDirectory->pwd();

      for (int iGrid = 0; iGrid < nGrids; iGrid++)
	{
	  mygrid[iGrid]->getGridReference()->Print("all");
	}
    }

}
//
//
//
user_base_hhc* userfunc() {return new UserHHC;}
//
//
//
void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  nj = nlogrid::numberOfJets;
  // number of UP quarks
  nu = 2U;
  // number of DOWN quarks
  nd = 3U;
}
void psinput(phasespace_hhc *ps, double& s)
{
  //  total c.m. energy square
  s = sqrts * sqrts;                  // unit:GeV
  
  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 

void UserHHC::initfunc(unsigned int)
{
  pdf = new pdf_genlha("PDFsets/cteq6mE.LHgrid",0);
}

void UserHHC::clusterJets(const event_hhc& p)
{
  //
  // initialisation
  //
  x1 = x2 = 0.;
  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      SCALE2[iGrid] = 0.;
      eventSelected[iGrid] = false;
      (pj[iGrid]).resize(1,0);
    }
  //
  //----- do the cluster analysis-----
  //

  cj=jetclus->operator()(p); 
  
  if(cj.upper() < nlogrid::numberOfJets) 
    { 
      return;                                         //  unable to resolve
    }
  std::sort(cj.begin(), cj.end(), pT_sort());
  //
  //  SELECT ONLY INTERESTING JETS
  //
  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      for (int iJet = 1; iJet <= nlogrid::numberOfJets; iJet++)
	{
	  if ((std::abs(cj[iJet].rapidity()) < eta_cut[iJet - 1][iGrid]) &&  (cj[iJet].perp() > pT_cut[iJet - 1][iGrid])) 
	    {
	      (pj[iGrid]).push_back(cj[iJet]);
	    }
	}
    }

  bool goodEvent = false;
  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      std::sort((pj[iGrid]).begin(), (pj[iGrid]).end(), pT_sort());
      if ((pj[iGrid]).upper() >= nlogrid::numberOfJets)
	{
	  eventSelected[iGrid] = true;
	  goodEvent = true;
	}
    }
  
  if (goodEvent)
    {
      x1 = p[-1].T()/(0.5*sqrts);
      x2 =  p[0].T()/(0.5*sqrts);

      for (int iGrid = 0; iGrid < nGrids; iGrid++)
	{
	  if (eventSelected[iGrid]) SCALE2[iGrid] = ((pj[iGrid])[1]).perp2();
	}
    }
}
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//
//   User analysis
//
void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{ 
  int iorder = amp.contrib()==0 ? 0 : 1 ;     // born contrib = 0 ; NLO = 1

  //  std::cout <<"\t\t !!!! \t\t"<<amp.contrib()<<std::endl;
  if( iorder == 0) 
    {
      eventNb++;
      // persistency
      // do some things every saveAfterEvents
      if(eventNb % saveAfterEvents == 1 && !( eventNb == 1)) 
	{
	  cout<<"saving grid after "<<(eventNb-1)<<" runs"<<endl;
	  for (int iGrid = 0; iGrid < nGrids; iGrid++)
	    {
	      mygrid[iGrid]->writeGrid(eventNb - 1);
	    }
	}

      if(eventNb == (maxNumberEvents + 1))
	{
	  cout<<"Maximum number of events = "<<eventNb<<" has been simulated..."<<endl;
	  for (int iGrid = 0; iGrid < nGrids; iGrid++)
	    {
	      cout <<"Has been selected "<<selectedEventsNb[iGrid]<<" events"<<endl;
	      cout <<"grid # "<<iGrid<<" reference histogram contains "<<mygrid[iGrid]->getGridReference()->GetEntries()/7<<" number of entries per subprocess"<<endl;
	      mygrid[iGrid]->writeGrid(eventNb - 1);
	    }
	  exit(0);
	}
    }
  //
  // find jets in event satisfying all the cuts
  // and calculated amplitudes
  //
  clusterJets(p);
  
  amp.pdf_and_qcd_coupling(0, 1.000);
  
  double weight[nGrids][7];

  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      if (eventSelected[iGrid])
	{
	  if (iorder == 0) selectedEventsNb[iGrid]++;
	  amp(SCALE2[iGrid], SCALE2[iGrid]);
	  for (int iSubProcess = 0; iSubProcess <= 6; iSubProcess++)
	    {
	      weight[iGrid][iSubProcess] = amp(SCALE2[iGrid], SCALE2[iGrid])[iSubProcess];
	      weight[iGrid][iSubProcess] *= pb_fac;
	    }
	}
    }
  //
  //  grid works...
  //
  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      if ( eventSelected[iGrid])
	{ 
	  // calculate the observable
	  double observable = 0.;
	  if (gridName[iGrid].find("_mass_") != string::npos)
	    {
	      _Lv jSum = pj[iGrid][1] + pj[iGrid][2];
	      for (int i = 3; i <= nlogrid::numberOfJets; i++)
		{
		  jSum += pj[iGrid][i];
		}
	      observable = jSum.mag();
	    }
	  else if (gridName[iGrid].find("_pt_") != string::npos)
	    {
	      observable = pj[iGrid][nlogrid::numberOfJets].perp();
	    }
	  else if (gridName[iGrid].find("_phi") != string::npos)
	    {
	      if (pj[iGrid].upper() < 2)
		{
		  std::cout <<"something wrong in the code "<<std::endl;
		  std::cout <<p<<std::endl;
		  exit(0);
		}
	      observable = std::abs(pj[iGrid][1].phi() - pj[iGrid][2].phi());
	      if (observable > M_PI) observable = 2.0*M_PI - observable;
	    }
	  else 
	    {
	      observable = -10e10;
	    }

	  if( mygrid[iGrid]->getMode() == 0 )
	    {
	      //find interesting phase space   
	      // called after jets have been resolved
	      // no calculations ONLY fills interesting phase space with weight 1
	      if (debug) cout<<"\t GRID( "<<iGrid<<" ) phaseSpace is filling..."<<endl;

	      mygrid[iGrid]->fillPhaseSpace(x1, x2, SCALE2[iGrid], observable, weight[iGrid], iorder);
	      //mygrid[iGrid]->fillReferenceHistograms(iorder, SCALE2[iGrid], observable, amp);
		  
	      if (debug) 
		{
		  cout <<"\t\t reference filled GRID( "<<iGrid<<" )\t  x1 = "<<x1<<" x2 = "<<x2<<" Scale = "<<SCALE2[iGrid]<<endl;
		  cout <<"\t\t grid = "<<gridName[iGrid]<<" obs = "<<observable<<endl;
		  cout <<" total number of selected events = "<<selectedEventsNb[iGrid]<< endl;
		  mygrid[iGrid]->getGridReference()->Print();
		}
	    }
	  else
	    {
	      // calculate the weights for the grid
	      // and fill grid with weights
	      amp.pdf_and_qcd_coupling(pdf, pb_fac);
	      amp(SCALE2[iGrid], SCALE2[iGrid]);
	      
	      mygrid[iGrid]->fillWeights(x1, x2, SCALE2[iGrid], observable, weight[iGrid], iorder);
	      mygrid[iGrid]->fillReferenceHistograms(iorder, SCALE2[iGrid], observable, amp);
	      if (debug) 
		{
		  cout <<"\t\t weights filled GRID( "<<iGrid<<" )\t  x1 = "<<x1<<" x2 = "<<x2<<" Scale = "<<SCALE2[iGrid]<<endl;
		  cout <<"\t\t grid = "<<gridName[iGrid]<<" obs = "<<observable<<endl;
		  cout <<" \t\t total number of selected events = "<<selectedEventsNb[iGrid]<< endl;
		  mygrid[iGrid]->getGridReference()->Print();
		  cout<<"events in refHist = "<<mygrid[iGrid]->getGridReference()->GetEntries()/7<<endl;
		  cout<<"--------------------"<<endl;
		}
	    }         // else
	}             // if (eventSelected)
    }                 // for (;;iGrid++)
}                     // end of UserHHC::userfunc
