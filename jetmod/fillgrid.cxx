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


// THIS is already defined in Makefile
//#define USELHAPDF          // if defined use PDF from LHAPDF, if not use NLOJET PDF interface
//#define ADDREFERENCE       // add reference histograms per default only m_obs_bins is filled
//#define ADD_REN_REFERENCE  // add additional reference histos
//#define CONTACT            // Contact interactions 

#include <algorithm>
#include <kT_clus.h>

#ifdef USELHAPDF 
#include "pdf-genlha.h"
#else
#include "pdf-cteq6.h"
#endif

#include <string.h>
#include <stdio.h>
#include <time.h>
#include <bits/phys-cone_seedless.h>

#include "nlogrid.h"
#include <stdlib.h>
//
//
//
//const  long int saveAfterEvents = 1001;
//const long int maxNumberEvents = 1000;
// GPS tmp reduction
const  long int saveAfterEvents = 10001;
//const  long int saveAfterEvents = 100001;
const long int maxNumberEvents = 1000000;

const int debug = 0;
const double pb_fac = 3.89379656e8 ;    // conversion GeV^2 -> pb  
const double pb_fac_Unity = 1.0;        // unit factor

//const int nGrids = 2;
//std::string gridName[nGrids] = {"weight_c", "weight_f"};
const int nGrids = 1;
std::string gridName[nGrids] = {"weight_c"};
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
  
  static const int numberOfJets = 1;// MIN number of jets in event
  static const int sqrts = 14000;   // 14 TeV initial energy
  // static const int sqrts = 10000;   // 10 TeV initial energy

  double x1;                        // initial parton1 momentum fraction
  double x2;                        // initial parton2 momentum fraction
  double SCALE2[nGrids];            // highest PT^2 of jets in event
  // cut for events
  double pT_min[nGrids], eta_min[nGrids], eta_max[nGrids];  // min PT and eta range of resolved jets
  //
  bool eventSelected[nGrids];       // event has passed the cuts
  nlogrid* mygrid[nGrids];          // grid for coeficients
  
private:

  cone_seedless* jetclus;           // jet clustering algorithm
  void clusterJets(const event_hhc&);
  

#ifdef USELHAPDF                    // PDF part
  pdf_and_coupling_hhc * pdf;       // LHA
#else
  pdf_cteq6 pdf;                    // ZN
#endif
  
  
  typedef lorentzvector<double> _Lv;
  bounded_vector<_Lv> cj, pj[nGrids];                         // the jet structure in lab. frame
  bounded_vector<unsigned int> jet;
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
  cout<<" \t\t\t UserHHC destructor..."<<endl;
  delete jetclus;
  cout<<"\t\t\t\t jet clustering algorythm deleted"<<endl;
#ifdef USELHAPDF
  delete pdf;
  cout<<"\t\t\t\t pdf deleted"<<endl;
#endif
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
  jetclus=new cone_seedless(0,1);

  pT_min[0]  = 100.0;                       //
  eta_min[0] = 0.0;                         // interesting jet parameters 
  eta_max[0] = 1.0;                         //  fro observable #1

  if (nGrids>1) {
   pT_min[1]  = 100.0;                       //
   eta_min[1] = 2.0;                         // interesting jet parameters 
   eta_max[1] = 3.0;                         // for observable #2
  } 
 
  // create grid structure
  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      mygrid[iGrid] = new nlogrid(gridName[iGrid]);
      mygrid[iGrid]->getGridObject()->setCMSScale(sqrts);
    }
  //Zero the counter for number of events
  eventNb = 0;
  if (debug) 
    {
      cout<<"UserHHC::UserHHC() \t\t gDirectory = ";
      gDirectory->pwd();
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
  nj = UserHHC::numberOfJets;
  // number of UP quarks
  nu = 2U;
  // number of DOWN quarks
  nd = 3U;
}

void psinput(phasespace_hhc *ps, double& s)
{
  //  total c.m. energy square
  s = UserHHC::sqrts*UserHHC::sqrts;                  // unit:GeV
  
  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 


#include <stdlib.h>

void UserHHC::initfunc(unsigned int)
{
#ifdef USELHAPDF 

  //  system("setenv PDFSetPath `lhapdf-config --pdfsets-path`");
  //  std::string PDFSetPath = getenv("PDFSetPath");

  //  std::cout << "PDFSetPath " << PDFSetPath << std::endl;;

  //  pdf = new pdf_genlha("PDFsets/cteq6mE.LHgrid",0);
  pdf = new pdf_genlha("PDFsets/cteq6mE.LHgrid",0);
#endif
}

void UserHHC::clusterJets(const event_hhc& p)
{
  //  std::cout << "UserHHC::clusterJets(const event_hhc& p)" << std::endl;

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
  
  if(cj.upper() < numberOfJets) 
    { 
      return;                                         //  unable to resolve
    }
  std::sort(cj.begin(), cj.end(), pT_sort());
  //
  //  SELECT ONLY INTERESTING JETS
  //
  for(int i = 1; i <= cj.upper(); i++)
    {
      for (int iGrid = 0; iGrid < nGrids; iGrid++)
	{
	  double candidateRap = abs(cj[i].rapidity());
	  double candidatePT = cj[i].perp();
          if (debug)cout<<iGrid<<" eta_max= "<< eta_max[iGrid]
                               <<" eta_min= "<< eta_min[iGrid]
			<<" pT_min= "<<pT_min[iGrid]<<endl;

          if (debug)cout<<"  candidateRap= "<< candidateRap<<"  candidatePT= "<<candidatePT<<endl;
	  if(( candidateRap < eta_max[iGrid]) && (candidateRap > eta_min[iGrid]) && ( candidatePT > pT_min[iGrid]) ) 
	    {
	      (pj[iGrid]).push_back(cj[i]);
	    }
	}
    }

  if (debug) {
    cout<<" after selection of jets "<<endl;
      for (int iGrid = 0; iGrid < nGrids; iGrid++)
	{
	  for(int i = 1; i <= cj.upper(); i++)
	    {
	      cout<<"cj[ "<<i<<" ]= "<<cj[i]<<endl;
	    }
	  for(int i = 1; i <= (pj[iGrid]).upper(); i++)
	    {
	      cout<<"pj["<<iGrid<<"][ "<<i<<" ]= "<<(pj[iGrid])[i]<<endl;
	    }
	  cout<<"x1 = "<<x1<<" x2 = "<<x2<<" Scale2 = "<<SCALE2[iGrid]<<endl;
	}
      cout<<"---------------------------------"<<endl;
  }

  bool goodEvent = false;
  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      std::sort((pj[iGrid]).begin(), (pj[iGrid]).end(), pT_sort());
      if ((pj[iGrid]).upper() >= numberOfJets)
	{
	  eventSelected[iGrid] = true;
	  goodEvent = true;
	}
    }
 
  if (debug) {
   cout<< " goodEvent= "<<goodEvent<<endl;
   for (int iGrid = 0; iGrid < nGrids; iGrid++) {
     cout<<iGrid<<" eventSelected= "<<eventSelected[iGrid]<<endl;
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
  if (debug)
    {
      cout<<goodEvent<<"\t--------------------"<<endl;
      for (int iGrid = 0; iGrid < nGrids; iGrid++)
	{
	  for(int i = 1; i <= cj.upper(); i++)
	    {
	      cout<<"cj[ "<<i<<" ]= "<<cj[i]<<endl;
	    }
	  for(int i = 1; i <= (pj[iGrid]).upper(); i++)
	    {
	      cout<<"pj["<<iGrid<<"][ "<<i<<" ]= "<<(pj[iGrid])[i]<<endl;
	    }
	  cout<<"x1 = "<<x1<<" x2 = "<<x2<<" Scale2 = "<<SCALE2[iGrid]<<endl;
	}
      cout<<"---------------------------------"<<endl;
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
#if 0
  static int cock=0;

  if ( cock>9900 ) 
    std::cout << "UserHHC::userfunc()" << std::endl;


  if ( cock%100==0 ) std::cout << "UserHHC::userfunc() " << cock << std::endl;
 
  cock++;
#endif

  int iorder = amp.contrib()==0 ? 0 : 1 ;     // born contrib = 0 ; NLO = 1
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
	      long int simEvents = eventNb - 1;
	      mygrid[iGrid]->writeGrid(simEvents);
	    }
	}

      if(eventNb == (maxNumberEvents + 1))
	{
	  cout<<"Maximum number of events = "<<eventNb<<" has been simulated..."<<endl;
	  for (int iGrid = 0; iGrid < nGrids; iGrid++)
	    {
	      long int simEvents = eventNb - 1;
	      mygrid[iGrid]->writeGrid(simEvents);
	    }
	  exit(0);
	}
    }
  if (debug) cout<<" evtNb = "<<eventNb<<"\t"<<p<<endl;

  //
  // find jets in event satisfying all the cuts
  // and calculated amplitudes
  //
  clusterJets(p);

  amp.pdf_and_qcd_coupling(0, pb_fac_Unity);
  
  double weight[nGrids][7];

  for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
      if (eventSelected[iGrid])
	{
	  amp(SCALE2[iGrid], SCALE2[iGrid]);
	  for (int iSubProcess = 0; iSubProcess <= 6; iSubProcess++)
	    {
	      weight[iGrid][iSubProcess] = amp(SCALE2[iGrid], SCALE2[iGrid])[iSubProcess];
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
	  if( mygrid[iGrid]->getMode() == 0 )
	    {
	      //find interesting phase space   
	      // called after jets have been resolved
	      // no calculations ONLY fills interesting phase space with weight 1
	      if (debug) cout<<"\t GRID( "<<iGrid<<" ) phaseSpace is filling..."<<endl;
	      for(int jet = 1; jet <= (pj[iGrid]).upper(); jet++) 
		{ //loop over jets
		  double pTjet  = ((pj[iGrid])[jet]).perp();
		  mygrid[iGrid]->fillPhaseSpace(x1, x2, SCALE2[iGrid], pTjet, weight[iGrid], iorder);
		  
		  if (debug) cout <<"\t\t GRID( "<<iGrid<<" )\t jet # = "<<jet<<": pt = "<<pTjet<<" x1 = "<<x1<<" x2 = "<<x2<<" Scale = "<<SCALE2[iGrid]<<endl;
		}
	    }
	  else
	    { 
	      // calculate the weights for the grid
	      // and fill grid with weights
	      amp.pdf_and_qcd_coupling(pdf, pb_fac);
	      //TC not needed here:  amp(SCALE2[iGrid], SCALE2[iGrid]);
              //done in fillReferenceHistogram
	      if (debug) cout<<"\t GRID( "<<iGrid<<" ) weight filling..."<<endl;	      
	      for(int jet = 1; jet <= (pj[iGrid]).upper(); jet++) 
		{ //loop over jets
		  double pTjet  = ((pj[iGrid])[jet]).perp();
		  if (debug) cout <<"\t\t GRID( "<<iGrid
                                  <<" )\t jet # = "<<jet<<": pt = "<<pTjet
                                  <<" x1 = "<<x1<<" x2 = "<<x2
                                  <<" Scale = "<<SCALE2[iGrid]<<endl;

		  mygrid[iGrid]->fillWeights(x1, x2, SCALE2[iGrid], pTjet, weight[iGrid], iorder);
		  mygrid[iGrid]->fillReferenceHistograms(iorder, SCALE2[iGrid], pTjet, amp);
		}
	    }         // else
	}             // if (eventSelected)
    }                 // for (;;iGrid++)

#if 0
  if ( cock>9900 ) 
    std::cout << "UserHHC::userfunc() returning" << std::endl;
#endif

}                    // end of UserHHC::userfunc
