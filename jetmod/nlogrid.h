#ifndef NLOGRID__HH
#define NLOGRID_HH 1

#include <string>
#include "appl_grid/appl_grid.h"
#include "TH1D.h"
#include <bits/hhc-process.h>
#include <algorithm>

#ifdef USELHAPDF 
#include "pdf-genlha.h"
#else
#include "pdf-cteq6.h"
#endif

#include "scales.h"

class nlogrid
{
 public:
  nlogrid(std::string inputName = "test");
  ~nlogrid();

 public:

  inline int getMode(){return mode;};

  TH1D* getGridReference(){return gridObject->getReference();};

  void writeGrid(long int&);
  
  void fillPhaseSpace(
		      const double &x1, 
		      const double &x2, 
		      const double &SCALE2, 
		      const double &obs,
		      const double *weight,
		      const int    &iorder
		      );
  void fillWeights(
		   const double &x1, 
		   const double &x2, 
		   const double &SCALE2, 
		   const double &obs,
		   const double *weight,
		   const int    &iorder
		   );
  
#ifdef REN_REFERENCE

  void bookReferenceHistograms();
  void deleteReferenceHistograms();
  void writeReferenceHistograms();
#endif
  void fillReferenceHistograms(const int &iorder,
			       const double &SCALE2,
			       const double &obs, 
			       const nlo::amplitude_hhc& amp
			       );
  
  

 private:
  appl::grid* gridObject;
  std::string fullFileName;      // name of the file to store grid
  std::string observableName;    // name of the file to store grid
  int mode;                      // mode of grid (0 = new ; 1 = optimisation)
  long int numOfEvents;

  std::string refDirName;

 private:


#ifdef REN_REFERENCE

  // reference histograms for 2 orders and 7 subprocesses
  TH1D* soft_sub [7];
  // reference histogram for 2 orders, 7 subprocesses and  scale variations
  TH1D* soft_subscale [7][Nscales];            
  // reference histograms to test renormalisation scale dependence
  TH1D* soft_scale [Nscales];

#endif

};



#endif
