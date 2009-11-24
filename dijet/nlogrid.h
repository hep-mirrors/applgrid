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


class nlogrid
{
 public:
  nlogrid(std::string inputName = "test");
  ~nlogrid();

 public:

  inline int getMode(){return mode;};

  TH1D* getGridReference(){return gridObject->getReference();};

  void writeGrid(long int);
  void writeHistToGridFile(TH1D*);
  
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
  
  void fillReferenceHistograms(const int &iorder,
			       const double &SCALE2,
			       const double &obs, 
			       const nlo::amplitude_hhc& amp
			       );
  
  
 public:

  static const int numberOfJets;

 private:
  appl::grid* gridObject;
  std::string fullFileName;      // name of the file to store grid
  std::string observableName;    // name of the file to store grid
  int mode;                      // mode of grid (0 = new ; 1 = optimisation)
  long int numOfEvents;

  std::string refDirName;

};



#endif
