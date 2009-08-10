//   
//   fappl_grid.cxx        
//     
//      fortran callable wrapper functions for the c++  
//      appl grid project.
//                   
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: fappl.cxx, v1.0   Wed May 21 14:31:36 CEST 2008 sutt

#include <iostream>
using std::cout;
using std::endl;

#include "appl_grid/appl_grid.h"




// externally defined alpha_s and pdf routines for fortran 
// callable convolution wrapper
extern "C" double fnalphas_(const double& Q); 
extern "C" void   fnpdf_(const double& x, const double& Q, double* f);

// create the grid
extern "C" void bookgrid_(const int& Nobs, const double* binlims);

// delete the gread
extern "C" void releasegrid_();

// read from a file
extern "C" void readgrid_(const char* s);

// write to a file 
extern "C" void writegrid_(const char* s);

// add an entry 
extern "C" void fillgrid_(const int& ix1, const int& ix2, const int& iQ, 
			  const int& iobs, 
			  const double* w,
			  const int& iorder );  

// redefine the grid dimensions
extern "C" void redefine_(const int& iobs, const int& iorder, 
			  const int& NQ2, const double& Q2min, const double& Q2max, 
			  const int& Nx,  const double&  xmin, const double&  xmax); 

// do the convolution!! hooray!!
extern "C" void convolute_(const int& nloops, double* data);

// print the grid
extern "C" void printgrid_();




appl::grid* _grid = NULL;

void bookgrid_(const int& Nobs, const double* binlims) 
{
  cout << "bookgrid_() creating grid" << endl; 

  if ( !_grid ) _grid = new appl::grid( Nobs, binlims,
					2,    10, 1000, 1,
					12,  1e-5, 1, 3, 
					"nlojet", 1, 3, "f3");
  //  _grid->symmetrise(true);
  
}


void readgrid_(const char* s) { 
  if ( !_grid ) _grid = new appl::grid(s);
}
  
  
void printgrid_() { 
  cout << *_grid << endl;
}


void releasegrid_() { 
  if ( _grid ) { delete _grid; _grid = NULL; }
}


void redefine_(const int& iobs, const int& iorder, 
	       const int& NQ2, const double& Q2min, const double& Q2max, 
	       const int& Nx,  const double&  xmin, const double&  xmax) 
{
  if ( _grid ) _grid->redefine(iobs, iorder, 
			       NQ2, Q2min, Q2max, 
			       Nx,   xmin,  xmax); 
} 


void convolute_(const int& nloops, double* data) { 
  cout << "convolute_() nloops=" << nloops << endl; 
  if ( _grid ) { 
    //    cout << "   calling _grid->convolute()" << endl;
    TH1D* h = _grid->convolute(fnpdf_, fnalphas_, nloops);
    for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) { 
      data[i-1] = h->GetBinContent(i);      
      cout << "convolute_() data[" << i-1 << "]=" << data[i-1] << endl; 
    }
  }
}



void writegrid_(const char* s) { 
  cout << "writegrid_() writing " << s << endl;
  
  if ( _grid ) {
    _grid->trim();
    //    _grid->print();
    _grid->Write(s);
  }
}



void fillgrid_(const int& ix1, const int& ix2, const int& iQ,  
	       const int& iobs, 
	       const double* w, 
	       const int& iorder ) { 
  
  //  cout << "ix " << ix1 << " " << ix2 << "  iQ" << iQ << " " << iobs << "  " << iorder << endl;  
  
  if ( _grid ) _grid->fill_index(ix1, ix2, iQ, iobs, w, iorder);
}  
  

