// emacs: this is -*- c++ -*-
//
//   @file    mcfm_grid.h        
//
//                   
//  
//   Copyright (C) 2011 M.Sutton (sutt@cern.ch)    
//
//   $Id: mcfm_grid.h, v0.0   Fri 12 Aug 2011 07:50:44 CEST sutt $


#ifndef __MCFM_GRID_H
#define __MCFM_GRID_H

#include <iostream>

#include "appl_grid/appl_grid.h"

namespace appl { 

class mcfm_grid : public  grid {

public:

  //  mcfm_grid() { } 

  mcfm_grid(int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0,  int Q2order=5,  
	    int Nx=50,   double xmin=1e-5,     double xmax=0.9,          int xorder=5,
	    int Nobs=20, double obsmin=100.0,  double obsmax=7000.0, 
	    string genpdf="mcfm_pdf", 
	    int leading_order=0, int nloops=1, 
	    string transform="f2")  
    : grid( NQ2,  Q2min,  Q2max,  Q2order, Nx, xmin, xmax, xorder,
	    Nobs, obsmin, obsmax, genpdf,  leading_order,  nloops, 
	    transform )  { } 

  mcfm_grid( int Nobs, const double* obsbins,
	     int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0, int Q2order=5,
	     int Nx=50,   double xmin=1e-5,     double xmax=0.9,         int xorder=5, 
	     string genpdf="mcfm_pdf",
	     int leading_order=0, int nloops=1, 
	     string transform="f2" )
    : grid( Nobs, obsbins, NQ2,  Q2min,  Q2max,  Q2order, Nx, xmin, xmax, xorder,
	    genpdf,  leading_order,  nloops, 
	    transform )  { } 
  
  mcfm_grid( const vector<double> obs,
	     int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0,   int Q2order=5, 
	     int Nx=50,   double xmin=1e-5,     double xmax=0.9,           int xorder=5, 
	     string genpdf="mcfm_pdf", 
	     int leading_order=0, int nloops=1, 
	     string transform="f2" )
    : grid( obs, NQ2,  Q2min,  Q2max,  Q2order, Nx, xmin, xmax, xorder,
	    genpdf,  leading_order,  nloops, 
	    transform )  { } 
  
  // build a grid but don't build the internal igrids - these can be added later
  mcfm_grid( const vector<double> obs,
	     string genpdf="nlojet_pdf", 
	     int leading_order=0, int nloops=1, 
	     string transform="f2" )
    : grid( obs,  genpdf,  leading_order,  nloops, transform )  { } 
  
  
  mcfm_grid(const grid& g) : grid(g) { } 
  
  mcfm_grid(const string& filename="./grid.root", const string& dirname="grid") : grid(filename, dirname) { } 
  

  // destructor
  virtual ~mcfm_grid() { } 

  
  //  fill weight from MCFM common block  
  void fillMCFM(double obs);
  void collectWeight(const int& order, const int& id, double* wt);
  void decideSubProcess(const int& iflav1, const int& iflav2, int & iProcess , double &factor );

};


};



#endif  // __MCFM_GRID_H 










