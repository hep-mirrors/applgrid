// emacs: this is -*- c++ -*-
//
//   @file    fappl_grid.h        
//
//                   
//  
//   Copyright (C) 2013 M.Sutton (sutt@cern.ch)    
//
//   $Id: fappl_grid.h, v0.0   Fri 25 Oct 2013 00:52:32 CEST sutt $


#pragma once


/// externally defined alpha_s and pdf routines for fortran 
/// callable convolution wrapper
extern "C" double fnalphas_(const double& Q);


extern "C" void   fnpdf_(const double& x, const double& Q, double* xf);

/// create a grid
extern "C" void bookgrid_(int& id, const int& Nobs, const double* binlims);

/// delete a grid
extern "C" void releasegrid_(int& id);

/// delete all grids
extern "C" void releasegrids_();

/// read a grid from a file
extern "C" void readgrid_(int& id, const char* s);

/// write to a file 
extern "C" void writegrid_(int& id, const char* s);

/// add an entry 
extern "C" void fillgrid_(int& id, 
			  const int& ix1, const int& ix2, const int& iQ, 
			  const int& iobs, 
			  const double* w,
			  const int& iorder );  

/// redefine the grid dimensions
extern "C" void redefine_(int& id, 
			  const int& iobs, const int& iorder, 
			  const int& NQ2, const double& Q2min, const double& Q2max, 
			  const int& Nx,  const double&  xmin, const double&  xmax); 


/// get number of observable bins for a grid 
extern "C" int getnbins_(int& id);

/// get the number of a bin based on the ordonate 
extern "C" int getbinnumber_(int& id, double& x);

/// get low edge for a grid/bin number
extern "C" double getbinlowedge_(int& id, int& bin);

/// get width for a grid/bin number
extern "C" double getbinwidth_(int& id, int& bin);

/// do the convolution!! hooray!!

extern "C" void convolute_(int& id, double* data);

extern "C" void convoluteorder_(int& id, int& nloops, double* data);

extern "C" void convolutewrap_(int& id, double* data, 
			       void (*pdf)(const double& , const double&, double* ),
			       double (*alphas)(const double& ) );

extern "C" void fullconvolutewrap_(int& id, double* data, 
				   void (*pdf)(const double& , const double&, double* ),
				   double (*alphas)(const double& ),
				   int& nloops,
				   double& rscale, double& fscale  );



extern "C" void escaleconvolute_(int& id, double* data, double& Escale);

extern "C" void escaleconvolutewrap_(int& id, double* data, 
				     void (*pdf)(const double& , const double&, double* ),
				     double (*alphas)(const double& ), double& Escale );

extern "C" void escalefullconvolutewrap_(int& id, double* data, 
					 void (*pdf)(const double& , const double&, double* ),
					 double (*alphas)(const double& ),
					 int& nloops,
					 double& rscale, double& fscale,
					 double& Escale );


/// set the ckm matrix - a flat vector of 9 doubles, Vud, Vus, Vub, Vcd ...
extern "C" void setckm_( int& id, const double* ckm );

/// get the ckm matrix - a flat vector of 9 doubles, Vud, Vus, Vub, Vcd ...
extern "C" void getckm_( int& id, double* ckm );


/// print a grid
extern "C" void printgrid_(int& id);

/// print all grids
extern "C" void printgrids_();

/// print the grid documentation
extern "C" void printgriddoc_(int& id);

/// create grids from fastnlo
extern "C" void readfastnlogrids_(  int* ids, const char* s );


/// grid std::map management

extern "C" void ngrids_(int& n);

extern "C" void gridids_(int* ids);










