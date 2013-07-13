// emacs: this is -*- c++ -*-
//
//   basic_pdf.cxx        
//
//   pdf transform function for basic 121 ( 11 x 11 ) pdf combinations                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: basic_pdf.cxx, v1.0   Mon Dec 10 01:36:04 GMT 2007 sutt $



#include "basic_pdf.h"


// fortran callable wrapper
extern "C" void fbasic_pdf__(const double* fA, const double* fB, double* H) { 
  static basic_pdf pdf;
  pdf.evaluate(fA, fB, H);
}
