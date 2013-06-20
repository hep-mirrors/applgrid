// emacs: this is -*- c++ -*-
//
//   nlojet_pdf.cxx        
//
//   pdf transform function for nlojet                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: nlojet_pdf.cxx, v1.0   Mon Dec 10 01:36:04 GMT 2007 sutt $



#include "appl_grid/appl_pdf.h" 
using namespace appl;

#include "appl_grid/vrapz_pdf.h"


extern "C" void fvrapzLO_pdf__(const double* fA, const double* fB, double* H) { 
  static vrapzLO_pdf pdf;
  pdf.evaluate(fA, fB, H);
}


extern "C" void fvrapzNLO_pdf__(const double* fA, const double* fB, double* H) { 
  static vrapzNLO_pdf pdf;
  pdf.evaluate(fA, fB, H);
}


extern "C" void fvrapzNNLO_pdf__(const double* fA, const double* fB, double* H) { 
  static vrapzNNLO_pdf pdf;
  pdf.evaluate(fA, fB, H);
}




