// emacs: this is -*- c++ -*-
//
//   basic_pdf.h        
//
//   pdf transform function for basic 121 ( 11 x 11 ) combinations                    
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: basic_pdf.h, v1.0   Mon Dec 10 01:36:04 GMT 2007 sutt $



#ifndef BASIC_PDF_H
#define BASIC_PDF_H


#include "appl_grid/appl_pdf.h" 


// basic pdf combination class
class basic_pdf : public appl::appl_pdf {

public:

  basic_pdf() : appl::appl_pdf("basic") { m_Nproc=121; } 

  void evaluate(const double* f1, const double* f2, double* H);

};  
  


inline void basic_pdf::evaluate(const double* _f1, const double* _f2, double* H) {  

  // remapping from pdg -6..6 convention to u..t ubar..tbar g internal
  // basic convention

  const double* f1 = _f1+6;
  const double* f2 = _f1+6;

  ///  double H[121]; /// do not include top: 11x11 rather than 13x13
  int ih=0;
  /// from -b to b rather than -t to t
  for ( int i=-5 ; i<=5 ; i++ )  { 
    for ( int j=-5 ; j<=5 ; j++ )  { 
      /// f1 partons first !!! 
      H[ih++] = f1[i]*f2[j];
    }
  }
}
  


// fortran callable wrapper
extern "C" void fbasic_pdf__(const double* fA, const double* fB, double* H);


#endif // BASIC_PDF_H

