// emacs: this is -*- c++ -*-
//
//   mcfmw_pdf.cxx        
//
//   pdf transform functions                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: mcfmw_pdf.cxx, v1.0   Mon Dec 10 01:36:04 GMT 2007 sutt $



#include "appl_grid/appl_pdf.h" 
using namespace appl;

#include "appl_grid/mcfmw_pdf.h"

// get ckm related information

void mcfmw_pdf::make_ckmsum() { 
  // cout << "make_ckmsum() initialising" << endl;
  m_ckmsum = new double[13];

  double ckmsum[13] = { 
    0.000000000000000000000000000000000000 , 
    0.000000000000000000000000000000000000 , 
    0.000000000000000000000000000000000000 , 
    0.999908999999999936747485662635881454 , 
    0.000000000000000000000000000000000000 , 
    0.999908999999999936747485662635881454 , 
    0.000000000000000000000000000000000000 , 
    0.000000000000000000000000000000000000 , 
    0.999908999999999936747485662635881454 , 
    0.000000000000000000000000000000000000 , 
    0.999908999999999936747485662635881454 , 
    0.000000000000000000000000000000000000 , 
    0.000000000000000000000000000000000000 
  }; 
    
  for ( int i=0 ; i<13 ; i++ ) m_ckmsum[i] = ckmsum[i];

}


void mcfmw_pdf:: make_ckm() {  
  // cout << "make_ckm() initialising" << endl;
  m_ckm2 = new double*[13];

  for ( int i=0 ; i<13 ; i++ ) {
    m_ckm2[i] = new double[13]; 
    for ( int j=0 ; j<13 ; j++ ) m_ckm2[i][j] = 0;
  }

  m_ckm2[3][8]  =   0.049284000000000001417976847051249933 ;
  m_ckm2[8][3]  =   0.049284000000000001417976847051249933 ;

  m_ckm2[5][8]  =   0.950624999999999942268402719491859898 ;
  m_ckm2[8][5]  =   0.950624999999999942268402719491859898 ;

  m_ckm2[5][10] =   0.049284000000000001417976847051249933 ;
  m_ckm2[10][5] =   0.049284000000000001417976847051249933 ;

  m_ckm2[3][10] =   0.950624999999999942268402719491859898 ;
  m_ckm2[10][3] =   0.950624999999999942268402719491859898 ;

}


// fortran callable wrapper

extern "C" void fmcfmw_pdf__(const double* fA, const double* fB, double* H) { 
  static mcfmw_pdf pdf;
  pdf.evaluate(fA, fB, H);
}




