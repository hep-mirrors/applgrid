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


#include "TFile.h"
#include "TVectorT.h"
#include "TMatrixT.h"



// get ckm related information

void mcfmw_pdf::make_ckmsum() { 
  // cout << "make_ckmsum() initialising" << endl;
  m_ckmsum = new double[13];
  
  std::string ckmfile = std::string(_DIR_)+"/ckm.root";
  TFile f(ckmfile.c_str());
  std::cout << "ckmfile " << ckmfile << std::endl; 
  TVectorT<double>* _ckmsum = (TVectorT<double>*)f.Get("ckmsum");
  for ( int i=0 ; i<13 ; i++ ) { 
    m_ckmsum[i] = (*_ckmsum)(i);
    // _ckm[i] = 1;
  }
  f.Close();
}


void mcfmw_pdf:: make_ckm() {  
  // cout << "make_ckm() initialising" << endl;
  m_ckm2 = new double*[13];
  std::string ckmfile = std::string(_DIR_)+"/ckm.root";
  TFile f(ckmfile.c_str());
  TMatrixT<double>* _ckm2 = (TMatrixT<double>*)f.Get("ckm2");
  for ( int i=0 ; i<13 ; i++ ) { 
    m_ckm2[i] = new double[13];
    for ( int j=0 ; j<13 ; j++ )  m_ckm2[i][j] = (*_ckm2)(i,j);
    //    _ckm[i][i] = 1;
  }
  f.Close();
}


// fortran callable wrapper

extern "C" void fmcfmw_pdf__(const double* fA, const double* fB, double* H) { 
  static mcfmw_pdf pdf;
  pdf.evaluate(fA, fB, H);
}




