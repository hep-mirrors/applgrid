// emacs: this is -*- c++ -*-
//
//   appl_pdf.cxx        
//
//   pdf transform functions                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: appl_pdf.cxx, v1.0   Mon Dec 10 01:36:04 GMT 2007 sutt $


#include "appl_grid/appl_pdf.h" 
using namespace appl;

#include "appl_grid/mcfmz_pdf.h"
#include "appl_grid/mcfmw_pdf.h"
#include "appl_grid/mcfmQQ_pdf.h"
#include "appl_grid/nlojet_pdf.h"
#include "appl_grid/jetrad_pdf.h"
#include "appl_grid/nlojetpp_pdf.h"
#include "appl_grid/dis_pdf.h"
#include "appl_grid/vrapz_pdf.h"
// #include "appl_grid/generic_pdf.h"


namespace appl { 

// initialise the map with some default instances
// although the user could create these themselves
// if they wanted
pdfmap appl_pdf::m_pdfmap; 

mcfmz_pdf     mcfmzpdf;
mcfmwp_pdf    mcfmwppdf;
mcfmwm_pdf    mcfmwmpdf;

mcfmCC_pdf    mcfmccpdf;
mcfmBB_pdf    mcfmbbpdf;
mcfmTT_pdf    mcfmttpdf;

nlojet_pdf    nlojetpdf;
nlojetpp_pdf  nlojetpppdf;
jetrad_pdf    jetradpdf;
dis_pdf       dispdf;
vrapzLO_pdf   vrapzLOpdf;
vrapzNLO_pdf  vrapzNLOpdf;
vrapzNNLO_pdf vrapzNNLOpdf;

// generic_pdf   genericpdf;



/// get ckm related information 


void appl_pdf::make_ckm( bool Wp ) { 
  
  // cout << "make_ckm() initialising" << endl;
  m_ckm2 = std::vector<std::vector<double> >(13, std::vector<double>(13,0));

  if ( Wp ) { 
    
    std::cout << "creating ckm matrix terms for Wplus production" << std::endl;
    
    m_ckm2[3][8]  =   0.049284000000000001417976847051249933 ;
    m_ckm2[8][3]  =   0.049284000000000001417976847051249933 ;
    
    m_ckm2[5][8]  =   0.950624999999999942268402719491859898 ;
    m_ckm2[8][5]  =   0.950624999999999942268402719491859898 ;
    
    m_ckm2[5][10] =   0.049284000000000001417976847051249933 ;
    m_ckm2[10][5] =   0.049284000000000001417976847051249933 ;
    
    m_ckm2[3][10] =   0.950624999999999942268402719491859898 ;
    m_ckm2[10][3] =   0.950624999999999942268402719491859898 ;
    
  }
  else { 
    
    std::cout << "creating ckm matrix terms for Wminus production" << std::endl;
    
    m_ckm2[4][9] =   0.049284000000000001417976847051249933 ;
    m_ckm2[9][4] =   0.049284000000000001417976847051249933 ;
    
    m_ckm2[7][4] =   0.950624999999999942268402719491859898 ;
    m_ckm2[4][7] =   0.950624999999999942268402719491859898 ;
    
    m_ckm2[7][2] =   0.049284000000000001417976847051249933 ;
    m_ckm2[2][7] =   0.049284000000000001417976847051249933 ;
    
    m_ckm2[9][2] =   0.950624999999999942268402719491859898 ;
    m_ckm2[2][9] =   0.950624999999999942268402719491859898 ;
    
  }  
 

  m_ckmsum = std::vector<double>(m_ckm2.size(),0);
  for ( unsigned i=0 ; i<m_ckm2.size() ; i++ ) { 
    for ( unsigned j=0 ; j<m_ckm2[i].size() ; j++ ) m_ckmsum[i] += m_ckm2[i][j]; 
  }

}


void appl_pdf::setckm2( const std::vector<std::vector<double> >& ckm2 ) { 
    m_ckm2 = ckm2; 
    m_ckmsum = std::vector<double>(m_ckm2.size(),0);
    for ( unsigned i=0 ; i<m_ckm2.size() ; i++ ) { 
      for ( unsigned j=0 ; j<m_ckm2[i].size() ; j++ ) m_ckmsum[i] += m_ckm2[i][j]; 
    }
} 


};




