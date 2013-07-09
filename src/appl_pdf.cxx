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

#include "mcfmz_pdf.h"
#include "mcfmw_pdf.h"
#include "mcfmwc_pdf.h"
#include "mcfmQQ_pdf.h"
#include "nlojet_pdf.h"
#include "jetrad_pdf.h"
#include "nlojetpp_pdf.h"
#include "dis_pdf.h"
#include "vrapz_pdf.h"
// #include "generic_pdf.h"


namespace appl { 


// initialise the map with some default instances
// although the user could create these themselves
// if they wanted
pdfmap appl_pdf::m_pdfmap; 


/// constructor and destructor
appl_pdf::appl_pdf(const string& name) : m_Nproc(0), m_name(name) { 
   if ( m_name!="" ) addtopdfmap(m_name, this);
}
  
appl_pdf::~appl_pdf() { 
  // when I'm destroyed, remove my entry from the map 
  pdfmap::iterator mit = m_pdfmap.find(m_name);
  if ( mit!=m_pdfmap.end() ) m_pdfmap.erase(mit);
} 


/// retrieve an instance from the map 
appl_pdf* appl_pdf::getpdf(const string& s, bool printout) {
  /// initialise the factory
  if ( m_pdfmap.size()==0 ) appl::appl_pdf::create_map(); 
  pdfmap::iterator itr = m_pdfmap.find(s);
  if ( itr!=m_pdfmap.end() ) return itr->second; 
  /// not found in map
  throw exception( std::cerr << "getpdf() " << s << " not instantiated in map " );
}


bool appl_pdf::create_map() { 

  std::cout << "appl_pdf::create_map() creating pdf combination factory" << std::endl;

  if ( m_pdfmap.size()==0 ) { 
    
    /// the appl_pdf add their own pointers to the 
    /// pdf map so we don;t need to remember their 
    /// pointers ourselves
    new  mcfmz_pdf;
    new  mcfmwp_pdf;
    new  mcfmwm_pdf;
    
    new  mcfmwpc_pdf;
    new  mcfmwmc_pdf;
    
    new  mcfmCC_pdf;
    new  mcfmBB_pdf;
    new  mcfmTT_pdf;

    new nlojet_pdf;
    new nlojetpp_pdf;
    new jetrad_pdf;
    new dis_pdf;

    new vrapzLO_pdf;
    new vrapzNLO_pdf;
    new vrapzNNLO_pdf;

    //    printmap( std::cerr );

  }

  return true;
}




/// get ckm related information 


void appl_pdf::make_ckm( bool Wp ) { 
  
  // cout << "make_ckm() initialising" << endl;
  m_ckm2 = std::vector<std::vector<double> >(13, std::vector<double>(13,0));

  if ( Wp ) { 
    
    //  std::cout << "creating ckm matrix terms for Wplus production" << std::endl;
    
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
    
    //    std::cout << "creating ckm matrix terms for Wminus production" << std::endl;
    
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

  /*
  for ( int i=0 ; i<13 ; i++ ) {
    for ( int j=0 ; j<13 ; j++ ) {
      if (m_ckm2[i][j]!=0)
	std::cout << " ckm[" << i << "][" << j << "]\t =\t " << m_ckm2[i][j] << std::endl;
    }
  }
  */
}


void appl_pdf::setckm2( const std::vector<std::vector<double> >& ckm2 ) { 
    m_ckm2 = ckm2; 
    m_ckmsum = std::vector<double>(m_ckm2.size(),0);
    for ( unsigned i=0 ; i<m_ckm2.size() ; i++ ) { 
      for ( unsigned j=0 ; j<m_ckm2[i].size() ; j++ ) m_ckmsum[i] += m_ckm2[i][j]; 
    }
} 


};




