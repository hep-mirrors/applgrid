//
//   @file    lumi_pdf.cxx         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2013 M.Sutton (sutt@cern.ch)    
//
//   $Id: lumi_pdf.cxx, v0.0   Tue  9 Jul 2013 08:14:47 CEST sutt $


#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>


#include "appl_grid/lumi_pdf.h"



lumi_pdf::lumi_pdf(const std::string& s, const std::vector<int>& combinations ) :  //, bool amcflag ) : 
  appl_pdf(s), m_filename(s) //,  m_amcflag(amcflag)
{

  std::cout << "lumi_pdf::lumi_pdf() " << s << "\tv size " << combinations.size() << std::endl; 

  /// need to decode the input std::vector

  if ( combinations.size() ) { 
    /// std::vector initialised from serialised std::vector
    unsigned iv = 0;
    unsigned nproc = combinations[iv++];
    for ( unsigned i=0 ; i<nproc && iv<combinations.size() ; i++ ) {
      int index  = combinations[iv++];
      int npairs = combinations[iv++];
      std::vector<int> v(npairs*2+2);
      v[0] = index;
      v[1] = npairs;
      
      for ( int j=0 ; j<npairs ; j++ ) { 
	v[j*2+2] = combinations[iv++];
	v[j*2+3] = combinations[iv++];
      }

      combination c(v);
      if ( c.size() ) add(c); 
    } 
  }
  else { 
    /// else read from file ...
    
    std::ifstream infile(m_filename.c_str());

    std::string   line;

    while (std::getline(infile, line)) {
      //    std::cout << "line: " << line << std::endl;
      combination c( line );
      if ( c.size() ) add(c); 
    }
    
  }
 
  // some checking

  //  for ( int i=0 ; i<m_combinations.size() ; i++ ) { 
  //    if ( m_combinations[i].
  //  }

  m_Nproc = m_combinations.size();

  std::cout << *this << std::endl;

}



void lumi_pdf::evaluate(const double* xfA, const double* xfB, double* H) { 
  for ( unsigned i=size() ; i-- ; ) { 
    H[i] = m_combinations[i].evaluate( xfA, xfB ); 
  }
}




int  lumi_pdf::decideSubProcess(const int iflav1, const int iflav2) { 
  /// do we really need this? in any case, need to sort though each subprocess
  /// checking if the provided flavours are in any of the pairs
  /// so quicker to set up a reverse 13x13 lookup when needed, and not 
  /// otherwise
  if ( m_lookup.size()==0 ) { 
    /// create a 13 x 13 lookup table 
    m_lookup = std::vector<std::vector<int> >(13, std::vector<int>(13, -1) ); 
    for ( unsigned i=size() ; i-- ; ) { 
      const combination& c = m_combinations[i];
      for ( unsigned j=c.size() ; j-- ; ) m_lookup[ c[j].first+6 ][ c[j].second+6 ] = i;
    } 
  }

  return m_lookup[iflav1+6][iflav2+6];
}



std::vector<int> lumi_pdf::serialise() const { 

  std::vector<int> v;

  v.push_back( Nproc() );
  for ( int i=0 ; i<Nproc() ; i++ ) { 
    v.push_back(i);
    const combination& c = m_combinations[i];
    v.push_back( c.size() );
    for ( unsigned j=0 ; j<c.size() ; j++ ) { 
      v.push_back( c[j].first );
      v.push_back( c[j].second );
    }
  }  
  
  return v;
}


