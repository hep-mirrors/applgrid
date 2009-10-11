//
//   @file    hoppet_init.cxx         
//   
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: hoppet_init.cxx, v0.0   Sat 26 Sep 2009 13:24:31 BST sutt $

#include "appl_grid/hoppet_init.h"

#include <iostream>
#include <vector>
#include <cmath>


hoppet_init::hoppet_init() {

  double dy = 0.1;   
  int nloop = 2;         

  std::cout << "hoppet_init::hoppet_init()  dy = " << dy << "\tnloop = " << nloop << std::endl; 
  // the internal grid spacing (smaller->higher accuarcy)
  // 0.1 should provide at least 10^{-3} accuracy
  // the number of loops to initialise (max=3!)  

#ifdef HOPPET
  hoppetstart_( dy, nloop );
#else
  std::cerr << " hoppet_init::hoppet_init() called but hoppet not included" << std::endl;
  exit(0);
#endif
}


hoppet_init::~hoppet_init() {
  //  std::cout << "hoppet_init::~hoppet_init()" << std::endl;
}


void hoppet_init::fillCache( void (*pdf)(const double&, const double&, double* )  ) {

  //  fill the cache

  //  hoppetassign_( pdf );

  m_cache.clear();

  for ( double lQ=2 ; lQ<=6 ; lQ+=2 ) { 
    double Q = std::pow(10,lQ);
    for ( double lx=-5 ; lx<0 ; lx+=1 ) { 
      double x = std::pow(10,lx);
      double xf[13];
      pdf( x, Q, xf ); 
      for ( int i=0 ; i<13 ; i++ ) m_cache.push_back( xf[i] ); 

      // for ( int i=0 ; i<13 ; i++ ) std::cout << "\t" << xf[i];
      // std::cout << std::endl;

      //      if ( lQ==2 || lQ==5 ) { 
      //	std::cout << Q << "\t" << x << "\t\t";
      //	for ( int i=0 ; i<13 ; i++ ) std::cout << xf[i] << "\t";
      //	std::cout << std::endl;
      //      }

    }
  }  

}



bool hoppet_init::compareCache( void (*pdf)(const double&, const double&, double* )  ) {
 
  // fill a seperate cache for comparison

  // flag whether the cache has changed
  bool changed = false;

  // is the cahche empty?
  if ( m_cache.size()==0 ) {
#     ifdef HOPPET 
      hoppetassign_( pdf );
#     endif
      fillCache( pdf );
      return changed = true;  
  }

  // fill the new cache for comparison
  //  std::cout << "ccache.size() = " << ccache.size() << "\tm_cache.size() = " << m_cache.size() << std::endl;  

  std::vector<double> ccache;
    
  for ( double lQ=2 ; lQ<=6 ; lQ+=2 ) { 
    double Q = std::pow(10,lQ);
    for ( double lx=-5 ; lx<0 ; lx+=1 ) { 
      double x = std::pow(10,lx);
      double xf[13];
      pdf( x, Q, xf ); 
      for ( int i=0 ; i<13 ; i++ ) ccache.push_back( xf[i] ); 
      
      //      if ( lQ==2 || lQ==5 ) { 
      //	std::cout << Q << "\t" << x << "\t\t";
      //	for ( int i=0 ; i<13 ; i++ ) std::cout << xf[i] << "\t";
      //	std::cout << std::endl;
      //      }

    }
  }
  

  // now compare with existing cache
  
  //  std::cout << "ccache.size() = " << ccache.size() << "\tm_cache.size() = " << m_cache.size() << std::endl;  
  
  if ( ccache.size()!= m_cache.size() ) changed = true;  
  for ( unsigned i=0 ; i<ccache.size() ; i++ ) if ( ccache[i]!=m_cache[i] ) changed = true;
  
  // cache (and hence the pdf) has changed, so replace the old cached 
  // values with the new and reinitialise hoppet with this new pdf 
  if ( changed ) { 
#   ifdef HOPPET 
    hoppetassign_( pdf );
#   endif
    m_cache = ccache;
  }
  
  //  std::cout << "hoppet_init::compareCache() changed = " << changed << std::endl; 

  return changed;
  
}

