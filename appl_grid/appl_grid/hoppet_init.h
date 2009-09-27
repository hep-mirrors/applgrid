// emacs: this is -*- c++ -*-
//
//   @file    hoppet_init.h        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: hoppet_init.h, v0.0   Sat 26 Sep 2009 13:24:26 BST sutt $


#ifndef __HOPPETINIT_H
#define __HOPPETINIT_H

#include <iostream>
#include <vector>
#ifdef HOPPET
#include "hoppet_v1.h"
#endif

class hoppet_init {

public:

  hoppet_init();

  virtual ~hoppet_init();

  void fillCache( void (*pdf)(const double&, const double&, double* )  );
  
  bool compareCache( void (*pdf)(const double&, const double&, double* )  );

private:

  std::vector<double> m_cache;

};



inline std::ostream& operator<<( std::ostream& s, const hoppet_init& _h ) { 
  return s;
}


#endif  // __HOPPETINIT_H 










