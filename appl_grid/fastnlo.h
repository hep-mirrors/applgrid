// emacs: this is -*- c++ -*-
//
//   @file    fastnlo.h        
//
//            native fastnlo interface for appl_grid                    
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: fastnlo.h, v0.0   Thu  8 Oct 2009 23:04:35 BST sutt $


#ifndef __FASTNLO_H
#define __FASTNLO_H

#include <iostream>
#include <vector>
#include <string>

#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_pdf.h"
#include "appl_grid/appl_timer.h"


class fastnlo {
  
public:
  
  fastnlo(const std::string& filename ); 
  
  virtual ~fastnlo() { 
    if ( m_manage_grids ) { 
      for ( int i=m_grid.size() ; i-- ; ) { delete m_grid[i]; m_grid[i]=NULL; } 
    }
  } 

  std::vector<appl::grid*> grids() { return m_grid; }
  
  void manageGrids(bool b=true) { m_manage_grids=b; }  

private:
  
  bool m_manage_grids;

  std::vector<appl::grid*> m_grid;
  
};


inline std::ostream& operator<<( std::ostream& s, const fastnlo& _f ) { 
  return s;
}


#endif  // __FASTNLO_H 










