//
//   TFileString.cxx        
//
//    root TObject string vector class for writing string vectors
//    to root files               
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: TFileString.cxx, v0.0   Sat Mar 15 19:50:15 GMT 2008 sutt



#include "appl_grid/TFileString.h"

ClassImp(TFileString)


ostream& operator<<(ostream& s, const TFileString& fs) { 
  s << fs.name() << ":";
  for ( unsigned i=0 ; i<fs.size() ; i++ ) s << "\t" << fs[i];
  return s;
}



