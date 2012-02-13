//
//   TFileVector.cxx        
//
//    root TObject string vector class for writing string vectors
//    to root files               
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: TFileString.cxx, v0.0   Sat Mar 15 19:50:15 GMT 2008 sutt



#include "appl_grid/TFileVector.h"

ClassImp(TFileVector)


ostream& operator<<(ostream& s, const TFileVector& fs) { 
  return s;
}



