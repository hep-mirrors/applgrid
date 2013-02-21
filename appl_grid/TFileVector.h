// emacs: this is -*- c++ -*-
//
//   TFileVector.h        
//
//    root TObject string vector class for writing string vectors
//    to root files               
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: TFileVector.h, v0.0   Sat Mar 15 19:49:16 GMT 2008 sutt


#ifndef  TFILESTRING_H
#define  TFILESTRING_H

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "TObjString.h"
#include "TObject.h"


// template<class T>
class TFileVector : public TObjString { 

public:
  
  TFileVector(const string& name="") : TObjString(name.c_str()) { } 

  std::vector<std::vector<double> >&       histos() { return mv; }
  const std::vector<std::vector<double> >& histos() const { return mv; }

  // get the name
  string  name() const { return GetName(); } 

  // get a value 
  std::vector<double>&       operator[](int i) { return mv[i]; }
  const std::vector<double>& operator[](int i) const { return mv[i]; }
  
  // get the size
  unsigned size()  const { return mv.size(); } 

  // add an element
  void add(std::vector<double>& s) { 
    mv.push_back(s); 
  }

private:
  
  std::vector<std::vector<double> > mv;

  ClassDef(TFileVector, 1)

}; 


ostream& operator<<(ostream& s, const TFileVector& fv);



#endif  //  TFILEVECTOR_H 










