// emacs: this is -*- c++ -*-
//
//   TFileString.h        
//
//    root TObject string vector class for writing string vectors
//    to root files               
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: TFileString.h, v0.0   Sat Mar 15 19:49:16 GMT 2008 sutt


#ifndef __TFILESTRING_H
#define __TFILESTRING_H

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


class TFileString : public TObjString { 

public:
  
  TFileString(const string& name="") : TObjString(name.c_str()) { } 

  TFileString(const string& name, const string& tag) : 
    TObjString(name.c_str())
  { mstring.push_back(tag.c_str()); } 
 
  vector<string>&       tags()       { return mstring; }
  const vector<string>& tags() const { return mstring; }

  // get the name
  string  name() const { return GetName(); } 

  // get a value 
  string& operator[](int i)       { return mstring[i]; }
  string  operator[](int i) const { return mstring[i]; }
  
  // get the size
  unsigned size()           const { return mstring.size(); } 

  // add an element
  void add(const string& s) { mstring.push_back(s); }

private:
  
  vector<string> mstring;

  ClassDef(TFileString, 1)

}; 


ostream& operator<<(ostream& s, const TFileString& fs);



#endif  // __TFILESTRING_H 










