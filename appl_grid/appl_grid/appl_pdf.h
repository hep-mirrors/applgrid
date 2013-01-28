
// emacs: this is -*- c++ -*-
//
//   appl_pdf.h        
//
//   pdf transform functions header                  
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: appl_pdf.h, v   Fri Dec 21 22:19:50 GMT 2007 sutt $


#ifndef __APPL_PDF_H
#define __APPL_PDF_H

#include <iostream>
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

// #include <sstream>
// using std::ostringstream;
// using std::stringstream;

#include <vector> 
using std::vector;

#include <map> 
using std::map;

#include <string> 
using std::string;


namespace appl { 


class appl_pdf;
typedef map <const string, appl_pdf*> pdfmap;


// this is a *maybe* nice class, a base class for pdf 
// functions
//
// it has a virtual evaluate() method to be definied in 
// the derived class, and a static map of all the names
// of instances of the derived classes
//
// when a new instance of the class is created, it 
// automatically adds it's name to the map, so the user 
// doesn't need to worry about consistency, and removes 
// itself when the derived instance is deleted

class appl_pdf { 

public:
  appl_pdf(string name) : m_Nproc(0), m_name(name) { 
    addtopdfmap(m_name, this);
  }
  
  virtual ~appl_pdf() { 
    // when I'm destroyed, remove my entry from the map 
    m_pdfmap.erase(m_pdfmap.find(m_name));
  } 

  virtual void evaluate(const double* f1, const double* f2, double* H) = 0; 

  int     Nproc() const { return m_Nproc; } 
  string   name() const { return m_name;  }

  string  rename(const std::string& name) { 
    /// remove my entry from the map, and add me again with my new name
    m_pdfmap.erase(m_pdfmap.find(m_name));
    m_name = name;
    addtopdfmap(m_name, this);
    return m_name;
  }

  // retrieve an instance from the map 
  static appl_pdf* getpdf(const string& s) { 
    if ( m_pdfmap.find(s)!=m_pdfmap.end() ) return m_pdfmap.find(s)->second;
    else  throw exception( cerr << "getpdf() " << s << " not instantiated in map" );
  }

  // pdf error exception
  class exception { 
  public: 
    exception(const string& s="") { cerr << s << endl; }; 
    exception(ostream& s)         { cerr << s << endl; }; 
    //    exception(ostringstream& s)   { cerr << s << endl; }; 
    //    exception(stringstream& s)    { cerr << s << endl; }; 
  };
  


  /// code to allow optional vector of subprocess contribution names

  const vector<string>& subnames() const { return m_subnames; }

  void addSubnames( const vector<string>& subnames ) { m_subnames = subnames; }

  void  addSubname( const string& subname ) { 
    if ( int(m_subnames.size())<m_Nproc-1 ) m_subnames.push_back(subname); 
  }


private:

  static void addtopdfmap(const string& s, appl_pdf* f) { 
    if ( m_pdfmap.find(s)==m_pdfmap.end() ) { 
      m_pdfmap[s] = f;
      //      cout << "appl_pdf::addtomap() registering " << s << " in map addr \t" << f << endl;
    }
    else { 
      throw exception( cerr << "appl_pdf::addtopdfmap() " << s << " already in map\t0x" << m_pdfmap.find(s)->second  );
    }
  }

protected:

  int    m_Nproc;
  string m_name;

  vector<string> m_subnames;
  
  static pdfmap m_pdfmap;

};



};


#endif  // __APPL_PDF_H 










