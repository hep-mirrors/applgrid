
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


#include <exception> 


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

  // pdf error exception
  class exception : public std::exception { 
  public: 
    exception(const string& s="") { cerr << what() << " " << s << endl; }; 
    exception(ostream& s)         { cerr << what() << " " << s << endl; }; 
    virtual const char* what() const throw() {  return "appl::appl_pdf::exception "; }
  };
  
public:

  // retrieve an instance from the map 
  static appl_pdf* getpdf(const string& s) {
    pdfmap::iterator itr = m_pdfmap.find(s);
    if ( itr!=m_pdfmap.end() ) return itr->second;    
    throw exception( cerr << "getpdf() " << s << " not instantiated in map " );
  }

  
  static void printmap() {
    pdfmap::iterator itr = m_pdfmap.begin();
    while ( itr!=m_pdfmap.end() )  {
      std::cout << "pdfmap " << itr->first << "\t\t" << itr->second << std::endl;
      itr++;
    } 
  }
  
public:

  appl_pdf(const string& name) : m_Nproc(0), m_name(name) { 
    if ( m_name!="" ) addtopdfmap(m_name, this);
  }
  
  virtual ~appl_pdf() { 
    // when I'm destroyed, remove my entry from the map 
    pdfmap::iterator mit = m_pdfmap.find(m_name);
    if ( mit!=m_pdfmap.end() ) m_pdfmap.erase(mit);
  } 

  virtual void evaluate(const double* f1, const double* f2, double* H) = 0; 

  int     Nproc() const { return m_Nproc; } 
  string   name() const { return m_name;  }

  string  rename(const std::string& name) { 
    /// remove my entry from the map, and add me again with my new name
    if ( m_pdfmap.find(m_name)!=m_pdfmap.end() ) { 
      m_pdfmap.erase(m_pdfmap.find(m_name));
    }
    else { 
      std::cout << "appl_pdf::rename() " << m_name << " not in map" << std::endl;
    }
    m_name = name;
    addtopdfmap(m_name, this);
    return m_name;
  }


  /// code to allow optional vector of subprocess contribution names

  const vector<string>& subnames() const { return m_subnames; }

  void addSubnames( const vector<string>& subnames ) { m_subnames = subnames; }

  void  addSubname( const string& subname ) { 
    if ( int(m_subnames.size())<m_Nproc-1 ) m_subnames.push_back(subname); 
  }

  /// access the ckm matrices - if no matrices are required these vectors have 
  /// zero size

  const std::vector<double>&               getckmsum() const { return m_ckmsum; }
  const std::vector<std::vector<double> >& getckm2()   const { return m_ckm2; }

  /// set the ckm matrices from external values
  void setckm2( const std::vector<std::vector<double> >& ckm2 ); 

  /// code to create the ckm matrices using the hardcoded default 
  /// values if required
  /// takes a bool input - if true creates the ckm for Wplus, 
  /// false for Wminus
  void make_ckm( bool Wp=true );

  void SetNProc(int Nsub){ m_Nproc=Nsub; return;};
  
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

  // ckm matrix related information 
  std::vector<double>               m_ckmsum;
  std::vector<std::vector<double> > m_ckm2;

};



};


#endif  // __APPL_PDF_H 










