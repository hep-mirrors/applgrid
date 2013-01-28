// emacs: this is -*- c++ -*-
//
//   @file    generic_pdf.h        
//
//            a generic pdf type - to read in the combinations
//            for the subprocesses from a file, with the file 
//            format determined by Tancredi (details will be 
//            filled in as the implementation becomes more complete)                       
//  
//   Copyright (C) 2013 M.Sutton (sutt@cern.ch)    
//
//   $Id: generic_pdf.h, v0.0   Mon 28 Jan 2013 15:41:10 GMT sutt $


#ifndef  GENERIC_PDF_H
#define  GENERIC_PDF_H

#include <iostream>

#include "appl_grid/appl_pdf.h" 
using namespace appl;




class generic_pdf : public appl_pdf {

public:

  generic_pdf() : appl_pdf("generic"), m_initialised(false) {
    /// set up later generic pdf here ...
  } 

  virtual ~generic_pdf() { } 

  void evaluate(const double* fA, const double* fB, double* H);

  /// additional user defined functions to actually initialise 
  /// based on the input file

  void initialise(const std::string& file);

  bool initialised() const { return m_initialised; }

private:

  string m_filename;  /// this might eventually become a string encoding the grid

  bool   m_initialised;

};



inline void  generic_pdf::evaluate(const double* fA, const double* fB, double* H) {  
  ///  fill this in with tancredi's code ...
  if ( !m_initialised ) return;

}



inline std::ostream& operator<<( std::ostream& s, const generic_pdf& _g ) { 
  return s;
}


#endif  // GENERIC_PDF_H 










