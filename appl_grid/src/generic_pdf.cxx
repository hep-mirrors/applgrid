//
//   @file    generic_pdf.cxx         
//   
//            a generic pdf type - to read in the combinations
//            for the subprocesses from a file, with the file 
//            format determined by Tancredi (details will be 
//            filled in as the implementation becomes more complete)
//
//   @author M.Sutton
// 
//   Copyright (C) 2013 M.Sutton (sutt@cern.ch)    
//
//   $Id: generaic_pdf.cxx, v0.0   Mon 28 Jan 2013 15:40:45 GMT sutt $



#include "appl_grid/appl_pdf.h" 
using namespace appl;

#include "appl_grid/generic_pdf.h"


/// don't want a fortran implementation
//  extern "C" void fgeneric_pdf__(const double* fA, const double* fB, double* H) { 
//    static generic_pdf pdf;
//    pdf.evaluate(fA, fB, H);
//  }

generic_pdf::generic_pdf(const std::string& s) 
  : appl_pdf("generic-"+s), m_initialised(false) {
  /// set up later generic pdf here ...
  
  if ( s!="" ) initialise( s );
  
} 



void generic_pdf::initialise(const std::string& file) { 
  /// first rename me with an appropriate file name
  /// we chose the structure "generic-filename"
  /// so the grid can tell it is generic, and can also 
  /// store the filename 
  /// NB: at some point we will encode the file contents in 
  ///     the grid itself  
  rename("generic-"+file);

  /// add tancredi's code here ...

  std::cout << "generic_pdf::initialise() " << name() << std::endl; 

  m_initialised = true;
} 


void  generic_pdf::evaluate(const double* fA, const double* fB, double* H) {  
  ///  fill this in with tancredi's code ...
  if ( !m_initialised ) return;

}

