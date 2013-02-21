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


  generic_pdf(const std::string& s="");

  virtual ~generic_pdf() { 
    if ( H )        delete[] H; 
    if ( m_ckmsum ) delete[] m_ckmsum;
    if ( m_ckm2 ) { 
      for ( int i=0 ; i<13 ; i++ ) if ( m_ckm2[i] ) delete[] m_ckm2[i];
      delete[] m_ckm2;
    }
  } 

  void evaluate(const double* fA, const double* fB, double* H);

  /// additional user defined functions to actually initialise 
  /// based on the input file

  void initialise(const std::string& filename);

  bool initialised() const { return m_initialised; }


public:

  void ReadSubprocessSteering(const std::string& fname);

  void Print_ckm();  

  void PrintSubprocess();

  int  GetSubProcessNumber(){ return procname.size(); }

  double* GetGeneralisedPdf(const double *,const double *); 

  int  decideSubProcess(const int iflav1, const int iflav2, const int nproc);

  void SetSubCurrentProcess(int mypro) { currentsubprocess=mypro; }

  int GetCurrentSubProcess() { 
    if ( currentsubprocess==-1 ) std::cout << " MySubProcess: current subprocess not defined ! " << std::endl;
    return currentsubprocess;
  }

  void SetCurrentProcess(int mypro){ currentprocess=mypro; }
  
  int GetCurrentProcess() { 
    if (currentprocess==-1) cout<<" MySubProcess: current process not defined ! "<<endl;
    return currentprocess;
  }
  
  int GetnQuark() { return nQuark; }

  void SetnQuark(int nq) { nQuark=nq; }

  void PrintFlavourMap() {
    cout<<" print out flavour map "<<endl;
    std::map<int,int>::iterator imap;
    for (imap = flavourtype.begin(); imap!=flavourtype.end(); ++imap){
      cout<<" Flavourtype["<<imap->first<<"]= "<<imap->second<<" "<<endl;
    }
  }


  
private:

  /// this might eventually become a string encoding the grid
  string m_filename;  

  /// has this been initialised yet?
  bool m_initialised;

  /// ckm matrices should they be needed ...
  double*  m_ckmsum;
  double** m_ckm2; 

  bool    debug;

  int     nQuark; //Number of Quarks not yet implemented

  double* H; // generalised PDF defined in GetSubprocess

  // map name to index
  // maps flavour type names to iflavour <-> -6...6
  std::map<string, int> iflavour; 
  
  // maps iflavour to flavour typs names -6...6 <-> name 
  std::map<int, string> flavname;  

  // maps iflavour to types -6...6 <->  -2, -1, 0, 1, 2 
  std::map<int,int> flavourtype;

  // sum of quark flavour belonging to one flavour type up,down,gluon
  std::map<int,double> pdfsumtypes1;
  std::map<int,double> pdfsumtypes2;
  
  // maps iflavour -2, -1, 0, 1, 2 of first or second proton to iprocess
  std::map<int,int> Flav1; 
  std::map<int,int> Flav2;

  std::vector<string> procname; // names of subprocesses

  int currentsubprocess;
  int currentprocess;


};



inline std::ostream& operator<<( std::ostream& s, const generic_pdf& _g ) { 
  return s;
}



#endif  // GENERIC_PDF_H 












