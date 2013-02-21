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


#include <iostream>
#include <fstream>

#include "appl_grid/appl_pdf.h" 
using namespace appl;

#include "appl_grid/generic_pdf.h"


/// don't want a fortran implementation
//  extern "C" void fgeneric_pdf__(const double* fA, const double* fB, double* H) { 
//    static generic_pdf pdf;
//    pdf.evaluate(fA, fB, H);
//  }

generic_pdf::generic_pdf(const std::string& s) 
  : appl_pdf(s), 
    m_initialised(false), 
    m_ckmsum(0), m_ckm2(0), H(0)
 {
   /// need to check has the appropriate form, 
   if ( s!="" ) { 
     if ( s.find(".dat")==std::string::npos ) throw exception( std::cerr << "generic_pdf() file " << s << " does not have .dat extension - will not be able to save this grid to file" << std::endl );  
     else initialise( s );
   }
} 



void generic_pdf::initialise(const std::string& filename) { 

  if ( m_initialised ) { 
    std::cerr << "generic_pdf::initialise() already initialised" << std::endl;
    return;
  }


  /// first rename me with an appropriate file name
  /// we chose the structure "filename.dat"
  /// so the grid can tell it is generic, and can also 
  /// store the filename 
  /// NB: at some point we will encode the file contents in 
  ///     the grid itself  
  rename(filename);

  /// add tancredi's code here ...

  std::cout << "generic_pdf::initialise() " << name() << std::endl; 

  m_initialised = true;


  /// tancredi's code 

  debug = true;
  
  currentprocess=-1;
  currentsubprocess=-1;
  nQuark=6;  
  
  make_ckmsum( m_ckmsum ); 
  /// for Wp production 
  make_ckm( m_ckm2, true );  
  /// for Wm production 
  //  make_ckm( m_ckm2, false );  
  
  Print_ckm();  
  
  std::vector<std::string> names;
  names.push_back("tbar");
  names.push_back("bbar");
  names.push_back("cbar");
  names.push_back("sbar");
  names.push_back("ubar");
  names.push_back("dbar");
  names.push_back("g");
  names.push_back("d");
  names.push_back("u");
  names.push_back("s");
  names.push_back("c");
  names.push_back("b");
  names.push_back("t");
  
  std::string* name2 = (&names[0])+6; /// so we can use -6..6 indexing
  for (int i=-6; i<7; i++) {
    //  if (i>0)  name2 = names[i];
    //  if (i<0)  name2 = names[-i];
    //  if (i==0) name2 = names[i];
    //  cout<<i<<" name2= "<<name2<<endl;
    flavname[i] = name2[i];
    iflavour[name2[i]]=i;    
  }
  
  pdfsumtypes1.clear();
  pdfsumtypes2.clear();
  
  Flav1.clear();
  Flav2.clear();
  
  ReadSubprocessSteering(filename);
  
  PrintSubprocess();
  
  flavourtype.clear();
  int ifltype = -99;  
  for(int ifl = -nQuark; ifl <= nQuark; ifl++) {
    if ((ifl== 2)||(ifl== 4)||(ifl== 6))  ifltype =  2; 
    if ((ifl==-2)||(ifl==-4)||(ifl==-6))  ifltype = -2; 
    if ((ifl== 1)||(ifl== 3)||(ifl== 5))  ifltype =  1;
    if ((ifl==-1)||(ifl==-3)||(ifl==-5))  ifltype = -1;
    
    if (ifl==21) ifltype= 0;
    flavourtype[ifl]=ifltype;
  }
  
  if (debug) PrintFlavourMap();
  
  const int nsub = GetSubProcessNumber();
  
  //if (debug) 
  std::cout << "MySubProcess nsub = " << nsub << std::endl;
  H = new double[nsub];
  
} 


void  generic_pdf::evaluate(const double* fA, const double* fB, double* H) {  
  ///  fill this in with tancredi's code ...
  if ( !m_initialised ) return;

}






void generic_pdf::ReadSubprocessSteering(const std::string& fname){
  
  // 
  // read steering deciding on subprocesses
  //
  
  std::cout << "generic_pdf::ReadSubprocessSteering: read subprocess configuration file: " << fname << endl; 
  
  std::ifstream infile(fname.c_str(), std::ios::in);
  if ( !infile ) { // Check open
    std::cerr << "Can't open " << fname << std::endl;
    infile.close();
    exit(1);
  } else {
    std::cout <<" MyData:ReadData: read data file: " << fname << std::endl;
  }
  
  
  char line[256];
  int iproc=0;

  while (infile.good()) {

    if (debug) std::cout << " good: " << infile.good() << " eof: " << infile.eof() << std::endl;

    if (!infile.good()) break;

    infile.getline(line,sizeof(line),'\n');
    std::string mytest = std::string(line);
    //std::cout<< "sizof test= "<< mytest.size() << std::endl;
    
    if (debug) std::cout<< "line= "<< line << "\n";

    if(line[0] != '%'  && mytest.size()!=0) {
      char flav1[100];char flav2[100];
      sscanf(line,"%s %s",flav1,flav2);
      std::cout << " ReadSubProcesses: flav1 = " << flav1 << " flav2 = " << flav2 << endl;
      
      int ifl1=iflavour[flav1];
      int ifl2=iflavour[flav2];

      // std::cout<<"ReadSubProcesses: ifl1 = " << ifl1 << " ifl2 = " << ifl2 << std::endl;
      Flav1[iproc]=ifl1;
      Flav2[iproc]=ifl2;
      iproc++;
      std::string myprocname=string(flav1)+"-"+string(flav2);
      //std::cout << iproc << " process name = " << myprocname << endl;;
      procname.push_back(myprocname);

    }

  }

}



int generic_pdf::decideSubProcess(const int iflav1, const int iflav2, const int nproc)
{
  //
  // 
  // iflav1 change from 0 to 21 (convention for gluons in sherpa)
  // assume that ckm comes with weights

  int    iProcess = -1;

  // std::cout << " iflav1 = " << iflav1 << " iflav2 = " << iflav2 << std::endl;
  int ifl1=flavourtype[iflav1];
  int ifl2=flavourtype[iflav2];
  
  for ( unsigned i=0 ; i<procname.size() ; i++ ) {
    if (iProcess!=-1) continue;
    if (debug) std::cout << " " << i << " name= " << procname[i]
			 << " Flav1, Flav2 = " << Flav1[i] << " " << Flav2[i]
			 << std::endl; 

    if ( Flav1[i]==ifl1 && Flav2[i]==ifl2 ) {
      iProcess=i;
      // std::cout << " iProcess found " << iProcess << endl;
    }
  }
  
  if (debug) std::cout << " iProcess found " << iProcess << std::endl;
  
  if ( iProcess==-1 ) { 
    std::cout << " ***decideSubProcess " << iflav1 << " <> " << iflav2 << " nproc = " << nproc << std::endl; 
  }
  // else
  //    std::cout << " ***decideSubProcess " << iflav1 << " <> " << iflav2 << " iProcess = " << iProcess << std::endl; 
  return iProcess;
}




// get ckm related information

void generic_pdf::Print_ckm() {  
  //
  // print ckm matrix for W+ and W-
  // 
  
  cout << "generic_pdf::Print_ckm" << endl;
  //
  for ( int i=0 ; i<13 ; i++ ) {
    for ( int j=0 ; j<13 ; j++ ) {
      if (m_ckm2[i][j]!=0)
	std::cout << " ckm[" << i << "][" << j << "]\t =\t " << m_ckm2[i][j] << std::endl;
    }
  }
  
  return;
}



void generic_pdf::PrintSubprocess() { 
  
  std::cout << "generic_pdf::PrintSubprocess:" << std::endl;

  // cout<<" size of Flav1: "<<  Flav1.size()<<endl;
  // std::map<int,int>::iterator imap;
  // for (imap = Flav1.begin(); imap!=Flav1.end(); ++imap){
  //  cout<<" "<<imap->first
  //           <<" flav1, flav2= "<<imap->second<<" "<<imap->second<<endl;
  // }

  std::cout << "\t Number of subprocesses: " << procname.size() << std::endl;
  
  for (unsigned i=0; i< procname.size(); i++) {
    std::cout << "\t" << i << " Flav1, Flav2 = " << Flav1[i] << " " << Flav2[i]
	      << "  " << procname[i]
	      << std::endl;
  };

  return;

}



double *generic_pdf::GetGeneralisedPdf(const double *f1, const double *f2){

  //
  // provides generalised PDF per subprocess
  // 
  
  //if (debug) cout<<" GetGeneralisedPdf: currentprocess= "<<currentprocess<<endl;
  if (debug) cout << " GetGeneralisedPdf: " <<endl;
    
  /// shift pointers by nQuark so we can index them using -nQuark .. nQuark  
  f1 += nQuark;
  f2 += nQuark;

  double* ckmsum = m_ckmsum + nQuark;
  
  for(int i = -nQuark; i <= nQuark; i++) {
    
    std::cout << " f1[" << i << "]= " << f1[i]
	      << " f2[" << i << "]= " << f2[i] << " ckmsum[" << i << "]= " << ckmsum[i] << std::endl;
    
    int ifl=flavourtype[i];
    //if (debug) cout<<" iflavour["<<i<<"]= "<<ifl<<endl;
    pdfsumtypes1[ifl] += f1[i]*ckmsum[i];
    pdfsumtypes2[ifl] += f2[i]*ckmsum[i];
    //if (debug) cout<<" pdfsumtypes1["<<ifl<<"]= "<<pdfsumtypes1[ifl]
    //               <<" pdfsumtypes2["<<ifl<<"]= "<<pdfsumtypes2[ifl]<<endl;
  }

  if (debug) std::cout << " Number of subprocesses: " <<  procname.size() << std::endl;

  for ( unsigned iproc=0 ; iproc<procname.size() ; iproc++ ) {

    int ifl1=Flav1[iproc];  
    int ifl2=Flav2[iproc];
   
    H[iproc]=pdfsumtypes1[ifl1]*pdfsumtypes2[ifl2];
    
    if (debug) cout<<iproc<<" name= "<<procname[iproc]<<" ifl1= "<<ifl1<<" ifl2= "<<ifl2
		   <<" pdfsumtypes1["<<ifl1<<"]= "<<pdfsumtypes1[ifl1]
		   <<" pdfsumtypes1["<<ifl2<<"]= "<<pdfsumtypes2[ifl2]
		   <<" H= "<<H[iproc]
		   <<endl;
  };
  
  return H;

}; 




