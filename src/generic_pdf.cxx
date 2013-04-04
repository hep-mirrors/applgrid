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
//   $Id: generic_pdf.cxx, v0.0   Mon 28 Jan 2013 15:40:45 GMT sutt $


#include <stdlib.h>

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
    H(0)
 {

   /// need to check has the appropriate form, 
   if ( s!="" ) { 
     if ( s.find(".dat")==std::string::npos ) throw exception( std::cerr << "generic_pdf() file " << s << " does not have .dat extension - will not be able to save this grid to file" << std::endl );  
     else initialise( s );
   }

   //debug=false;
  cout<<" initialize generic pdf "<<s<<" debug= "<<debug<<endl;
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

  debug=false;
  
  currentprocess=-1;
  currentsubprocess=-1;
  nQuark=6; // offset to go from 0,..,13 to  -6,..,0.,..,6

  
  //make_ckmsum( m_ckmsum ); 
  /// for Wp production
  // change this to appropriate valuse 
  //make_ckm(true );  
  /// for Wm production 
  //  make_ckm( m_ckm2, false );  
  //Print_ckm();  
  
  std::vector<std::string> names;
  names.push_back("topbar");
  names.push_back("beautybar");
  names.push_back("charmbar");
  names.push_back("strangebar");
  names.push_back("upbar");
  names.push_back("downbar");
  names.push_back("gluon");
  names.push_back("down");
  names.push_back("up");
  names.push_back("strange");
  names.push_back("charm");
  names.push_back("beauty");
  names.push_back("top");
  
  std::string* name2 = (&names[0])+6; /// so we can use -6..6 indexing
  for (int i=-nQuark; i<=nQuark; i++) {
    //  if (i>0)  name2 = names[i];
    //  if (i<0)  name2 = names[-i];
    //  if (i==0) name2 = names[i];
    //cout<<i<<" name2= "<<name2<<endl;
    flavname[i] = name2[i];
    iflavour[name2[i]]=i;    
    if (debug)
    cout<<"generic_pdf "<<" iflavour["<<name2[i]<<"]= "<<iflavour[name2[i]]<<endl;
  }
  
  //pdfsumtypes1.clear();
  //pdfsumtypes2.clear();
  
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
    
    if (ifl==0) ifltype= 0;
    flavourtype[ifl]=ifltype;
  }

  //this->MakeCkm();
 //if (debug) PrintFlavourMap();
  
  const int nsub = GetSubProcessNumber();
  
  if (debug) 
   std::cout << "generic_pdf::initialize nsub = " << nsub << std::endl;
  H = new double[nsub];
  
} 


void  generic_pdf::evaluate(const double* fA, const double* fB, double* H) {  
 //  fill this in with tancredi's code ...
 if ( !m_initialised ) {
  cout<<"  generic_pdf::evaluate not initialized "<<endl;
  return; 
 }

 //need to move this out to initialize 
 this->MakeMyCkm();

 /*
 if (debug){
  cout<<"generic_pdf print my ckm matrix "<<endl;
  this->Print_myckm();  
 }
 */
 if (debug) cout << "generic_pdf:evaluate " <<endl;
 // shift pointers by nQuark so we can index them using -nQuark .. nQuark  
 //fA += nQuark;
 //fB += nQuark;
 //double* ckmsum = m_ckmsum + nQuark;
  
 //this->PrintFlavourMap();

 /*
 if (debug) {
  for(int i = -nQuark; i <= nQuark; i++) {
   std::cout << " fA[" << i << "]= " << fA[i+nQuark]
             << " fB[" << i << "]= " << fB[i+nQuark] 
             << std::endl;
  }
 }  
 */
 // reset pdf sums per flavour -2,-1,0,1,2 downbar, upbar, gluon, up, down
 pdfA.clear();
 pdfB.clear();
 for(int i = -2; i <= 2; i++) {
  pdfA[i]=0.; pdfB[i]=0.;
 }
 pdfA[0]=fA[0+nQuark]; // gluon
 pdfB[0]=fA[0+nQuark];

/*
 if (debug) {
  for(int i = -2; i <= 2; i++) {
   cout<<" pdfA["<<i<<"]= "<<pdfA[i]
       <<" pdfB["<<i<<"]= "<<pdfB[i]
       <<" m_ckmsum["<<i<<"]= "<<m_ckmsum[nQuark+i]
       <<endl;
  }
 }
*/
 std::vector<double> myckmsum = std::vector<double>(myckm2.size(),0);
 for ( unsigned i=0 ; i<myckm2.size() ; i++ ) { 
  for ( unsigned j=0 ; j<myckm2[i].size() ; j++ ) myckmsum[i] += myckm2[i][j]; 
 }


 /*
 for(int i=-6; i <=6; i++) {
  for(int i=-6; i <=6; i++) {
   myckm2[i+nQuark][j+nQuark]
  }
 }
 */
 for(int i=-6; i <=6; i++) {
  int j=flavourtype[i];
  if (j==0) continue;
  pdfA[j] += fA[nQuark + i]*myckmsum[nQuark+j];
  pdfB[j] += fB[nQuark + i]*myckmsum[nQuark+j];
  if (debug)
    cout<<" i= "<<i<<" j= "<<j
      <<" fA["<<i<<"]= "<<fA[nQuark + i]
      <<" fB["<<i<<"]= "<<fB[nQuark + i]
      <<" pdfA["<<j<<"]= "<<pdfA[j]<<" pdfB["<<j<<"]= "<<pdfB[j]
      <<" myckmsum["<<nQuark+j<<"]= "<<myckmsum[nQuark+j]
      <<endl;
 }

 /*
 if (debug) {
  for(int i = -2; i <= 2; i++) {
   cout<<" pdfA["<<i<<"]= "<<pdfA[i]
       <<" pdfB["<<i<<"]= "<<pdfB[i]
       <<endl;
  }
 }
 */
 /*
 double S12=0., S21=0.;
 for (int i1 = -3; i1 <= -1; i1 += 2){
  for(int i2 =  2; i2 <=  4; i2 += 2){
    cout<<" i1= "<<i1+nQuark<<" i2= "<<i2+nQuark<<endl;
    S12 += fA[i1+nQuark]*fB[i2+nQuark]*m_ckm2[i1+nQuark][i2+nQuark];
    S21 += fA[i2+nQuark]*fB[i1+nQuark]*m_ckm2[i2+nQuark][i1+nQuark];
   }
 }    
 */
 for ( unsigned iproc=0 ; iproc<procname.size() ; iproc++ ) {
  int ifl1=Flav1[iproc];  
  int ifl2=Flav2[iproc];

  H[iproc]=pdfA[ifl1]*pdfB[ifl2];
  if (ifl1==ifl2)   H[iproc]*=2.; // symetric contributions are counted twice

  if (debug) cout<<iproc<<" ifl1= "<<ifl1<<" ifl2= "<<ifl2
         	 <<" pdfA["<<ifl1<<"]= "<<pdfA[ifl1]
	         <<" pdfB["<<ifl2<<"]= "<<pdfB[ifl2]
	         <<" H= "<<H[iproc]
                 <<" name= "<<procname[iproc]
	         <<endl;
 }



 return;

}; 

void generic_pdf::ReadSubprocessSteering(const std::string& fname){
  
  // 
  // read steering deciding on subprocesses
  //
  
  if (debug)
  std::cout << "generic_pdf::ReadSubprocessSteering: read subprocess configuration file: " << fname << endl; 
  
  std::ifstream infile(fname.c_str(), std::ios::in);
  if ( !infile ) { // Check open
    std::cerr << "Can't open " << fname << std::endl;
    infile.close();
    exit(1);
  } else {
    std::cout <<" generic_pdf:ReadData: read data file: " << fname << std::endl;
  }
  
  
  char line[256];
  int iproc=0;

  while (infile.good()) {

    //if (debug) std::cout << " good: " << infile.good() << " eof: " << infile.eof() << std::endl;

    if (!infile.good()) break;

    infile.getline(line,sizeof(line),'\n');
    std::string mytest = std::string(line);
    //std::cout<< "sizof test= "<< mytest.size() << std::endl;
    
    //if (debug) std::cout<< "line= "<< line << "\n";

    if(line[0] != '%'  && mytest.size()!=0) {
      char flav1[100];char flav2[100];
      sscanf(line,"%s %s",flav1,flav2);
      std::cout << " ReadSubProcesses: flav1 = " << flav1 << " flav2 = " << flav2 << endl;
      
      int ifl1=iflavour[flav1];
      int ifl2=iflavour[flav2];

      if (debug)
      std::cout<<"ReadSubProcesses: ifl1 = " << ifl1 << " ifl2 = " << ifl2 << std::endl;
      Flav1[iproc]=ifl1;
      Flav2[iproc]=ifl2;
      iproc++;
      std::string myprocname=string(flav1)+"-"+string(flav2);
      if (debug) std::cout << iproc << " process name = " << myprocname << endl;;
      procname.push_back(myprocname);

    }

  }

}



int generic_pdf::decideSubProcess(const int iflav1, const int iflav2)
{
  // 
  // iflav1 change from 0 to 21 (convention for gluons in sherpa)
  // assume that ckm comes with weights

  int    iProcess = -1;
 
  if (debug) std::cout << "generic_pdf::decideSubProces: " << std::endl; 
  if (debug) std::cout << " iflav1 = " << iflav1 << " iflav2 = " << iflav2 << std::endl;
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
  
  if (debug) std::cout << "generic_pdf:decideSubprocess iProcess found " << iProcess << std::endl;
  
  if ( iProcess==-1 ) { 
    std::cout << "generic_pdf:decideSubprocess " << iflav1 << " <> " << iflav2 << std::endl; 
  }
  // else
  //    std::cout << " ***decideSubProcess " << iflav1 << " <> " << iflav2 << " iProcess = " << iProcess << std::endl; 
  currentsubprocess=iProcess;
  return iProcess;
}


void generic_pdf::MakeMyCkm() { 
  //
  // make by hand CKM matrix
  //
    myckm2 = std::vector<std::vector<double> >(13, std::vector<double>(13,0));

    myckm2[3][8]  =   0.049284000000000001417976847051249933 ;
    myckm2[8][3]  =   0.049284000000000001417976847051249933 ;
    
    myckm2[5][8]  =   0.950624999999999942268402719491859898 ;
    myckm2[8][5]  =   0.950624999999999942268402719491859898 ;
    
    myckm2[5][10] =   0.049284000000000001417976847051249933 ;
    myckm2[10][5] =   0.049284000000000001417976847051249933 ;
    
    myckm2[3][10] =   0.950624999999999942268402719491859898 ;
    myckm2[10][3] =   0.950624999999999942268402719491859898 ;

    myckm2[4][9] =   0.049284000000000001417976847051249933 ;
    myckm2[9][4] =   0.049284000000000001417976847051249933 ;
    
    myckm2[7][4] =   0.950624999999999942268402719491859898 ;
    myckm2[4][7] =   0.950624999999999942268402719491859898 ;
    
    myckm2[7][2] =   0.049284000000000001417976847051249933 ;
    myckm2[2][7] =   0.049284000000000001417976847051249933 ;
    
    myckm2[9][2] =   0.950624999999999942268402719491859898 ;
    myckm2[2][9] =   0.950624999999999942268402719491859898 ;

  return;
}

// get ckm related information

/*
void generic_pdf::Print_ckm() {  
  //
  // print ckm matrix for W+ and W-
  // 
  
  cout << "generic_pdf::Print_ckm = " << endl;
  cout<<" ckm size= "<<m_ckm2.size()<<endl;
  if (m_ckm2.size()<=0) return;

  //
  for ( int i=0 ; i<13 ; i++ ) {
    for ( int j=0 ; j<13 ; j++ ) {
     if (m_ckm2[i][j]!=0)
      std::cout << " ckm[" << i << "][" << j << "]\t =\t " << m_ckm2[i][j] << std::endl;
    }
  }
  
  return;
}
*/
void generic_pdf::Print_myckm() {  
  // 
  
  cout << "generic_pdf::Print_myckm = " << endl;
  cout<<" ckm size= "<<myckm2.size()<<endl;
  if (myckm2.size()<=0) return;

  //
  for ( int i=0 ; i<13 ; i++ ) {
    for ( int j=0 ; j<13 ; j++ ) {
     if (myckm2[i][j]!=0)
      std::cout << " myckm2[" << i-nQuark << "][" << j-nQuark << "]\t =\t " << myckm2[i][j] << std::endl;
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






