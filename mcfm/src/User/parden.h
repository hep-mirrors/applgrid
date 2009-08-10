// emacs: this is -*- c++ -*-
/**************************************************************************
 **
 **   File:         parden.h  
 **
 **   Description:    
 **                   
 **                   
 ** 
 **   Author:       M.Sutton  
 **
 **   Created:      Thu Oct  4 12:24:26 BST 2007
 **   Modified:     
 **                   
 **                   
 **
 **************************************************************************/ 


#ifndef __PARDEN_H
#define __PARDEN_H

extern "C" { 
void setup_();
void initpdfset_(const char* name);
void initpdfsetbyname_(const char* name);
void initpdfsetbynamem_(int& , const char* name);
void initpdf_(int& set);
void numberpdf_(int& set);
void getrenfac_(double& Q);
void getalfas_(int& n, double& alfas, double& Q);
double alphaspdf_(const double& );
void evolvepdf_(double& x, double& Q, double* f);
void getthreshold_(int& imem, double& Q);

// void parden_(int& nloop);
// void pdflo_(const double& x, const double& Q, double* f);
// void crossing_(const double& xk, double* f);
}



#endif  /* __PARDEN_H */










