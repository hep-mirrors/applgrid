#ifndef __pdf_genlha_h__
#define __pdf_genlha_h__ 1


#include <bits/hhc-process.h>
#include "LHAPDF.h"


//----- used namespaces -----
using namespace std;
using namespace nlo;

class pdf_genlha
  : public pdf_and_coupling_hhc
{
public:
  //   constructor
  inline explicit pdf_genlha(const string & name, unsigned int mem = 0)
    { initPDFset(name);
      initPDF(mem);
    }
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
    return alphasPDF(sqrt(mr2))/6.28318530717958647692;
  }
  
  //   the parton distribution function
  void hadronA(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13]; 
    evolvePDF(x, sqrt(Q2), __f);
    for(int i=-6; i <= 6; i++) f[i] = __f[6+i]/x;
  }
  
  void hadronB(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13];
    evolvePDF(x, sqrt(Q2), __f);
    for(int i=-6; i <= 6; i++) {
      f[i] = __f[6+i]/x; //p
      //f[i] = __f[6-i]/x; //pbar
    }
  }
};


#endif // __pdf_genlha_h__
