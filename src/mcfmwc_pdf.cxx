
#include "appl_grid/appl_pdf.h" 
using namespace appl;

#include "mcfmwc_pdf.h"


// fortran callable wrapper

extern "C" void fmcfmwpc_pdf__(const double* fA, const double* fB, double* H) { 
  static mcfmwpc_pdf pdf;
  pdf.evaluate(fA, fB, H);
}


extern "C" void fmcfmwmc_pdf__(const double* fA, const double* fB, double* H) { 
  static mcfmwmc_pdf pdf;
  pdf.evaluate(fA, fB, H);
}




