#ifndef __FERMI2002_H__
#define __FERMI2002_h__ 1

// lhpdf includes 
#include <lhpdf_fermi2002.h>

// lhpdf++ includes 
#include <lhpdf++.h>


namespace lhpdf {
  
  typedef pdfset<&::lhpdf_pdfset_fermi2002_100>  pdfset_fermi2002_100;
  typedef pdfset<&::lhpdf_pdfset_fermi2002_1000> pdfset_fermi2002_1000;


  class fermi2002_100 : public pdf
  {
  public:
    fermi2002_100(unsigned int mem) 
      : pdf(::lhpdf_pdfset_fermi2002_100, mem) {}
    
    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_fermi2002_100);
    }
  };
  
  
  class fermi2002_1000 : public pdf
  {
  public:
    fermi2002_1000(unsigned int mem) 
      : pdf(::lhpdf_pdfset_fermi2002_1000, mem) {}
    
    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_fermi2002_1000);
    }
  };

}
#endif
