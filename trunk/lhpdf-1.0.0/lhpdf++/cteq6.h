#ifndef __CTEQ6_H__
#define __CTEQ6_H__ 1

// lhpdf includes 
#include <lhpdf_cteq6.h>

// lhpdf++ includes 
#include <lhpdf++.h>


namespace lhpdf {
  
  typedef pdfset<&::lhpdf_pdfset_cteq6>  pdfset_cteq6;

  class cteq6 : public pdf
  {
  public:
    cteq6(unsigned int mem) 
      : pdf(::lhpdf_pdfset_cteq6, mem) {}
    
    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_cteq6);
    }
  };
}
#endif
