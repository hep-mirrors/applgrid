#ifndef __MRST2001_H__
#define __MRST2001_H__ 1

// lhpdf includes 
#include <lhpdf_mrst2001.h>

// lhpdf++ includes 
#include <lhpdf++.h>


namespace lhpdf {
  
  typedef pdfset<&::lhpdf_pdfset_mrst2001> pdfset_mrst2001;

  class mrst2001 : public pdf
  {
  public:
    mrst2001(unsigned int mem) 
      : pdf(::lhpdf_pdfset_mrst2001, mem) {}

    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_mrst2001);
    }
  };
}
#endif
