#ifndef __BOTJE_H__
#define __BOTJE_h__ 1

// lhpdf includes 
#include <lhpdf_botje.h>

// lhpdf++ includes 
#include <lhpdf++.h>


namespace lhpdf {
  
  typedef pdfset<&::lhpdf_pdfset_botje_100>  pdfset_botje_100;
  typedef pdfset<&::lhpdf_pdfset_botje_1000> pdfset_botje_1000;


  class botje_100 : public pdf
  {
  public:
    botje_100(unsigned int mem) 
      : pdf(::lhpdf_pdfset_botje_100, mem) {}
    
    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_botje_100);
    }
  };
  
  
  class botje_1000 : public pdf
  {
  public:
    botje_1000(unsigned int mem) 
      : pdf(::lhpdf_pdfset_botje_1000, mem) {}
    
    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_botje_1000);
    }
  };

}
#endif
