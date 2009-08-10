#ifndef __ALEKHIN_H__
#define __ALEKHIN_H__ 1

// lhpdf includes 
#include <lhpdf_alekhin.h>

// lhpdf++ includes 
#include <lhpdf++.h>


namespace lhpdf {
  
  typedef pdfset<&::lhpdf_pdfset_alekhin_100>  pdfset_alekhin_100;
  typedef pdfset<&::lhpdf_pdfset_alekhin_1000> pdfset_alekhin_1000;


  class alekhin_100 : public pdf
  {
  public:
    alekhin_100(unsigned int mem) 
      : pdf(::lhpdf_pdfset_alekhin_100, mem) {}
    
    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_alekhin_100);
    }
  };
  
  
  class alekhin_1000 : public pdf
  {
  public:
    alekhin_1000(unsigned int mem) 
      : pdf(::lhpdf_pdfset_alekhin_1000, mem) {}
    
    //  Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { 
      return ::lhpdf_pdfset_numberpdf(::lhpdf_pdfset_alekhin_1000);
    }
  };

}
#endif
