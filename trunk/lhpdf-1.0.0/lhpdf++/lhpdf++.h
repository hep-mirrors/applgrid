#ifndef __LHPDFXX_H__
#define __LHPDFXX_H__ 1

#include <lhpdf.h>


namespace lhpdf {


  template<const ::lhpdf_pdfset_t **_PdfSet>
  struct pdfset 
  {
    ///   conversion
    static const ::lhpdf_pdfset_t * parton_function_set;
    
    ///   Returns the number of the available pdf in the pdf set   
    static unsigned int number_of_pdf() { return ::lhpdf_pdfset_numberpdf(*_PdfSet);}

    ///   Returns the evolution order of the pdf set  
    static unsigned int order_pdf() { return ::lhpdf_pdfset_orderpdf(*_PdfSet);}
    
    ///   Returns the evolution order of the \f$\alpha_s\f$ 
    static unsigned int order_as() { return ::lhpdf_pdfset_orderas(*_PdfSet);}
    
    ///   Returns the number of the active flavours  
    static unsigned int nfmax() { return ::lhpdf_pdfset_nfmax(*_PdfSet);}
    
    ///   Returns the ratio of the renormalization and factorization scales 
    static double renfac() { return ::lhpdf_pdfset_renfac(*_PdfSet);}
    
    ///   Returns the fitting scale  
    static double scale() { return ::lhpdf_pdfset_scale(*_PdfSet);}
    
    ///   Returns the alphas at the fitting scale  
    static double alfas(unsigned int mem) { return ::lhpdf_pdfset_alfas(*_PdfSet, mem);}
    
    ///   Returns the description  
    static const char **desc() { return ::lhpdf_pdfset_desc(*_PdfSet);}
    
    ///   Returns the quark masses  
    static double qmass(int i) { return ::lhpdf_pdfset_qmass(*_PdfSet, i);}
    
    ///   Returns the flavour threshold in the evolution 
    static double threshold(int i) { return ::lhpdf_pdfset_threshold(*_PdfSet, i);}
    
    ///   Returns weight of the pdfs 
    static double weightpdf(unsigned int mem) { return ::lhpdf_pdfset_weightpdf(*_PdfSet, mem);} 
  };

  template<const ::lhpdf_pdfset_t **_PdfSet>
  const ::lhpdf_pdfset_t * pdfset<_PdfSet>::parton_function_set = *_PdfSet;
  

  class pdf
  {
  public:    
    ///   constructors
    pdf(const ::lhpdf_pdfset_t *s, unsigned mem) 
      : _M_pdf(::lhpdf_pdf_alloc(s, mem)) {}

    ///   destructor
    ~pdf() { ::lhpdf_pdf_free(_M_pdf);}
    
    ///  Evolve the parton distribution function 
    void operator()(double x, double q, double *f) const {
      ::lhpdf_pdf_evolvepdf(_M_pdf, x, q, f);
    }
    
    /// Evolve the \f$\alpha_s\f$  
    double operator()(double q) const {
      return ::lhpdf_pdf_evolveas(_M_pdf, q);
    }
    
    ///   Returns the evolution order of the pdf set 
    unsigned int order_pdf() const { return ::lhpdf_pdf_orderpdf(_M_pdf);}
    
    ///   Returns the evolution order of the \f$\alpha_s\f$
    unsigned int order_as() const { return ::lhpdf_pdf_orderas(_M_pdf);}
    
    ///   Returns the number of the active flavours
    unsigned int nfmax() const { return ::lhpdf_pdf_nfmax(_M_pdf);}
    
    ///   Returns the ratio of the renormalization and factorization scales
    double renfac() const { return ::lhpdf_pdf_renfac(_M_pdf);}
    
    ///   Returns the fitting scale
    double scale() const { return ::lhpdf_pdf_scale(_M_pdf);}
    
    ///   Returns the alphas at the fitting scale
    double alfas() const { return ::lhpdf_pdf_alfas(_M_pdf);}
    
    ///   Returns the description
    const char ** desc() const { return ::lhpdf_pdf_desc(_M_pdf);}
    
    ///   Returns the quark masses
    double qmass(int i) const { return ::lhpdf_pdf_qmass(_M_pdf, i);}
    
    ///   Returns the flavour threshold in the evolution 
    double threshold(int i) const { return ::lhpdf_pdf_threshold(_M_pdf, i);}
    
    ///   Returns weight of the pdfs 
    double weightpdf() const { return ::lhpdf_pdf_weightpdf(_M_pdf);} 
    
  private:
    ::lhpdf_pdf_t *_M_pdf;

    //  no copy and assignment
    pdf(const pdf&) {}
    void operator=(const pdf&) {}
  };
  
  
}

#endif
