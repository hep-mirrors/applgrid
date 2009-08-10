
/* Standard C includes */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/* lhpdf includes */
#include "lhpdf.h"


double lhpdf_log_pade(double x, const LHAPDF_FLOAT *fp, unsigned int fn, ...)
{
  if(x < 0.9999999) {
    va_list ap;
    double poly;
    unsigned int i, m[6];
    
    va_start(ap, fn);
    if(fn != 6) {
      fprintf(stderr, "LHPDF : wrong input paramets of the function log_pade\n"); 
      exit(-1);
    }
    
    for(i = 0; i < fn; i++)
      m[i] = va_arg(ap, unsigned int);
    va_end(ap);
    
    poly = fp[m[1]]*log(x)+fp[m[2]]*log(1.0-x)
      +    fp[m[3]]*x + fp[m[5]]*log(1.0+x*exp(fp[m[4]]));
    
    return fp[m[0]]*exp(poly);
  } else return 0.0;
}

double lhpdf_x_taylor(double x, const LHAPDF_FLOAT *fp, double fpow, unsigned int fn, ...)
{
  va_list ap;
  double poly = 1.0;
  unsigned int j, m1, m2, m3, m4;
  
  va_start(ap, fn);
  m1 = va_arg(ap, unsigned int);
  m2 = va_arg(ap, unsigned int);
  m3 = va_arg(ap, unsigned int);
  
  for(j=4; j <= fn; j++) {
    m4 = va_arg(ap, unsigned int);
    poly += fp[m4]*pow(x, (j-3)/fpow);
  }
  
  va_end(ap);
  poly *= fp[m1];
  
  return pow(x, fp[m2])*pow(1.0-x, fp[m3])*poly;
}

double lhpdf_cteq6_ratio(double x, const LHAPDF_FLOAT *fp,
			  unsigned int m1, unsigned int m2, 
			  unsigned int m3, unsigned int m4, 
			  unsigned int m5)
{
  double poly, b0 = 10.0;
  
  poly  = exp(fp[m1])*pow(x, fp[m2]-1.0)*pow(1.0-x, fp[m3]);
  poly += (1.0 + fp[m4]*x)*pow(1.0-x, fp[m5]);

  if(poly > b0) return poly;
  else if(poly < -b0) return 0.0;
  else return poly + log(1.0+exp(-b0*poly)-exp(-b0))/b0;
}
 

lhpdf_pdf_t *
lhpdf_pdf_alloc(const lhpdf_pdfset_t *pdfs, unsigned int mem)
{
  lhpdf_pdf_t *s = malloc(sizeof(lhpdf_pdf_t));
  s->pdfset = pdfs;

  if(mem >= pdfs->nlist) {
    fprintf(stderr, "Maximum number of PDFs in list exceeded: %u >= %u\n", mem,pdfs->nlist);
    fprintf(stderr, "Returning most likely PDF\n");
    mem = 0;
  }
  
  s->fparm = (s->pdfset)->fparml + mem*((s->pdfset)->nparm);
  
  switch((s->pdfset)->evlcode) {
  case 0: (*(lhpdf_evlcode_evlcteq->grid_alloc))(s); break;
  case 1: (*(lhpdf_evlcode_qcdnum->grid_alloc))(s); break;
  }

  return s;
}

void lhpdf_pdf_free(lhpdf_pdf_t *pdf)
{
  switch((pdf->pdfset)->evlcode) {
  case 0: (*(lhpdf_evlcode_evlcteq->grid_free))(pdf->grid);
  case 1: (*(lhpdf_evlcode_qcdnum->grid_free))(pdf->grid);
  }
  
  if(pdf) free(pdf);
}

#define b0 1.2202
#define b1 0.4897
#define b2 0.1913

double lhpdf_pdf_evolveas(lhpdf_pdf_t *pdf, double Q)
{
  double L, As, alphas = 0.0;
  const lhpdf_pdfset_t *ps = pdf->pdfset;

  if(ps->method == 0) {
    if((ps->etype == 1) || (ps->etype == 2)) {
      L = log(Q/ps->q0);
      
      if(ps->etype == 2) As = pdf->fparm[ps->parm];
      else As = ps->alfasQ;
      
      switch(ps->alford) { 
      case 0: L = b0*L; break;
      case 1: L = (b0 + As*b1)*L; break;
      case 2: L = (b0+As*b1+As*As*b2)*L - 0.25*As*As*b0*b1*L*L; break;
      default: 
	fprintf(stderr, "LHAPDF: The order of the alfas must be less than 2!\n");
	exit(-1);
      }
      alphas = As/(1.0 + As*L);
    }
  } else if(ps->method == 1) {
    switch((pdf->pdfset)->evlcode) {
    case 0: alphas = (*(lhpdf_evlcode_evlcteq->evolve_alfas))(Q);
    case 1: alphas = (*(lhpdf_evlcode_qcdnum->evolve_alfas))(Q);
    }
  } else {
    fprintf(stderr, "LHAPDF: Unknown evolution method for alpha_qcd!");
    exit(-1);
  }
  
  return alphas;
}

void lhpdf_pdf_evolvepdf(const lhpdf_pdf_t *pdf,
			  double x, double q, double *f)
{
  switch((pdf->pdfset)->evlcode) {
  case 0: (*(lhpdf_evlcode_evlcteq->evolve_pdf))(pdf, x, q, f); break;
  case 1: (*(lhpdf_evlcode_qcdnum->evolve_pdf))(pdf, x, q, f); break;
  }
}



double lhpdf_pdf_alfas(const lhpdf_pdf_t *pdf) {
  return (pdf->pdfset)->etype == 2 ? pdf -> fparm[(pdf->pdfset)->parm] : (pdf->pdfset)->alfasQ;
}

double lhpdf_pdf_scale(const lhpdf_pdf_t *pdf) {
  return (pdf->pdfset)->q0;
}


unsigned int lhpdf_pdf_nfmax(const lhpdf_pdf_t *pdf)
{
  unsigned int i, nfmax = 0U;
  
  for(i = 0; i < 13; i++)
    if((pdf->pdfset)->threshold[i] >= 0.0) nfmax++;
  return nfmax/2;
}

unsigned int lhpdf_pdf_orderpdf(const lhpdf_pdf_t *pdf) {
  return (pdf->pdfset)->evlord;
}

unsigned int lhpdf_pdf_orderas(const lhpdf_pdf_t *pdf) {
  return (pdf->pdfset)->alford;
}

double lhpdf_pdf_renfac(const lhpdf_pdf_t *pdf) {
  return (pdf->pdfset)->mur;
}

const char **lhpdf_pdf_desc(const lhpdf_pdf_t *pdf) {
  return (pdf->pdfset)->descr;
}

double lhpdf_pdf_qmass(const lhpdf_pdf_t *pdf, int nf) 
{
  switch(abs(nf)) {
  case 4: return (pdf->pdfset)->mc; break;
  case 5: return (pdf->pdfset)->mb; break;
  case 6: return (pdf->pdfset)->mt; break;
  default: return 0.0;
  }
}

double lhpdf_pdf_threshold(const lhpdf_pdf_t *pdf, int nf) {
  return (pdf->pdfset)->threshold[nf+6];
}

double lhpdf_pdf_weightpdf(const lhpdf_pdf_t *pdf) 
{
  if((pdf->pdfset)->fw >= 0) 
    return pdf->fparm[(pdf->pdfset)->fw];
  else return 1.0/((pdf->pdfset)->nlist);
}


