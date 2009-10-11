/* lhpdf includes */
#include "lhpdf.h"



unsigned int lhpdf_pdfset_nfmax(const lhpdf_pdfset_t *pdfs)
{
  unsigned int i, nfmax = 0U;
  
  for(i = 0; i < 13; i++)
    if(pdfs->threshold[i] >= 0.0) nfmax++;
  return nfmax/2;
}

double lhpdf_pdfset_scale(const lhpdf_pdfset_t *pdfs) {
  return pdfs->q0;
}
  
double lhpdf_pdfset_alfas(const lhpdf_pdfset_t *pdfs, unsigned int mem)
{
  return pdfs -> etype == 2 ? 
    *(pdfs->fparml + mem*(pdfs->nparm) + pdfs->parm) : pdfs->alfasQ;
}

unsigned int lhpdf_pdfset_numberpdf(const lhpdf_pdfset_t *pdfs) {
  return pdfs->nlist;
}

unsigned int lhpdf_pdfset_orderpdf(const lhpdf_pdfset_t *pdfs) {
  return pdfs->evlord;
}

unsigned int lhpdf_pdfset_orderas(const lhpdf_pdfset_t *pdfs) {
  return pdfs->alford;
}

double lhpdf_pdfset_renfac(const lhpdf_pdfset_t *pdfs) {
  return pdfs->mur;
}

const char **lhpdf_pdfset_desc(const lhpdf_pdfset_t *pdfs) {
  return pdfs->descr;
}

double lhpdf_pdfset_qmass(const lhpdf_pdfset_t *pdfs, int nf) 
{
  switch(abs(nf)) {
  case 4: return pdfs->mc; break;
  case 5: return pdfs->mb; break;
  case 6: return pdfs->mt; break;
  default: return 0.0;
  }
}

double lhpdf_pdfset_threshold(const lhpdf_pdfset_t *pdfs, int nf) {
  return pdfs->threshold[nf+6];
}

double lhpdf_pdfset_weightpdf(const lhpdf_pdfset_t *pdfs, unsigned int mem) 
{
  if(pdfs->fw >= 0) 
    return *(pdfs->fparml + mem*(pdfs->nparm) + pdfs->fw);
  else return 1.0/(pdfs->nlist);
}
