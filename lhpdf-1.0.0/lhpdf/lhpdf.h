/**
 * @file lhpdf.h
 */
#ifndef __lhpdf_h__
#define __lhpdf_h__ 1

/* Standard C includes */
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef double LHAPDF_FLOAT;   

/** 
 * \struct lhpdf_pdfset_t
 *   \brief Structure to store the static information of the pdfset
 *
 *   bla-bla-bla
 */
typedef struct {
  /* Description of the pdf set */
  const char **descr;

  /* Evolution code : evlcteq = 0, qcdnum = 1 */
  unsigned int evlcode;

  /* Order ov the evolution : lo = 0, nlo = 1, nnlo = 2 */
  unsigned int evlord;
  
  /* fitting scale & the ratio of the renorm. and fact. scales */
  double q2fit, mur;
  
  /* Name of the grid */
  const char *gridname;
  
  /* Limits */
  double xmin,xmax,qmin,qmax;
  
  /*   size of the grid */
  unsigned int nx, nq; 


  /* Order ov the evolution : lo = 0, nlo = 1, nnlo = 2 */
  unsigned int alford;
  
  /** \brief Type of the evolution : Fixed = 1, Variable = 2 */
  unsigned int etype;
  
  /** \brief Method of the evolution: Internal = 0, EvolCode = 1 */
  unsigned int method;
  
  /** \brief The \f$\alpha_s(Q_0)\f$ is defined by one of the parameter of the pdf set.
   *  This variable is the index of that parameter which gives the \f$\alpha_s\f$ at the fitting scale.
   */
  unsigned int parm;
  
  /** \brief The \f$\alpha_s\f$ strong coupling at the fitting scale */
  double alfasQ; 
 
  /** \brief The \f$Q_0\f$ fitting scale */
  double q0;
  
  /** charm mass */
  double mc, mb, mt;
  
  
  /* Functions */
  void (*init_pdf_func)(const LHAPDF_FLOAT *, double, double *);

  /* Threshold */
  const double threshold[13];
  
  /* Fw >= 0 if this pdf is a weighted pdf */
  int fw;

  
  /* Number of list */
  unsigned int nlist;
  
  /* Number of parameters */
  unsigned int nparm;
  
  /* Parameters */
  const LHAPDF_FLOAT *fparml;
} lhpdf_pdfset_t;

/**   Returns the number of the available pdf in the pdf set  */ 
unsigned int lhpdf_pdfset_numberpdf(const lhpdf_pdfset_t *);

/**   Returns the evolution order of the pdf set  */
unsigned int lhpdf_pdfset_orderpdf(const lhpdf_pdfset_t *);

/**   Returns the evolution order of the \f$\alpha_s\f$ */
unsigned int lhpdf_pdfset_orderas(const lhpdf_pdfset_t *);

/**   Returns the number of the active flavours  */
unsigned int lhpdf_pdfset_nfmax(const lhpdf_pdfset_t *);

/**   Returns the ratio of the renormalization and factorization scales */
double lhpdf_pdfset_renfac(const lhpdf_pdfset_t *);

/**   Returns the fitting scale  */
double lhpdf_pdfset_scale(const lhpdf_pdfset_t *);

/**   Returns the alphas at the fitting scale  */
double lhpdf_pdfset_alfas(const lhpdf_pdfset_t *, unsigned int);

/**   Returns the description  */
const char ** lhpdf_pdfset_desc(const lhpdf_pdfset_t *);

/**   Returns the quark masses  */
double lhpdf_pdfset_qmass(const lhpdf_pdfset_t *, int);

/**   Returns the flavour threshold in the evolution */
double lhpdf_pdfset_threshold(const lhpdf_pdfset_t *, int);

/**   Returns weight of the pdfs */
double lhpdf_pdfset_weightpdf(const lhpdf_pdfset_t *, unsigned int); 


typedef struct {
  /* Constant pointer to the pdf set */
  const lhpdf_pdfset_t *pdfset;

  /* The current pdf */
  const LHAPDF_FLOAT *fparm;
  
  /* Store the grid and the evolution code dependent parameters */
  void *grid;  
} lhpdf_pdf_t;


/**  Constructor and destructor */
lhpdf_pdf_t * lhpdf_pdf_alloc(const lhpdf_pdfset_t *pdfs, unsigned int mem);
void lhpdf_pdf_free(lhpdf_pdf_t *pdf);

/**  Evolve the \f$\alpha_s\f$  */
double lhpdf_pdf_evolveas(lhpdf_pdf_t *pdf, double);

/**  Evolve the parton distribution function  */
void lhpdf_pdf_evolvepdf(const lhpdf_pdf_t *, double, double, double *);


/**   Returns the evolution order of the pdf set  */
unsigned int lhpdf_pdf_orderpdf(const lhpdf_pdf_t *);

/**   Returns the evolution order of the \f$\alpha_s\f$ */
unsigned int lhpdf_pdf_orderas(const lhpdf_pdf_t *);

/**   Returns the number of the active flavours  */
unsigned int lhpdf_pdf_nfmax(const lhpdf_pdf_t *);

/**   Returns the ratio of the renormalization and factorization scales */
double lhpdf_pdf_renfac(const lhpdf_pdf_t *);

/**   Returns the fitting scale  */
double lhpdf_pdf_scale(const lhpdf_pdf_t *);

/**   Returns the alphas at the fitting scale  */
double lhpdf_pdf_alfas(const lhpdf_pdf_t *);

/**   Returns the description  */
const char ** lhpdf_pdf_desc(const lhpdf_pdf_t *);

/**   Returns the quark masses  */
double lhpdf_pdf_qmass(const lhpdf_pdf_t *, int);

/**   Returns the flavour threshold in the evolution */
double lhpdf_pdf_threshold(const lhpdf_pdf_t *, int);

/**   Returns weight of the pdfs */
double lhpdf_pdf_weightpdf(const lhpdf_pdf_t *); 





typedef struct {
  double (*evolve_alfas)(double);
  void (*evolve_pdf)(const lhpdf_pdf_t *, double, double, double *);
  void (*grid_alloc)(lhpdf_pdf_t *);
  void (*grid_free)(void *);
} lhpdf_evlcode_t;

extern const lhpdf_evlcode_t *lhpdf_evlcode_evlcteq;
extern const lhpdf_evlcode_t *lhpdf_evlcode_qcdnum;

double lhpdf_log_pade(double x, const LHAPDF_FLOAT *fp, unsigned int fn, ...);

double lhpdf_x_taylor(double x, const LHAPDF_FLOAT *fp, double fpow, 
		       unsigned int fn, ...);


double lhpdf_cteq6_ratio(double x, const LHAPDF_FLOAT *fp, 
			  unsigned int m1, unsigned int m2, 
			  unsigned int m3, unsigned int m4, 
			  unsigned int m5);
  


#ifdef __cplusplus
}
#endif


#endif
