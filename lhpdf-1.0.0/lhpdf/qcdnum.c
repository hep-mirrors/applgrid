/* Standard C includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* lhpdf includes */
#include "lhpdf.h"


typedef int ftnlen;

/* QCDNUM f77 functions & subroutines */
extern int qninit_(void);
extern int grxdef_(int *, double *);
extern int evlsea_(char *, int *, int *, int *, int *, ftnlen);
extern int evolbp_(char *, int *, int *, ftnlen);
extern int evolcp_(char *,  int *, int *, int *, ftnlen);
extern int qnbook_(int *, char *, ftnlen);
extern int evolsg_(int *, int *, int *);
extern int grqdef_(int *, double *, double *);
extern int grqinp_(double *, int *);
extern int qniset_(char *, int *, ftnlen);
extern int qthres_(double *, double *);
extern int qnrset_(char *, double *, ftnlen);
extern int qnlset_(char *, int *, ftnlen);
extern int qnfilw_(int *, int *);
extern int qnpset_(char *, int *, int *, double *, ftnlen); 
extern int evolnm_(char *, int *, int *, int *, ftnlen);
extern int opengridfile_(char *, ftnlen);
extern int nflget_(int *);
extern int iqfromq_(double *);
extern double xfromix_(int *);



static void qcdnum_initevolvecode(const lhpdf_pdfset_t *pdfs)
{
  static int c__0 = 0, c__1 = 1, c__2 = 2, c__3 = 3, c__5 = 5;
  static int c__4 = 4, c__6 = 6, c__7 = 7, c__8 = 8;
  static int c_true = 1, c_false = 0;

  int evlord, nx = pdfs->nx, nq = pdfs->nq;
  double d_1, mc2, tc2, mb2, tb2, mt2;
  double xmin = pdfs->xmin, qmin = pdfs->qmin, qmax = pdfs->qmax;
  double q2fit = pdfs->q2fit, mc = pdfs->mc, mb = pdfs->mb, mt = pdfs->mt;
  char *gridname = (char *) pdfs->gridname;

  qninit_();
  
  evlord = pdfs->evlord + 1;
  qniset_("ORDER", &evlord, (ftnlen) 5);
  grxdef_(&nx, &xmin);
  grqdef_(&nq, &qmin, &qmax);
  grqinp_(&q2fit, &c__1);
  
  mc2 = mc*mc;
  tc2 = (pdfs->threshold[10])*(pdfs->threshold[10]);
  d_1 = tc2 - 1e-4;
  qnrset_("CMASS", &mc, (ftnlen) 5);
  qnrset_("MCALF", &mc, (ftnlen) 5);
  grqinp_(&tc2, &c__1);
  grqinp_(&d_1, &c__1);
  
  mb2 = mb*mb;
  tb2 = (pdfs->threshold[11])*(pdfs->threshold[11]);
  d_1 = tb2 - 1e-4;
  qnrset_("BMASS", &mb, (ftnlen) 5);
  qnrset_("MBALF", &mb, (ftnlen) 5);
  grqinp_(&tb2, &c__1);
  grqinp_(&d_1, &c__1);
  
  mt2 = mt*mt;
  qnrset_("TMASS", &mt, (ftnlen) 5);
  qnrset_("MTALF", &mt, (ftnlen) 5);
  qthres_(&tc2, &tb2);
    
  qnbook_(&c__2, "dm", (ftnlen) 2);
  qnbook_(&c__3, "um", (ftnlen) 2);
  qnbook_(&c__4, "dp", (ftnlen) 2);
  qnbook_(&c__5, "up", (ftnlen) 2);
  qnbook_(&c__6, "sp", (ftnlen) 2);
  qnbook_(&c__7, "cp", (ftnlen) 2);
  qnbook_(&c__8, "bp", (ftnlen) 2);
  qnlset_("W1ANA", &c_true, (ftnlen) 5);
  qnlset_("W2NUM", &c_true, (ftnlen) 5);
  qnlset_("W2STF", &c_false,(ftnlen) 5);
  
  if(strcmp(gridname, "none") == 0) qnfilw_(&c__0, &c__0);
  else opengridfile_(gridname, (ftnlen) strlen(gridname));
}


static void qcdnum_initpdf(const lhpdf_pdf_t *pdf)
{
  double tc2, tb2, as0, q2, q2fit = (pdf->pdfset)->q2fit;
  int iq0, iqb, iqc, nf0, ix, c__1 = 1;
  int nq = (pdf->pdfset)->nq;
  double q0;
  
  tc2 = ((pdf->pdfset)->threshold[10])*((pdf->pdfset)->threshold[10]);
  tb2 = ((pdf->pdfset)->threshold[11])*((pdf->pdfset)->threshold[11]);

  iq0 = iqfromq_(&q2fit);
  iqc = iqfromq_(&tc2);
  iqb = iqfromq_(&tb2);
  nf0 = nflget_(&iq0);

  q0 = lhpdf_pdf_scale(pdf);
  as0 = lhpdf_pdf_alfas(pdf);
  q2 = q0*q0;
  
  qnrset_("ALFAS", &as0, (ftnlen) 5);
  qnrset_("ALFQ0", &q2, (ftnlen) 5);
  
  for(ix = 1; ix <= (int) ((pdf->pdfset)->nx); ix++) {
    int i;
    double f[13];
    double __nf0 = (double) nf0;
    double singlet = 0.0, gl, dm, um, dp, up, sp;
    
    (*((pdf->pdfset)->init_pdf_func))(pdf->fparm, (double) xfromix_(&ix), f+6);
    
    for(i = 1; i <= (int) nf0; i++)  
      singlet += f[i+6] + f[6-i];
    
    gl = f[6];
    dm = f[7] - f[5];
    um = f[8] - f[4];
    dp = f[7] + f[5] - singlet/__nf0;
    up = f[8] + f[4] - singlet/__nf0;
    sp = f[9] + f[3] - singlet/__nf0;
    
    qnpset_("SINGLET", &ix, &iq0, &singlet, (ftnlen) 7);
    qnpset_("GLUON",   &ix, &iq0, &gl,      (ftnlen) 5);
    qnpset_("DM",      &ix, &iq0, &dm,      (ftnlen) 2);
    qnpset_("UM",      &ix, &iq0, &um,      (ftnlen) 2);
    qnpset_("DP",      &ix, &iq0, &dp,      (ftnlen) 2);
    qnpset_("UP",      &ix, &iq0, &up,      (ftnlen) 2);
    qnpset_("SP",      &ix, &iq0, &sp,      (ftnlen) 2);
  }
  
  evolsg_(&iq0, &c__1, &nq);
  evolnm_("DM", &iq0, &c__1, &nq, (ftnlen) 2);
  evolnm_("UM", &iq0, &c__1, &nq, (ftnlen) 2);
  evlsea_("dp", &iq0, &iqc, &iqb, &nq, (ftnlen) 2);
  evlsea_("up", &iq0, &iqc, &iqb, &nq, (ftnlen) 2);
  evlsea_("sp", &iq0, &iqc, &iqb, &nq, (ftnlen) 2);

  /* --   Heavy quark evolution */
  evolcp_("cp", &iqc, &iqb, &nq, (ftnlen) 2);
  evolbp_("bp", &iqb, &nq, (ftnlen) 2);
}


/*   structure for storing the mutable variables  */
typedef struct {
  double q, x;
  int iq, ix, ng;
} qcdnum_last_call_t;


typedef struct {
  /*   x and q2 tables and related parameters  */
  int ngrver;
  unsigned int nxx, nq2, *nfmap;
  double *xxtab, *q2tab, xmicut, qmicut, qmacut, rs2cut, qminas;
  
  /*   The weight and pdf tables  */
  double pwgt[11][9][3], *pdfqcd[11];
  
  /* mutable variables */
  qcdnum_last_call_t *last;
} qcdnum_grid_t;


static int iqfromq(const qcdnum_grid_t *qcgrid, double q)
{
  static const double epsi = 1e-6;
  int i, iq;
  
  if((q == (qcgrid->last)->q) && (qcgrid->ngrver == (qcgrid->last)->ng))
    return (qcgrid->last)->iq;
  
  iq = (qcgrid->last)->iq = 0;
  (qcgrid->last)->ng = qcgrid->ngrver;
  (qcgrid->last)->q = q;
  
  if((qcgrid->nq2 == 0) || (q/qcgrid->q2tab[0] < 1.0 - epsi) ||
     (q/qcgrid->q2tab[qcgrid->nq2 - 1] > epsi + 1.0)) return 0;
  
  for(i = 0; i < qcgrid->nq2; i++)
    if(qcgrid->q2tab[i]/q <= epsi + 1.0) iq = i+1;
  
  (qcgrid->last)->iq = fabs(qcgrid->q2tab[iq-1]/q-1.0) <= epsi ? iq : -iq;
  
  return (qcgrid->last)->iq;
}


static int ixfromx(const qcdnum_grid_t *qcgrid, double x)
{
  static const double epsi = 1e-6;
  int i, ix;
  
  if((x == (qcgrid->last)->x) && (qcgrid->ngrver == (qcgrid->last)->ng))
    return (qcgrid->last)->ix;
 
  ix = (qcgrid->last)->ix = 0;
  (qcgrid->last)->ng = qcgrid->ngrver;
  (qcgrid->last)->x = x;
  
  if((x > 1.0) || (qcgrid->nxx <= 0) || (x/qcgrid->xxtab[0] < 1.0-epsi))
    return 0;
  
  for(i = 0; i < qcgrid->nxx; ++i)
    if(qcgrid->xxtab[i]/x <= epsi + 1.0) ix = i+1;
  
  (qcgrid->last)->ix = abs(qcgrid->xxtab[ix-1]/x-1.0) <= epsi ? ix : -ix;
  
  return (qcgrid->last)->ix;
}

static int icutxq(const qcdnum_grid_t *qcgrid, double x, double q)
{
  int i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0;
  
  if((qcgrid->nxx <= 0) || (qcgrid->nq2 <= 0) || (x > 1.0)) 
    return 11111;
  
  if(x < qcgrid->xxtab[0] || (x < qcgrid->xmicut && qcgrid->xmicut > 0.0)) i1 = 1;
  if(q < qcgrid->q2tab[0] || (q < qcgrid->qmicut && qcgrid->qmicut > 0.0)) i2 = 1;
  if(q > qcgrid->q2tab[qcgrid->nq2-1] || (q > qcgrid->qmacut && qcgrid->qmacut > 0.0)) i3 = 1;
  if(q > x * qcgrid->rs2cut && qcgrid->rs2cut > 0.0) i4 = 1;
  if(q < qcgrid->q2tab[0] || (q < qcgrid->qminas && qcgrid->qminas > 0.0)) i5 = 1;
  
  return i5*10000 + i4*1000 + i3*100 + i2*10 + i1;
}


static 
double get_pdfij(const qcdnum_grid_t *qcgrid, unsigned int id, unsigned int ix, unsigned int iq)
{
  double ret_val = 0.0;
  unsigned int i, nf = qcgrid->nfmap[iq-1];

  for(i = 0; i < 11; i++) {
    ret_val += qcgrid->pwgt[i][id][nf-3] 
      * qcgrid->pdfqcd[i][(ix-1)*(qcgrid->nq2) + iq-1];
  }
  
  return ret_val;
} 

static double get_pdfxq(const qcdnum_grid_t *qcgrid, int id, int ix, int iq, double tx, double tq)
{
  double f1, f2, f11, f12, f21, f22;
  
  f11 = get_pdfij(qcgrid, id, ix, iq);
  f12 = get_pdfij(qcgrid, id, ix, iq+1);
  f21 = get_pdfij(qcgrid, id, ix+1, iq);
  f22 = get_pdfij(qcgrid, id, ix+1, iq+1);
  f1 = (1.0-tq)*f11 + tq*f12;
  f2 = (1.0-tq)*f21 + tq*f22;
  
  return (1.0-tx)*f1 + tx*f2;
}


double qpdfxq(const qcdnum_grid_t *qcgrid, int id, double x, double q)
{
  int jfl, ix, iq, nf;
  double tx, tq;
  
  jfl = icutxq(qcgrid, x, q);
  if(jfl != 0) return 0.0;
  
  ix = abs(ixfromx(qcgrid, x));
  iq = abs(iqfromq(qcgrid, q));
  if(iq > qcgrid->nq2 - 1) iq = qcgrid->nq2 - 1;
  nf = qcgrid->nfmap[iq-1];

  tx = (x-qcgrid->xxtab[ix-1])/(qcgrid->xxtab[ix]-qcgrid->xxtab[ix-1]);
  tq = log(q/qcgrid->q2tab[iq-1])/log(qcgrid->q2tab[iq]/qcgrid->q2tab[iq-1]);
  return get_pdfxq(qcgrid, id, ix, iq, tx, tq);
}


static
void qcdnum_evolvepdf(const lhpdf_pdf_t *pdf, double x, double q, double *f)
{
  qcdnum_grid_t *grid = (qcdnum_grid_t *) (pdf->grid);
  double tc2 = ((pdf->pdfset)->threshold[10])*((pdf->pdfset)->threshold[10]);
  double tb2 = ((pdf->pdfset)->threshold[11])*((pdf->pdfset)->threshold[11]);
  double q2 = q*q, singlet, dm, um, dp, up, sp, cp, bp;
  unsigned int nf = 5U;
  
  if(q2 < tb2) nf = 4;
  if(q2 < tc2) nf = 3;
  
  f[0] = qpdfxq(grid, 0, x, q2);
  singlet = qpdfxq(grid, 1, x, q2);
  dm = qpdfxq(grid, 2, x, q2);
  um = qpdfxq(grid, 3, x, q2);
  dp = qpdfxq(grid, 4, x, q2);
  up = qpdfxq(grid, 5, x, q2);
  sp = qpdfxq(grid, 6, x, q2);
  
  f[-2] = 0.5*(up - um + singlet/nf);
  if(f[-2] < 0.0) f[-2] = 0.0;

  f[-1] = 0.5*(dp - dm + singlet/nf);
  if(f[-1] < 0.0) f[-1] = 0.0;
  
  
  f[3] = f[-3] = 0.5*(sp + singlet/nf);
  if(f[3] < 0.0) f[3] = f[-3] = 0.0;
  
  if (nf >= 4) {
    cp = qpdfxq(grid, 7, x, q2);
    f[4] = f[-4] = 0.5*(cp + singlet/nf);
    if(f[4] < 0.0) f[4] = f[-4] = 0.0;
  } else f[4] = f[-4] = 0.0;

  if (nf >= 5) {
    bp = qpdfxq(grid, 8, x, q2);
    f[5] = f[-5] = 0.5*(bp + singlet/nf);
    if(f[5] < 0.0) f[5] = f[-5] = 0.0;
  } else f[5] = f[-5] = 0.0;

  f[1] = dm + f[-1];
  if(f[1] < 0.0) f[1] = 0.0;
  
  f[2] = um + f[-2];
  if(f[2] < 0.0) f[2] = 0.0;
  
  f[6] = f[-6] = 0.0;
}

extern int getngrver_();
extern int getnxx_();
extern int getnq2_();
extern int getnfmap_(int *);
extern double getxxtab_(int *);
extern double getq2tab_(int *);
extern double getpdfqcd_(int *, int *, int *);
extern double getpwgt_(int *, int *, int *);
extern double getxmicut_();
extern double getqmicut_();
extern double getqmacut_();
extern double getqminas_();
extern double getrs2cut_();

static void qcdnum_grid_alloc(lhpdf_pdf_t *pdf)
{
  int i, ix, iq, id, nf;
  unsigned int nxq;
  qcdnum_grid_t *grid;
  
  /*  Initialize the the evolution code */
  qcdnum_initevolvecode(pdf->pdfset);
  
  /*  Do the evolution and build the grid */
  qcdnum_initpdf(pdf);
  
  /*  Allocate memory for the grid */
  grid = (qcdnum_grid_t *) (pdf->grid = malloc(sizeof(qcdnum_grid_t)));
  grid->last = (qcdnum_last_call_t *) malloc(sizeof(qcdnum_last_call_t));

  /*  Grid : Mutable variables  */
  (grid->last)->q = (grid->last)->x = 0.0;
  (grid->last)->iq = (grid->last)->ix = (grid->last)->ng = 0;
  
  /*  Get the size and some parameter of the grid  */
  grid->ngrver = (unsigned int) getngrver_();
  grid->nxx = (unsigned int) getnxx_();
  grid->nq2 = (unsigned int) getnq2_();
  grid->xmicut = (double) getxmicut_();
  grid->qmicut = (double) getqmicut_();
  grid->qmacut = (double) getqmacut_();
  grid->rs2cut = (double) getrs2cut_();
  grid->qminas = (double) getqminas_();

  
  /*  Get the x, q2 and nfmap arrays  */
  grid->xxtab = (double *) malloc((grid->nxx + 1)*sizeof(double));
  for(i = 1; i <= grid->nxx+1; i++) grid->xxtab[i-1] = getxxtab_(&i);
  
  grid->q2tab = (double *) malloc((grid->nq2)*sizeof(double));
  grid->nfmap = (unsigned int *) malloc((grid->nq2)*sizeof(unsigned int));
  for(i = 1; i <= grid->nq2; i++) {
    grid->q2tab[i-1] = getq2tab_(&i);
    grid->nfmap[i-1] = getnfmap_(&i);
  }
  
  /*  Grid : The grid  */
  nxq = (grid->nxx + 1)*(grid->nq2);
  for(i = 0; i < 11; i++) {
    grid->pdfqcd[i] = malloc(nxq*sizeof(double));
    
    for(ix = 1; ix <= grid->nxx+1; ix++) 
      for(iq = 1; iq <= grid->nq2; iq++)
	grid->pdfqcd[i][(ix-1)*(grid->nq2) + iq-1] = getpdfqcd_(&i,&ix,&iq);
  }
  
  /*  Get the weihghts  */
  for(i = 0; i < 11; i++)
    for(id = 0; id < 9; id++) 
      for(nf = 3; nf <= 5; nf++)
	grid->pwgt[i][id][nf-3] = getpwgt_(&i, &id, &nf);  
}

static void qcdnum_grid_free(void *s) 
{
  unsigned int i;
  qcdnum_grid_t *grid = (qcdnum_grid_t *) s;
  
  if(grid->last)  free(grid->last);
  if(grid->xxtab) free(grid->xxtab);
  if(grid->q2tab) free(grid->q2tab);
  if(grid->nfmap) free(grid->nfmap);
  
  for(i = 0; i < 11; i++)
    if(grid->pdfqcd[i]) 
      free(grid->pdfqcd[i]);

  if(grid) free(grid);
}


static const lhpdf_evlcode_t evlcode_qcdnum = {
  0, 
  qcdnum_evolvepdf,
  qcdnum_grid_alloc,
  qcdnum_grid_free
};

const lhpdf_evlcode_t *lhpdf_evlcode_qcdnum = &evlcode_qcdnum;

