/* Standard C includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* lhpdf includes */
#include "lhpdf.h"

extern int bldat1_();
extern int bldat2_();
extern int bldat3_();

typedef int ftnlen;

extern double alpi_(double *);
extern int evlpar_(int *, char *, double *, int *, ftnlen);
extern int parqcd_(int *, char *, double *, int *, ftnlen);
extern int parpdf_(int *, char *, double *, int *, ftnlen);
extern int evolve_(double (*)(), int *);
extern int setlam_(int *, double *, int *);


typedef struct {
  /* max Nf */
  unsigned int nfmx;
    
  /* Lambda(nf=5,MSbar), quark masses */
  double al, valqms[6];
  
  /* xv, tv arrays and the pdf grid */
  unsigned int nx, nt;
  double *xv, *xvpow, *tv, *upd;
} evlcteq_grid_t;



static double evlcteq_alfas(double Q)
{
  /* THIS IS JUST AN INTERFACE TO THE F77 FUNCTION. 
   * IT IS PDFSET DEPENDENT!
   */
  return 3.14159265358979323848*alpi_(&Q);
}


static void evlcteq_initevolvecode(const lhpdf_pdfset_t *pdfs)
{
  int c__1 = 1, j;
  double ordi, ahdn = 1.0, mass = 180.0, anf, qini;
  double qmax = pdfs->qmax, xmin = pdfs->xmin, m4 = pdfs->threshold[10],
    m5 = pdfs->threshold[11];
  
  bldat1_();
  bldat2_();
  bldat3_();

  parpdf_(&c__1, "IHDN", &ahdn, &j, (ftnlen) 4);
  parpdf_(&c__1, "QMAX", &qmax, &j, (ftnlen) 4);
  parpdf_(&c__1, "XMIN", &xmin, &j, (ftnlen) 4);
  
  ordi = pdfs->alford + 1.0;
  parqcd_(&c__1, "ORDR", &ordi, &j, (ftnlen) 4); 
  parqcd_(&c__1, "M4", &m4, &j, (ftnlen) 2);
  parqcd_(&c__1, "M5", &m5, &j, (ftnlen) 2);
  parqcd_(&c__1, "M6", &mass, &j, (ftnlen) 2);
  
  anf  = (double) lhpdf_pdfset_nfmax(pdfs);
  qini = sqrt(pdfs->q2fit);
  ordi = pdfs->evlord + 1.0;

  parqcd_(&c__1, "NFL",  &anf,  &j, (ftnlen) 3);
  parpdf_(&c__1, "QINI", &qini, &j, (ftnlen) 4);
  evlpar_(&c__1, "NFMX", &anf,  &j, (ftnlen) 4);
  parpdf_(&c__1, "IKNL", &ordi, &j, (ftnlen) 4);
}

static double astolam(double as, double q, unsigned int nloop, unsigned int nf)
{
  static const double pi = 3.141592653589793;
  double b, t, bp, ot, as0, as1, xlp, xlt;
  
  xlp = nloop - 1.0;
  b = (33.0 - 2.0*nf)/(12.0*pi);
  bp = (153.0 - 19.0*nf)/(2.0*pi*(33.0 - 2.0*nf))*xlp;
  t = 1.0/(b*as);

  /* Solve the equation */  
  do {
    xlt = log(t);
    ot = t;
    
    /* Solve the equation */
    /* Value and Derivative of alfa with respect to t */
    as0 = (1 - bp*xlt/(b*t))/(b*t);
    as1 = -1.0/(b*t*t) - bp/(b*b)*(1.0 - 2.0*xlt)/(t*t*t);
    t += (as - as0)/as1;
  } while(fabs(ot-t)/ot > 1e-5);
  
  return q/exp(t/2.0);
}

static const lhpdf_pdf_t *evlcteq_current_pdf;

static double evlcteq_parmflavor(int *i, double *x)
{
  int i0;
  double f[13];
  
  (*((evlcteq_current_pdf->pdfset)->init_pdf_func))(evlcteq_current_pdf->fparm, *x, f+6);
  switch(*i) {
  case -2: i0 =-1; break;
  case -1: i0 =-2; break;
  case  1: i0 = 2; break;
  case  2: i0 = 1; break;
  default: i0 = *i;
  }
  
  return f[i0+6]/(*x);
}

static void evlcteq_initpdf(const lhpdf_pdf_t *pdf)
{
  int c__1 = 1, j, i__1, nflavor;
  double alfas0, q0, blam, anx, anq;

  /* Initialize the evolution code */
  evlcteq_initevolvecode(pdf->pdfset);
  
  /* Calculate the grid */
  j = (pdf->pdfset)->alford; i__1 = j + 1;
  q0 = lhpdf_pdf_scale(pdf);
  alfas0 = lhpdf_pdf_alfas(pdf);
   
  nflavor = (int) lhpdf_pdfset_nfmax(pdf->pdfset);
  blam = astolam(alfas0, q0, i__1, nflavor);

  anx = (pdf->pdfset)->nx;
  anq = (pdf->pdfset)->nq - 1;
  
  setlam_(&nflavor, &blam, &i__1);
  parpdf_(&c__1, "NX", &anx, &j, (ftnlen) 2);
  parpdf_(&c__1, "NT", &anq, &j, (ftnlen) 2);

  evlcteq_current_pdf = pdf;
  evolve_((double (*)()) evlcteq_parmflavor, &j);
  
  if(j != 0) {
    fprintf(stderr, "EVLCTEQ Evolve Error code : %d", j);
    exit(-1);
  }
}

static
double evlcteq_polint(const double *xa, const double *ya, int n, double x)
{
  int m, i, ns = 1;
  static double __c[10], __d[10];
  double *c = __c - 1, *d = __d - 1;
  double den, ho, hp, w, y, dift, dif;
  
  //   adjust parameters
  --xa; --ya;      
  
  dif = fabs(x-xa[1]);
  for(i = 1; i <= n; i++) {
    if((dift = fabs(x-xa[i])) < dif) { 
      ns = i; dif = dift;
    }
    c[i] = d[i] = ya[i];
  }
  
  y = ya[ns--];
  for(m = 1; m < n; m++) {
    for(i = 1; i <= n-m; i++) {
      ho = xa[i] - x;
      hp = xa[i+m] - x; 
      w = c[i+1] - d[i];
      if((den = ho-hp) == 0.0) exit(-1);
      den = w/den; 
      d[i] = hp*den; 
      c[i] = ho*den;
    }
    y += 2*ns < n-m ? c[ns+1] : d[ns--];
  }
  
  return y;
}


static 
double evlcteq_pardis(const void *grid, int iprtn, double x, double q)
{
  int jm, jx, jlx = -1, ju, jq, jlq = -1, j1, ip, jtmp, it;
  double tt, ss, fij[4], fx, ff, fvec[4];
  evlcteq_grid_t *evl = (evlcteq_grid_t *) grid;
  unsigned int nx = evl->nx, nq = evl->nt;

  if(iprtn != 0)
    if(q <= evl->valqms[abs(iprtn)-1]) 
      return 0.0;

  ju = nx + 1;
  while(ju - jlx > 1) {
    jm = (ju + jlx)/2;
    if(x >= evl->xv[jm]) jlx = jm;
    else ju = jm;
  }
  
  if(jlx <= -1) {
    fprintf(stderr, "Severe error: x <= 0 in pardis! x = %g \n", x);
    exit(-1);
  } else if(jlx == 0) jx = 0;
  else if(jlx <= (int) nx-2) jx = jlx-1;
  else if(jlx == (int) nx-1 || x < 1.00001) jx = jlx-2;
  else {
    fprintf(stderr, "Severe error: x > 1 in pardis! x = %g", x);
    exit(-1);
  }

  tt = log(log(q/evl->al));
  ss = pow(x, 0.3);
    
  ju = nq + 1;
  while(ju - jlq > 1) {
    jm = (ju + jlq)/2;
    if (tt >= evl->tv[jm]) jlq = jm;
    else ju = jm;
  }
    
  if(jlq <= 0) jq = 0;
  else if(jlq <= (int) nq-2) jq = jlq-1;
  else jq = nq-3;
    
  ip = (iprtn >= 3 ? -iprtn : iprtn);
  jtmp = ((ip + evl->nfmx)*(nq+1) + (jq-1))*(nx+1) + jx + 1;
    
  for(it = 0; it < 4; it++) {
    j1 = jtmp + (it+1)*(nx+1);
    if(jx == 0) {
      fij[0] = 0.0;
      fij[1] = (evl->upd[j1  ])*(evl->xv[1])*(evl->xv[1]);
      fij[2] = (evl->upd[j1+1])*(evl->xv[2])*(evl->xv[2]);
      fij[3] = (evl->upd[j1+2])*(evl->xv[3])*(evl->xv[3]);

      fx = evlcteq_polint(evl->xvpow, fij, 4, ss);
      if (x > 0.0) fvec[it] = fx/(x*x);
    } else if(jlx == (int) nx - 1)
      fvec[it] = evlcteq_polint(evl->xvpow + nx-3, evl->upd + j1-1, 4, ss);
    else {
      double *svec = evl->xvpow + jx-1;
      double s12 = svec[1]-svec[2], s13 = svec[1]-svec[3],
	s23 = svec[2]-svec[3], s24 = svec[2]-svec[4],
	s34 = svec[3]-svec[4], sy2 = ss-svec[2], sy3 = ss-svec[3];
	
      double s1213 = s12+s13, s2434 = s24+s34;	
      double sf2 = evl->upd[j1], sf3 = evl->upd[j1+1];
	
      fvec[it] = (sf2*sy3-sf3*sy2)/s23
	+ ((s34*sy2-s2434*sy3)*(evl->upd[j1-1] - sf2*s13/s23+sf3*s12/s23)/s12 
	   - (s12*sy3-s1213*sy2)*(evl->upd[j1+2] + sf2*s34/s23-sf3*s24/s23)/s34
	   )*sy2*sy3/(s23*(s12*s34 - s1213*s2434));
    }
  }
    
  if(jlq <= 0) ff = evlcteq_polint(evl->tv, fvec, 4, tt);
  else if(jlq >= (int) nq-1)
    ff = evlcteq_polint(evl->tv + nq-3, fvec, 4, tt);
  else {
    double *tvec = evl->tv + jq-1;
    double t12 = tvec[1]-tvec[2], t13 = tvec[1]-tvec[3], tf3 = fvec[2],
      t23 = tvec[2]-tvec[3], t24 = tvec[2]-tvec[4], tf2 = fvec[1], h00,
      t34 = tvec[3]-tvec[4], ty2 = tt-tvec[2], ty3 = tt-tvec[3];
      
    h00 = (t34*ty2-(t24+t34)*ty3)*(fvec[0]-(tf2*t13-tf3*t12)/t23)/t12 
      -   (t12*ty3-(t12+t13)*ty2)*(fvec[3]-(tf3*t24-tf2*t34)/t23)/t34;
      
    ff = (h00*ty2*ty3/(t12*t34-(t12+t13)*(t24+t34))+tf2*ty3-tf3*ty2)/t23;
  }
    
  return ff < 0.0 ? 0.0 : ff;
}
  
extern int getnx_(void);
extern int getxv_(double *, double *);
extern int getnt_(void);
extern int gettv_(double *);
extern double getal_(void);
extern int getvalqms_(int *, double *);
extern int getupd_(double *);


static void evlcteq_grid_alloc(lhpdf_pdf_t *pdf)
{
  int nfmx, nupd;
  evlcteq_grid_t *evl;

  evlcteq_initpdf(pdf);
  pdf->grid = malloc(sizeof(evlcteq_grid_t));
  evl = (evlcteq_grid_t *) (pdf->grid);
  
  evl->nx = getnx_();
  evl->xv = (double *) malloc((evl->nx + 1)*sizeof(double));
  evl->xvpow = (double *) malloc((evl->nx +1)*sizeof(double));
  getxv_(evl->xv, evl->xvpow);
  
  evl->nt = getnt_();
  evl->tv = (double *) malloc((evl->nt + 1)*sizeof(double));
  gettv_(evl->tv);
  
  getvalqms_(&nfmx, evl->valqms);
  evl->nfmx = nfmx;
  evl->al = getal_();
  
  nupd = (evl->nx+1)*(evl->nt+1)*(evl->nfmx+3);
  evl->upd = (double *) malloc(nupd*sizeof(double));
  getupd_(evl->upd);
}

static void evlcteq_grid_free(void *s) 
{
  evlcteq_grid_t *evl = (evlcteq_grid_t *) s;
  
  if(evl == 0) return;
  
  if(evl->xv) free(evl->xv);
  if(evl->tv) free(evl->tv);
  if(evl->xvpow) free(evl->xvpow);
  if(evl->upd) free(evl->upd);
  
  free(evl);
}

static
void evlcteq_evolvepdf(const lhpdf_pdf_t *pdf, double x, double q, double *f)
{
  unsigned int nf = lhpdf_pdfset_nfmax(pdf->pdfset);

  f[ 2] = x*evlcteq_pardis(pdf->grid,  1, x, q);
  f[ 1] = x*evlcteq_pardis(pdf->grid,  2, x, q);
  f[-2] = x*evlcteq_pardis(pdf->grid, -1, x, q);
  f[-1] = x*evlcteq_pardis(pdf->grid, -2, x, q);
  f[ 0] = x*evlcteq_pardis(pdf->grid,  0, x, q);
  
  if(nf >= 3)
    f[3] = f[-3] = x*evlcteq_pardis(pdf->grid, -3, x, q);
  else f[3] = f[-3] = 0.0;
  
  if(nf >= 4)
    f[4] = f[-4] = x*evlcteq_pardis(pdf->grid, -4, x, q);
  else f[4] = f[-4] = 0.0;
  
  if(nf >= 5)
    f[5] = f[-5] = x*evlcteq_pardis(pdf->grid, -5, x, q);
  else f[5] = f[-5] = 0.0;

  
  if(nf == 6)
    f[6] = f[-6] = x*evlcteq_pardis(pdf->grid, -6, x, q);
  else f[6] = f[-6] = 0.0;
}

static const lhpdf_evlcode_t evlcode_evlcteq = {
  evlcteq_alfas,
  evlcteq_evolvepdf,
  evlcteq_grid_alloc,
  evlcteq_grid_free
};


const lhpdf_evlcode_t *lhpdf_evlcode_evlcteq = &evlcode_evlcteq;

