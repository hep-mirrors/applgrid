/* lhpdf includes */
#include "lhpdf.h"


static const char *mrst2001_describtion[] = {
  "MRST fit: hep-ph/0110215",
  "This set has 4 member PDF's. However, the publication only contains",
  "parametric information of 2 members.",
  "    mem=1 --> MRST2001: Best fit (alpha_S=0.119)",
  "    mem=2 --> MRST2001J: Jet data fit (alpha_S=0.212)",
  "Note that the best fit (mem=0) has been set equal to MRST2001",
  0
};

static void mrst2001_parm_pdf(const LHAPDF_FLOAT *fp, double x, double *f)
{
  double uv, dv, g1, g2, S, D;

  uv = lhpdf_x_taylor(x, fp, 2.0, 5, 1,2,3,4,5);
  dv = lhpdf_x_taylor(x, fp, 2.0, 5, 6,7,8,9,10);
  g1 = lhpdf_x_taylor(x, fp, 2.0, 5, 11,12,13,14,15);
  g2 = lhpdf_x_taylor(x, fp, 1.0, 3, 16,17,18);
  S  = lhpdf_x_taylor(x, fp, 2.0, 5, 19,20,21,22,23);
  D  = lhpdf_x_taylor(x, fp, 1.0, 5, 24,25,26,27,28);

  f[-6] = f[6] = 0.0;
  f[-5] = f[5] = 0.0;
  f[-4] = f[4] = 0.0;
  f[-3] = f[3] = 0.1*S;
  
  f[-2] = 0.2*S-0.5*D;
  f[-1] = 0.2*S+0.5*D;
  f[0] = g1-g2;
  f[1] = dv+f[-1];
  f[2] = uv+f[-2];
}

static const LHAPDF_FLOAT mrst2001_list[3][29] = { 
  { 0.119, 0.158, 0.25, 3.33, 5.61, 55.49, 0.040, 0.27, 3.88, 52.73, 30.65, 1.90, 0.09, 3.70, 1.26, -1.43, 0.21,  -0.33, 10.0, 0.222, -0.26, 7.10, 3.42, 10.30, 1.195, 1.24, 9.10, 14.05, -45.52},
  { 0.119, 0.158, 0.25, 3.33, 5.61, 55.49, 0.040, 0.27, 3.88, 52.73, 30.65, 1.90, 0.09, 3.70, 1.26, -1.43, 0.21,  -0.33, 10.0, 0.222, -0.26, 7.10, 3.42, 10.30, 1.195, 1.24, 9.10, 14.05, -45.52},
  { 0.121, 0.158, 0.25, 3.33, 5.61, 55.49, 0.040, 0.27, 3.88, 52.73, 30.65, 123.5,1.16, 4.69,-3.57,  3.41, 0.038, -0.5,  10.0, 0.222, -0.26, 7.10, 3.42, 10.30, 1.195, 1.24, 9.10, 14.05, -45.52}
};

static const lhpdf_pdfset_t mrst2001_pdf = { 
  /* Description of the pdf set */
  mrst2001_describtion, 

  /* Evolution */
  1, 1, 1.0, 1.0, "large.grid", 1e-6,1.0, 1.0,1e10, 400U,112U,

  /* Alpha QCD */
  1, 2, 0, 0, 0.119, 91.71, 1.43, 4.3, 180.0,

  /* Parametrization */
  &mrst2001_parm_pdf,
  {-1.0,4.3,1.43,0.0,0.0,0.0, -1.0, 0.0,0.0,0.0,1.43,4.3,-1.0}, -1,
  
  /* List of the parameters */
  3, 29, (const LHAPDF_FLOAT *) mrst2001_list
};


const lhpdf_pdfset_t *lhpdf_pdfset_mrst2001 = &mrst2001_pdf;
