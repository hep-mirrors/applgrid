
//   genpdf.cxx        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: genpdf.cxx, v0   Thu Mar  6 15:16:43 CET 2008 sutt


#include <TFile.h>
#include <TVectorT.h>
#include <TMatrixT.h>


double* make_ckmsum() { 
  double* _ckm = new double[13];
  TFile f("ckm.root");
  TVectorT<double>* _ckmsum = (TVectorT<double>*)f.Get("ckmsum");
  for ( int i=0 ; i<13 ; i++ ) { 
    _ckm[i] = (*_ckmsum)(i);
    // _ckm[i] = 1;
  }
  f.Close();
  return _ckm;
}

double** make_ckm() { 
  double** _ckm = new double*[13];
  TFile f("ckm.root");
  TMatrixT<double>* _ckm2 = (TMatrixT<double>*)f.Get("ckm2");
  for ( int i=0 ; i<13 ; i++ ) { 
    _ckm[i] = new double[13];
    for ( int j=0 ; j<13 ; j++ )  _ckm[i][j] = (*_ckm2)(i,j);
    //    _ckm[i][i] = 1;
  }
  f.Close();
  return _ckm;
}

double*  ckmsum  = make_ckmsum();
double** ckm2    = make_ckm();


 

void  mcfm_pdf__(const double* fA, const double* fB, double* H) { 

  const int nQuark = 6;
  const int iQuark = 4; 

  // MCFM W production
  //
  //  double* pdf::Hpdf_W_only(const double & x1, const double & x2, const double & fScale) 
  
  //   doDGLAP(x1,x2,fScale);

  // offset so we can use fA[-6]
  // fA += 6;
  // fB += 6;

  
  double GA=fA[6];
  double GB=fB[6];
  double QA=0; double QB=0; double QbA=0; double QbB=0;
  
  for(int i = 1; i <= iQuark; i++) 
    {
      QA += fA[nQuark + i]*ckmsum[nQuark + i];
      QB += fB[nQuark + i]*ckmsum[nQuark + i];
    }
  for(int i = -iQuark; i < 0; i++) 
    {
      QbA += fA[nQuark + i]*ckmsum[nQuark + i];
      QbB += fB[nQuark + i]*ckmsum[nQuark + i];
    }
  
  H[2]=QbA * GB;
  H[3]= QA * GB;
  H[4]= GA * QbB;
  H[5]= GA * QB;
  
  for (int i1 = 3; i1 <= 5; i1 += 2)
    {
      for(int i2 = 8; i2 <= 10; i2 += 2)
	{
	  H[0] += fA[i1]*fB[i2]*ckm2[i1][i2];
	  H[1] += fA[i2]*fB[i1]*ckm2[i2][i1];
	}
    }
  //    return H;     
  
}


