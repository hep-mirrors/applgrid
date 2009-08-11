//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#ifndef __NLO_NLO_FLKERN_H__
#define __NLO_NLO_FLKERN_H__ 1


#include <bits/nlo-color.h>



namespace nlo {

  //   Alterelli-Parisi splitting functions
  inline double Sqg(double x) {
    return Cf*(1.0 + (1.0 - x)*(1.0 - x))/x;
  }
  
  inline double Sgq(double x) {
    return 0.5 - x*(1.0 - x);
  }  
  
  inline double Sqq(double x) {
    return Cf*(1.0 + x*x)/(1.0 - x);
  }
  
  inline double Sgg(double x) {
    double omx = 1.0 - x;
    return 2.0*Ca*(x/omx + omx/x + x*omx);
  }
  
  //   Alterelli-Parisi probabilities
  inline void Pqg(double x, double *P) {
    P[0] = Sqg(x); P[1] = P[2] = 0.0;
  }
  
  inline void Pgq(double x, double *P) {
    P[0] = Sgq(x); P[1] = P[2] = 0.0;
  }  
  
  inline void Pqq(double x, double *P) {
    P[0] = Sqq(x); P[1] = -2.0*Cf/(1.0 - x); P[2] = Gq;
  }
  
  inline void Pgg(double x, unsigned int nf, double *P) {
    P[0] = Sgg(x); P[1] = -2.0*Ca/(1.0 - x); P[2] = Gg(nf); 
  }


  inline void Pqg(double x, double jac, double *P) {
    P[0] = jac*Sqg(x)/x; P[1] = 0.0;
  }

  inline void Pgq(double x, double jac, double *P) {
    P[0] = jac*Sgq(x)/x; P[1] = 0.0;
  }  
  
  inline void Pqq(double x, double jac, double *P) {
    P[0] = jac*Sqq(x)/x; P[1] = -2.0*jac*Cf/(1.0 - x) + Gq;
  }
  
  inline void Pgg(double x, double jac, unsigned int nf, double *P) {
    P[0] = jac*Sgg(x)/x; P[1] = -2.0*jac*Ca/(1.0 - x) + Gg(nf); 
  }

  //    Some other flavour kernel which help to implement 
  //    the Catani-Seymour substractin method for NLO calculations 
  inline void Kqg(double x, double al, double *K) {
    K[0] = Sqg(x)*std::log(al*(1.0 - x)/x) + Cf*x; K[1] = K[2] = 0.0;
  }
  
  inline void Kgq(double x, double al, double *K) {
    K[0] = Sgq(x)*std::log(al*(1.0 - x)/x) + x*(1.0 - x); K[1] = K[2] = 0.0;
  }

  void Kqq(double, double, double *); 
  void Kgg(double, unsigned int, double, double *); 


  inline void Kqg(double x, double jac, double al, double *K) {
    K[0] = jac*(Sqg(x)*std::log(al*(1.0 - x)/x) + Cf*x)/x; K[1] = 0.0;
  }
  
  inline void Kgq(double x, double jac, double al, double *K) {
    K[0] = jac*(Sgq(x)*std::log(al*(1.0 - x)/x) + x*(1.0 - x))/x; K[1] = 0.0;
  }

  void Kqq(double, double, double, double *); 
  void Kgg(double, double, unsigned int, double, double *); 


  inline void Tqg(double x, double al, double *K) 
  {
    K[0] = K[1] = K[2] = 0.0;
    if(x > 1.0 - al) K[0] = Sqg(x)*std::log((1.0-x)/al);
  }
  
  inline void Tgq(double x, double al, double *K)
  {
    K[0] = K[1] = K[2] = 0.0;
    if(x > 1.0 - al) K[0] = Sgq(x)*std::log((1.0-x)/al);
  }
  
  void Tqq(double, double, double *); 
  void Tgg(double, double, double *); 


  inline void Tqg(double x, double jac, double al, double *K) 
  {
    K[0] = K[1] = 0.0;
    if(x > 1.0 - al) K[0] = jac*Sqg(x)*std::log((1.0-x)/al)/x;
  }
  
  inline void Tgq(double x, double jac, double al, double *K)
  {
    K[0] = K[1] = 0.0;
    if(x > 1.0 - al) K[0] = jac*Sgq(x)*std::log((1.0-x)/al)/x;
  }
  
  void Tqq(double, double, double, double *); 
  void Tgg(double, double, double, double *); 
}   //   namespace nlo   

#endif
