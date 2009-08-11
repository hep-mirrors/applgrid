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

//   nlo includes
#include "bits/nlo-flkern.h"



namespace nlo {
  
#define pi56 8.22467033424113218234 


  void Kqq(double x, double al, double *K) 
  {
    double ox = 1.0 - x;
    double la = log(al), lx = log(ox/x);
    double ax = (al < ox ? ox/(ox + al) : (1.0 + ox)/(ox + al));
    
    K[0] = Sqq(x)*(la + lx) + Cf*(ox + 2.0*log(ax)/ox);
    K[1] = -2.0*Cf*lx/ox; K[2] = Cf*(la*la + pi56) - Gq - Kq(al);
  }
  
  void Kgg(double x, unsigned int nf, double al, double *K) 
  {
    double ox = 1.0 - x;
    double la = log(al), lx = log(ox/x);
    double ax = (al < ox ? ox/(ox + al) : (1.0 + ox)/(ox + al));
    
    K[0] = Sgg(x)*(la + lx) + 2.0*Ca*log(ax)/ox;
    K[1] = -2.0*Ca*lx/ox; K[2] = Ca*(la*la + pi56) - Gg(nf) - Kg(nf, al);
  }


  void Kqq(double x, double jac, double al, double *K) 
  {
    double ox = 1.0 - x;
    double la = log(al), lx = log(ox/x);
    double ax = (al < ox ? ox/(ox + al) : (1.0 + ox)/(ox + al));
    
    K[0] = jac*(Sqq(x)*(la + lx) + Cf*(ox + 2.0*log(ax)/ox))/x;
    K[1] = -2.0*jac*Cf*lx/ox + Cf*(la*la + pi56) - Gq - Kq(al);
  }
  
  void Kgg(double x, double jac, unsigned int nf, double al, double *K) 
  {
    double ox = 1.0 - x;
    double la = log(al), lx = log(ox/x);
    double ax = (al < ox ? ox/(ox + al) : (1.0 + ox)/(ox + al));
    
    K[0] = jac*(Sgg(x)*(la + lx) + 2.0*Ca*log(ax)/ox)/x;
    K[1] = -2.0*jac*Ca*lx/ox + Ca*(la*la + pi56) - Gg(nf) - Kg(nf, al);
  }


  void Tqq(double x, double al, double *K)
  {
    double ox = 1.0-x;
    double la = log(al);

    if(al > ox) K[0] = Sqq(x)*log(ox/al)+2.0*Cf*log((ox+al)/(1.0+ox))/ox;
    else K[0] = 2.0*Cf*log((ox+al)/ox)/ox;

    K[1] = -2.0*Cf*log(ox)/ox;
    K[2] = -Cf*(3.28986813369645287298+la*la);
  }

  void Tgg(double x, double al, double *K)
  {
    double ox = 1.0-x;
    double la = log(al);

    if(al > ox) K[0] = Sgg(x)*log(ox/al)+2.0*Ca*log((ox+al)/(1.0+ox))/ox;
    else K[0] = 2.0*Ca*log((ox+al)/ox)/ox;
    
    K[1] = -2.0*Ca*log(ox)/ox;
    K[2] = -Ca*(3.28986813369645287298+la*la);
  }

  void Tqq(double x, double jac, double al, double *K)
  {
    double ox = 1.0-x;
    double la = log(al);

    if(al > ox) K[0] = Sqq(x)*log(ox/al)+2.0*Cf*log((ox+al)/(1.0+ox))/ox;
    else K[0] = 2.0*Cf*log((ox+al)/ox)/ox;
    K[0] *= jac/x;
    K[1] = -2.0*jac*Cf*log(ox)/ox - Cf*(3.28986813369645287298+la*la);
  }
  
  void Tgg(double x, double jac, double al, double *K)
  {
    double ox = 1.0-x;
    double la = log(al);

    if(al > ox) K[0] = Sgg(x)*log(ox/al)+2.0*Ca*log((ox+al)/(1.0+ox))/ox;
    else K[0] = 2.0*Ca*log((ox+al)/ox)/ox;
    K[0] *= jac/x;   
    K[1] = -2.0*jac*Ca*log(ox)/ox - Ca*(3.28986813369645287298+la*la);
  }
}   //  namespace nlo

