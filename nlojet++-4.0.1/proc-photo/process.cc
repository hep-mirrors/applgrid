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
#include "bits/photo-process.h"


namespace nlo {
  
  
  weight_photo pdf_and_coupling_photo::
  pdf(double x1, double x2, double mf2, unsigned int nu, unsigned int nd)
  {
    weight_photo retval;
    static double __f[13];
    static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
    
    //----- calculat the pdf -----
    double *f = __f+6;
    this -> hadron(x2, mf2, nu, nd, f);
    
    //----- gluon pdfs -----
    retval[0] = f[0];
    
    //---- up type quarks -----
    for(unsigned int u = 0; u < nu && u < 3; u++)
      retval[1] += f[iu[u]] + f[-iu[u]];
    
    //----- down type quarks -----
    for(unsigned int d = 0; d < nd && d < 3; d++)
      retval[2] += f[id[d]] + f[-id[d]];
    
    return retval*(this -> photon(x1));
  }  
}
