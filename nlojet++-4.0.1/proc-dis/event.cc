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
#include "bits/dis-event.h"



namespace nlo {


  void lab_to_breit(event_dis& p)
  {
    lorentzvector<double> q(p[-1]-p[-2]);
    double xb = -0.5*q.mag2()/(p[0]*q);
    threevector<double> bVec = -((q + 2.0*xb*p[0]).boostVector());
    int i, low = p.lower(), up = p.upper();
    
    p[hadron(0)].boost(bVec);
    for(i = low; i <= up; i++) 
      p[i].boost(bVec);
    
    double phi = p[0].phi(), theta = p[0].theta();
    
    p[hadron(0)].rotateZ(-phi);
    p[hadron(0)].rotateY(-theta);  
    
    for(i = low; i <= up; i++) {
      p[i].rotateZ(-phi);
      p[i].rotateY(-theta);
    }
  }
  
  void breit_to_lab(double x, event_dis& p)
  {
    int i, low = p.lower(), up = p.upper();
    double bz = (1.0 - x)/(1.0 + x);
    threevector<double> bVec = -((p[-1] + p[hadron(0)]).boostVector());
    
    p[hadron(0)].boost(bVec);
    double phi = p[hadron(0)].phi(), theta = p[hadron(0)].theta();
    
    p[hadron(0)].rotateZ(-phi);
    p[hadron(0)].rotateY(-theta);  
    p[hadron(0)].boost(0.0, 0.0, bz);
    
    for(i = low; i <= up; i++) {
      p[i].boost(bVec);
      p[i].rotateZ(-phi);
      p[i].rotateY(-theta);
      p[i].boost(0.0, 0.0, bz);
    }
  }
}  // namespace nlo
