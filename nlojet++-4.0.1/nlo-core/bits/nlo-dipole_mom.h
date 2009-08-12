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
#ifndef __NLO_NLO_DIPOLE_MOM_H__
#define __NLO_NLO_DIPOLE_MOM_H__ 1


namespace nlo {


  template<class _Event> 
  bool dipole_mom_fff(double al, const _Event& p, int i, int j, int k, _Event& q) 
  {
    typedef typename _Event::lorentzvector_type::value_type value_type;

    value_type sij = p[i]*p[j], sik = p[i]*p[k], sjk = p[j]*p[k];
    
    value_type y  = sij/(sij+sik+sjk);

    if(y > (value_type) al) return false;
    
    value_type oy = (sik + sjk)/(sij + sik + sjk);
    int up = p.upper(), low = p.lower();
        
    for(int m = low; m < up; m++) q[m] = p[m];
    if(j != up) q[j] = p[up];
    
    q[i] += p[j] - (y/oy)*p[k];
    q[(k == up ? j : k)] /= oy;

    return true;
  }

  template<class _Event> 
  bool dipole_mom_ffi(double al, const _Event& p, int i, int j, int a, _Event& q)
  {
    typedef typename _Event::lorentzvector_type::value_type value_type;

    int up = p.upper(), low = p.lower();
    value_type sij = p[i]*p[j], sia = p[i]*p[a], sja = p[j]*p[a];
    
    value_type xija = (sia + sja - sij)/(sia + sja), 
      oxija = sij/(sia + sja); 
    
    if(oxija > al) return false;
    
    for(int m = low; m < up; m++) q[m] = p[m];
    if(j != up) q[j] = p[up];
    
    q[a]  = xija*p[a];
    q[i] += p[j] - oxija*p[a];
    
    return true;
  }

  template<class _Event> bool 
  dipole_mom_iff(double al, const _Event& p, int a, int i, int k, _Event& q)
  {
    typedef typename _Event::lorentzvector_type::value_type  value_type;

    int up = p.upper(), low = p.lower();
    value_type sai = p[i]*p[a], sak = p[k]*p[a], sik = p[i]*p[k];

    value_type xika = (sai + sak - sik)/(sai + sak),
      oxika = sik/(sai + sak), ui = sai/(sai + sak);
    
    if(ui > al) return false;   
    
    for(int m = low; m < up; m++) q[m] = p[m];
    if(i != up) q[i] = p[up];

    q[a]  = xika*p[a];
    q[(k == up ? i : k)] += p[i] - oxika*p[a];
    
    return true;
  }

  template<class _Event> bool 
  dipole_mom_ifi(double al, const _Event& p, int a, int i, int b, _Event& q)
  {
    typedef typename _Event::lorentzvector_type lorentzvector_type;
    typedef typename lorentzvector_type::value_type value_type;
    
    int up = p.upper(), low = p.lower();
    value_type sai = p[i]*p[a], sab = p[a]*p[b], sib = p[i]*p[b];
    
    value_type xiab = (sab - sai - sib)/sab, vi = sai/sab;
    
    if(vi > al) return false;
    
    q[a] = xiab*p[a];
    q[b] = p[b];
    
    lorentzvector_type K(p[a]+p[b]-p[i]), Kb(q[a]+p[b]);
    lorentzvector_type Q(K+Kb);
    
    value_type Q2 = Q.mag2(), K2 = K.mag2();
    
    for(int m = low; m < up; m++)
      if(m != i && m != a && m != b) 
        q[m] = p[m] - Q*(2.0*(Q*p[m])/Q2) + Kb*(2.0*(K*p[m])/K2);
    
    if(i != up)
      q[i] = p[up] - Q*(2.0*(Q*p[up])/Q2) + Kb*(2.0*(K*p[up])/K2);
    
    return true;
  }
}    //   namespace nlo


#endif
