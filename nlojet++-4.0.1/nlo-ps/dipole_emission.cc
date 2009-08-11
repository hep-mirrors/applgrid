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

// Standard C inludes
#include <cmath>

// nlo++ includes
#include "bits/psg-dipole.h"
#include "bits/hep-exception.h"


#define cutoff      1.0e-9
#define __NLO_PI    3.14159265358979323846
#define __NLO_2PI   6.28318530717958647692
#define __NLO_8PI2  78.95683520871486895054


namespace nlo {



  double dipole_emission::
  _S_Gfunc(double beta, double eps, double x) 
  {
    if(eps == 1.0) return std::pow(x, 1.0/beta)/(beta*x);
    else {
      double N = beta*std::pow(eps, 1.0/beta) - std::log(eps); 
      if(x <= eps) return std::pow(x, 1.0/beta)/(N*x);
      else return 1.0/(N*x);
    }
  }
  
  double dipole_emission::
  _S_gfunc(double beta, double eps, double r1, double r2)
  {
    if(eps == 1.0) return std::pow(r1, beta);
    else {
      double leps = - std::log(eps);
      double w = 1.0/(1.0 + leps/(beta*std::pow(eps, 1.0/beta))); 
      
      if(r2 <= w) return eps*std::pow(r1, beta);
      else return eps*std::exp(r1*leps);
    }
  }

  double dipole_emission::
  _S_Ffunc(double beta, double p, double x, double eta) 
  {
    double e = x*eta, ox = 1.0-x;
    double oe = 1.0-e, op = 1.0-p;
    
    if(p == 1.0) return x*std::pow(ox/oe, 1.0/beta)/(beta*ox);
    else {
      double N;
      if(p < oe) N = std::log(op*oe/(e*p)) + beta*std::pow(p, 1.0/beta);
      else N = beta*std::pow(oe, 1.0/beta);
      
      if(x >= 1.0-p) return x*std::pow(ox, 1.0/beta)/(N*ox);
      else return 1.0/(N*ox);
    }
  }
  
  double dipole_emission::
  _S_ffunc(double beta, double p, double r1, double r2, double e) 
  {
    double oe = 1.0-e, op = 1.0-p;
    if(p >= oe) return 1.0-oe*std::pow(r1, beta);
    else {
      double leps = std::log(op*oe/(e*p));
      double w = 1.0/(1.0 + leps/(beta*std::pow(p, 1.0/beta)));
      
      if(r2 <= w) return 1.0 - p*std::pow(r1, beta);
      else {
	double z = e*std::exp(r1*leps)/oe;
	return z/(1.0 + z);
      }
    }
  }


  void dipole_emission::
  gendip_fff(const random_generator& rng, double beta, double eps, reference pi, reference pj, reference pk) 
  {
    double y = _S_gfunc(beta, eps, rng(), rng());
    double z = 1.0 - _S_gfunc(beta, eps, rng(), rng());
    
    if(z < cutoff || 1.0 - z < cutoff) throw numeric_error(); 
    if(y < cutoff || 1.0 - y < cutoff) throw numeric_error(); 
       
    lorentzvector<double> qik   = pi + pk;
    lorentzvector<double> qk    = pk;
    threevector<double>   boost = qik.boostVector();
    
    double E  = std::sqrt(0.5*(pi*pk));
    double Ei = E*((1.0 - z)*y + z);
    double pz = E*((1.0 - z)*y - z);
    double pt = E*std::sqrt(4.0*z*(1.0 - z)*y);
    double ph = __NLO_2PI*rng();
    lorentzvector<double> qi(pt*std::cos(ph), pt*std::sin(ph), pz, Ei);
    
    qk.boost(-boost);
    qi.rotateY(qk.theta());
    qi.rotateZ(qk.phi());
    qi.boost(boost);
    
    pk *= 1.0 - y;
    pj  = qik - qi - pk;
    pi  = qi;
  }

  void dipole_emission::
  gendip_ffi(const random_generator& rng, double beta, double eps, double xmin, reference pi, reference pj, reference pa)
  {
    double x = _S_ffunc(beta, eps, rng(), rng(), xmin);
    double z = 1.0 - _S_gfunc(beta, eps, rng(), rng());
    
    if(z < cutoff || 1.0 - z < cutoff) throw numeric_error();
    if(x < cutoff || 1.0 - x < cutoff) throw numeric_error();
    
    lorentzvector<double> q = (pa /= x);
    threevector<double>   boost = (q + pi).boostVector();
    
    double E  = std::sqrt(0.5*(pi*q));
    double Ei = E*((1.0 - z)*(1.0 - x) + z);
    double pz = E*((1.0 - z)*(1.0 - x) - z);
    double pt = 2.0*E*std::sqrt(z*(1.0 - z)*(1.0 - x));
    double ph = __NLO_2PI*rng();
    lorentzvector<double> qi(pt*std::cos(ph), pt*std::sin(ph), pz, Ei); 
    
    q.boost(-boost);  
    qi.rotateY(q.theta());
    qi.rotateZ(q.phi());
    qi.boost(boost);
    
    pj = pa*(1.0-x) - qi + pi;
    pi = qi;
  }

  void dipole_emission::
  gendip_ifi(const random_generator& rng, double beta, double eps, double xmin, reference pa, 
	     reference pi, reference pb, boost_ifi& bst)
  {
    double u = _S_gfunc(beta, eps, rng(), rng());
    double x = _S_ffunc(beta, eps, rng(), rng(), xmin);
    
    if(u < cutoff || 1.0 - u < cutoff) throw numeric_error();
    if(x < cutoff || 1.0 - x < cutoff) throw numeric_error();
    
    lorentzvector<double> q  = (pa /= x);
    threevector<double>   boost((q + pb).boostVector());
    
    double E   = sqrt(0.5*(pa*pb));
    double Ei  = E*(1.0-x);
    double pz  = Ei*(1.0-2.0*u);
    double pt  = Ei*sqrt(4.0*u*(1.0-u));
    double phi = __NLO_2PI*rng();
    lorentzvector<double> qi(pt*cos(phi), pt*sin(phi), pz, Ei); 
    
    q.boost(-boost);  
    qi.rotateY(q.theta());
    qi.rotateZ(q.phi());
    qi.boost(boost);
    pi = qi;
    
    //  set up  the boost
    bst.set(pa+pb-pi, x*pa+pb);
  }

  double dipole_emission::
  jacobi_fff(double beta, double eps, const_reference pi, const_reference pj, const_reference pk)
  {
    double sij = pi*pj, sik = pi*pk, sjk = pj*pk;
    double y = sij/(sik + sjk + sij), z = sik/(sik + sjk);
    
    if(1.0 - z < cutoff || y < cutoff) throw numeric_error();
    return __NLO_8PI2*_S_Gfunc(beta, eps, 1.0-z)*_S_Gfunc(beta, eps, y)/(sik+sjk);
  }
  
  double dipole_emission::
  jacobi_ffi(double beta, double eps, double eta, const_reference pi, const_reference pj, const_reference pa)
  {
    double sai = pi*pa, saj = pj*pa;
    double x = 1.0 - (pi*pj)/(sai + saj), z = sai/(sai + saj);
    
    if(1.0 - z < cutoff || 1.0 - x < cutoff) throw numeric_error(); 
    return __NLO_8PI2*_S_Gfunc(beta, eps, 1.0-z)*_S_Ffunc(beta, eps, x, eta)/(sai+saj);
  }
  
  double dipole_emission::
  jacobi_ifi(double beta, double eps, double eta, const_reference pa, const_reference pi, const_reference pb)
  {
    double sab = pa*pb, sai = pi*pa, sbi = pi*pb;
    double x = 1.0 - (sai + sbi)/sab, u = sai/(sai + sbi);
    
    if(u < cutoff || 1.0 - x < cutoff) throw numeric_error(); 
    return __NLO_8PI2*_S_Gfunc(beta, eps, u)*_S_Ffunc(beta, eps, x, eta)/(sai+sbi);
  }
}
