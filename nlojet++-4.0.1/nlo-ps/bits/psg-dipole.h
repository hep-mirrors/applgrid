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
#ifndef __NLO_PSG_DIPOLE_H__
#define __NLO_PSG_DIPOLE_H__ 1

// nlo++ includes
#include <bits/hep-lorentzvector.h>
#include <bits/hep-exception.h>
#include <bits/hep-rng.h>



namespace nlo {
  
  /// \brief Class for generating dipole emmissons.
  /// 
  /// This class is a collection of some useful function to generate
  /// dipole emission. The algorithms are based on the dipole phase 
  /// space factorization formulae. Using these functions one can generate 
  /// extra particle in the final state emitting by an initial or final 
  /// particle. 
  class dipole_emission 
  {
  public:   
    /// Alias to the lorentz vector type.
    typedef lorentzvector<double> lorentzvector_type;
    
    /// Reference type of the lorentz vetor type.
    typedef lorentzvector_type& reference;
    
    /// Constant reference of the lorentz vector type.
    typedef const lorentzvector<double>& const_reference;
    
    /// \brief This class defines the lorentz boost in the case ....
    class boost_ifi 
    {
    public:
      boost_ifi() {}
      boost_ifi(const_reference k, const_reference kb) {
		this -> set(k, kb);
      }
      
      void set(const_reference k, const_reference kb) {
		K = k; Kb = kb; Q = k + kb; K2 = K.mag2(); Q2 = Q.mag2();
      }
      
      lorentzvector<double> operator()(const_reference pj) const {
		return pj - Q*(2.0*(Q*pj)/Q2) + K*(2.0*(Kb*pj)/K2);
      }
      
    private:
      lorentzvector<double> K, Kb, Q;
      double K2, Q2;
    };
    
    /// \fn static double gendip_fff(const random_generator&, double eps, reference pi, reference pj, reference pk);
    /// \brief This function generate dipole emission in the final state 
    /// when the recoiled parton (spectator) is also in the final state.
    static void gendip_fff(const random_generator&, double, double, reference, reference, reference);
    static void gendip_ffi(const random_generator&, double, double, double, reference, reference, reference); 
    static void gendip_ifi(const random_generator&, double, double, double, reference, reference, reference, boost_ifi&);
    
    static double jacobi_fff(double, double, const_reference, const_reference, const_reference);
    static double jacobi_ffi(double, double, double, const_reference, const_reference, const_reference);
    static double jacobi_ifi(double, double, double, const_reference, const_reference, const_reference);

    template<typename _InputIterator>
    static double jacobi_fff(double, double, _InputIterator, _InputIterator);

    template<typename _InputIterator>
    static double jacobi_ffi(double, double, double, const_reference, _InputIterator, _InputIterator);
    
    template<typename _InputIterator>
    static double jacobi_ifi(double, double, double, const_reference, const_reference, _InputIterator, _InputIterator);
    
  private:
    ///  private static members
    static double _S_Gfunc(double, double, double);
    static double _S_Ffunc(double, double, double, double);
    static double _S_gfunc(double, double, double, double);
    static double _S_ffunc(double, double, double, double, double);
  };

 
  template<typename _InputIterator> double dipole_emission::
  jacobi_fff(double beta, double eps, _InputIterator first, _InputIterator last)
  {
    static const double cutoff = 1e-9;
    double sij, sik, sjk, y, z, tmp, jac = 0.0;
    
    for(_InputIterator pi = first; pi < last; pi++)
      for(_InputIterator pj = pi+1; pj < last; pj++) {
	sij = (*pi)*(*pj);
	for(_InputIterator pk = first; pk < last; pk++)
	  if(pk != pj && pk != pi) {
	    sik = (*pi)*(*pk);
	    sjk = (*pj)*(*pk);
	    
	    y = sij/(sik + sjk + sij);
	    z = sik/(sik + sjk);
	    
	    if(z < cutoff || 1.0 - z < cutoff || y < cutoff) 
	      throw numeric_error();
	    
	    tmp = _S_Gfunc(beta, eps, 1.0-z) + _S_Gfunc(beta, eps, z);
	    jac += tmp*_S_Gfunc(beta, eps, y)/(sik+sjk);
	  }
      }
    
    return 78.95683520871486895054*jac;
  }
  
  template<typename _InputIterator> double dipole_emission::
  jacobi_ffi(double beta, double eps, double eta, const_reference pa,
	     _InputIterator first, _InputIterator last)
  {
    static const double cutoff = 1e-9;
    double sij, sai, saj, x, z, tmp, jac = 0.0;
    
    for(_InputIterator pi = first; pi < last; pi++) {
      sai = (*pi)*pa; 
      for(_InputIterator pj = pi + 1; pj < last; pj++) {
	saj = (*pj)*pa;
	sij = (*pi)*(*pj);
	
	x = 1.0 - ((*pi)*(*pj))/(sai + saj);
        z = sai/(sai + saj);
	
	if(z < cutoff || 1.0 - z < cutoff || 1.0 - x < cutoff) 
	  throw numeric_error(); 
	
	tmp = _S_Gfunc(beta, eps, 1.0-z) + _S_Gfunc(beta, eps, z);
	jac += tmp*_S_Ffunc(beta, eps, x, eta)/(sai+saj);
      }
    }  
    
    return 78.95683520871486895054*jac;
  }
  
  template<typename _InputIterator> double dipole_emission::
  jacobi_ifi(double beta, double eps, double eta, const_reference pa, const_reference pb,
	     _InputIterator first, _InputIterator last)
  {
    static const double cutoff = 1e-9;
    double sab = pa*pb, sai, sbi, x, u, jac = 0.0;
    
    for(_InputIterator pi = first; pi < last; pi++) {
      sai = (*pi)*pa; sbi = (*pi)*pb;
      x = 1.0 - (sai + sbi)/sab;
      u = sai/(sai + sbi);
      
      if(u < cutoff || 1.0 - u < cutoff || 1.0 - x < cutoff) 
	throw numeric_error(); 
      
      jac += _S_Gfunc(beta, eps, u)*_S_Ffunc(beta, eps, x, eta)/(sai+sbi);
    }
    
    return 78.95683520871486895054*jac;
  }
}

#endif
