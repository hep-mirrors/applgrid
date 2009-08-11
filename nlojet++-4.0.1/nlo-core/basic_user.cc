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

//   nlojet++ includes
#include <bits/nlo-basic_user.h>


namespace nlo {

	
  namespace __helper_basic_user {
       
    //   step function
    double _S_step_dbl_fn(double x, double y) {
      return x > y ? 1.0 : 0.0;
    }
  
    std::pointer_to_binary_function<double, double, double> _G_step_dbl(_S_step_dbl_fn);
    

    //   dirac delta function using with the histpoint1d
    double _S_dirac_fn(double x, const histpoint1d& y) {
      return y.xmin <= x && x < y.xmax ? 1.0/(y.xmax-y.xmin) : 0.0;
    }
    
    //   step function using with the histpoint1d
    double _S_step_fn(double x, const histpoint1d& y) {
      return x > y.xmin ? 1.0 : 0.0;
    }

    std::pointer_to_binary_function<double, const histpoint1d&, double> _G_dirac(_S_dirac_fn);
    std::pointer_to_binary_function<double, const histpoint1d&, double> _G_step(_S_step_fn);


    //  dirac and step distributions using with the histpoint2d
    double _S_dist2d_dd_fn(__helper_point2d p, const histpoint2d& y) 
    {
      if(!(y.xmin <= p.x && p.x < y.xmax)) return 0.0;
      if(!(y.ymin <= p.y && p.y < y.ymax)) return 0.0;
      return 1.0/(y.xmax-y.xmin)/(y.ymax-y.ymin);
    }

    double _S_dist2d_hd_fn(__helper_point2d p, const histpoint2d& y) 
    {
      if(p.x <= y.xmin) return 0.0;
      if(!(y.ymin <= p.y && p.y < y.ymax)) return 0.0;
      return 1.0/(y.ymax-y.ymin);
    }
    
    double _S_dist2d_dh_fn(__helper_point2d p, const histpoint2d& y) 
    {
      if(!(y.xmin <= p.x && p.x < y.xmax)) return 0.0;
      if(p.y <= y.ymin) return 0.0;
      return 1.0/(y.xmax-y.xmin);
    }

     double _S_dist2d_hh_fn(__helper_point2d p, const histpoint2d& y) 
    {
      if(p.x <= y.xmin) return 0.0;
      if(p.y <= y.ymin) return 0.0;
      return 1.0;
    }

    std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_dd(_S_dist2d_dd_fn);
    std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_dh(_S_dist2d_dh_fn);
    std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_hd(_S_dist2d_hd_fn);
    std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_hh(_S_dist2d_hh_fn);
  }  //  namespace __helper_basic_user
}
