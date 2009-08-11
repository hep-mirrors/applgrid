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
#ifndef __NLO_NLO_SPLIT_H__
#define __NLO_NLO_SPLIT_H__ 1

//   Standard includes
#include <complex>
#include <utility>

//   nlo includes
#include <lorentzvector.h>


namespace nlo {


  //   Splitting functions of the final state singularities.
  template<class _Lv>
  struct split_fin 
  {
    //  types
    typedef _Lv lorentzvector_type;
    typedef typename _Lv::value_type value_type;
    typedef std::pair<value_type, std::complex<value_type> > pair_type;
    
    //  destructor
    virtual ~split_fin() {}
    
    //  set up the parameters    
    virtual void set(const _Lv&, const _Lv&, const _Lv&) = 0;
    
    //  splitting functions
    virtual pair_type Vqg() const = 0;
    virtual pair_type Vqa() const = 0;
    virtual pair_type Vgg() const = 0;
  };

  //   Splitting functions of the initial state singularities.
  template<class _Lv>
  struct split_ini 
  {
    //  types
    typedef _Lv lorentzvector_type;
    typedef typename _Lv::value_type value_type;
    typedef std::pair<value_type, std::complex<value_type> > pair_type;
    
    //  destructor
    virtual ~split_ini() {}
    
    //  set up the parameters    
    virtual void set(const _Lv&, const _Lv&, const _Lv&) = 0;

    //  splitting functions
    virtual pair_type Vqg() const = 0;
    virtual pair_type Vga() const = 0;
    virtual pair_type Vqq() const = 0;
    virtual pair_type Vgg() const = 0;
  };

  
  
  class splitfff 
    : public split_fin<lorentzvector<double> > 
  {
    //    private types
    typedef lorentzvector<double> _Lv;
    
  public:
    //  destructor
    ~splitfff() {}
    
    //  set the data members
    void set(const _Lv&, const _Lv&, const _Lv&);

    //  calculate the splitting functions
    pair_type Vqg() const;
    pair_type Vqa() const;
    pair_type Vgg() const;

  private:
    double yijk, zi, sij;
    std::complex<double> hij;
  };
  

  class splitffi 
    : public split_fin<lorentzvector<double> > 
  {
    //    private types
    typedef lorentzvector<double> _Lv;
    
  public:
    //  destructor
    ~splitffi() {}

    //  set the data members
    void set(const _Lv&, const _Lv&, const _Lv&);

    //  calculate the splitting functions
    pair_type Vqg() const;
    pair_type Vqa() const;
    pair_type Vgg() const;

  private:
    double xija, zi, sij;
    std::complex<double> hij;
  };

  class splitiff 
    : public split_ini<lorentzvector<double> > 
  {
    //    private types
    typedef lorentzvector<double> _Lv;
    
  public:
    //  destructor
    ~splitiff() {}

    //  set the data members
    void set(const _Lv&, const _Lv&, const _Lv&);

    //  calculate the splitting functions
    pair_type Vqg() const;
    pair_type Vga() const;
    pair_type Vqq() const;
    pair_type Vgg() const;

  private:
    double xika, ui, sai;
    std::complex<double> hai;
  };

  class splitifi 
    : public split_ini<lorentzvector<double> >
  {
    //    private types
    typedef lorentzvector<double> _Lv;
    
  public:
    //  destructor
    ~splitifi() {}
    
    //  set the data members
    void set(const _Lv&, const _Lv&, const _Lv&);
    
    //  calculate the splitting functions
    pair_type Vqg() const;
    pair_type Vga() const;
    pair_type Vqq() const;
    pair_type Vgg() const;
    
  private:
    double xiab, sai;
    std::complex<double> hai;
  };
}  //  namespace nlo

#endif
