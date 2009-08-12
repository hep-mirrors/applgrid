//  Copyright (C) 2004 Zoltan Nagy
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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

//  This header file was generated automatically by mk-basic_user_result-header.
#ifndef __NLO_NLO_BASIC_USER_RESULT_H__
#define __NLO_NLO_BASIC_USER_RESULT_H__ 1

//   stdc++ includes
#include <string>
#include <iostream>


namespace nlo {

  class _NullDistBook{};

  template<class _DistBook0=_NullDistBook, class _DistBook1=_NullDistBook, 
           class _DistBook2=_NullDistBook, class _DistBook3=_NullDistBook, 
           class _DistBook4=_NullDistBook, class _DistBook5=_NullDistBook, 
           class _DistBook6=_NullDistBook, class _DistBook7=_NullDistBook, 
           class _DistBook8=_NullDistBook, class _DistBook9=_NullDistBook>
  class basic_user_result;


  template<class _DistBook0>
  class basic_user_result<_DistBook0, _NullDistBook, _NullDistBook, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
  };


  template<class _DistBook0>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _NullDistBook, _NullDistBook, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0;
    } else {
       read(is, x.book0);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1>
  class basic_user_result<_DistBook0, _DistBook1, _NullDistBook, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
  };


  template<class _DistBook0, class _DistBook1>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _NullDistBook, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1;
    } else {
       read(is, x.book0); read(is, x.book1);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2>
  class basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3>
  class basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                          _DistBook3, _NullDistBook, _NullDistBook, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;
    typedef _DistBook3 distbook3_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2; 
       book3 += x.book3;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       book3.assadd_weighted(x.book3);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
    distbook3_type book3;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _DistBook3, _NullDistBook, _NullDistBook, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook3).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2; is>>x.book3;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2); 
       read(is, x.book3);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4>
  class basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                          _DistBook3, _DistBook4, _NullDistBook, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;
    typedef _DistBook3 distbook3_type;
    typedef _DistBook4 distbook4_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2; 
       book3 += x.book3; book4 += x.book4;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       book3.assadd_weighted(x.book3);
       book4.assadd_weighted(x.book4);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
    distbook3_type book3;
    distbook4_type book4;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _DistBook3, _DistBook4, _NullDistBook, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook3).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook4).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2; is>>x.book3; is>>x.book4;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2); 
       read(is, x.book3); read(is, x.book4);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5>
  class basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                          _DistBook3, _DistBook4, _DistBook5, 
                          _NullDistBook, _NullDistBook, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;
    typedef _DistBook3 distbook3_type;
    typedef _DistBook4 distbook4_type;
    typedef _DistBook5 distbook5_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2; 
       book3 += x.book3; book4 += x.book4; book5 += x.book5;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       book3.assadd_weighted(x.book3);
       book4.assadd_weighted(x.book4);
       book5.assadd_weighted(x.book5);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
    distbook3_type book3;
    distbook4_type book4;
    distbook5_type book5;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _DistBook3, _DistBook4, _DistBook5, 
                   _NullDistBook, _NullDistBook, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook3).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook4).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook5).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2; is>>x.book3; is>>x.book4; 
       is>>x.book5;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2); 
       read(is, x.book3); read(is, x.book4); read(is, x.book5);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6>
  class basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                          _DistBook3, _DistBook4, _DistBook5, 
                          _DistBook6, _NullDistBook, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;
    typedef _DistBook3 distbook3_type;
    typedef _DistBook4 distbook4_type;
    typedef _DistBook5 distbook5_type;
    typedef _DistBook6 distbook6_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2; 
       book3 += x.book3; book4 += x.book4; book5 += x.book5; 
       book6 += x.book6;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       book3.assadd_weighted(x.book3);
       book4.assadd_weighted(x.book4);
       book5.assadd_weighted(x.book5);
       book6.assadd_weighted(x.book6);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
    distbook3_type book3;
    distbook4_type book4;
    distbook5_type book5;
    distbook6_type book6;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _DistBook3, _DistBook4, _DistBook5, 
                   _DistBook6, _NullDistBook, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook3).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook4).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook5).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook6).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2; is>>x.book3; is>>x.book4; 
       is>>x.book5; is>>x.book6;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2); 
       read(is, x.book3); read(is, x.book4); read(is, x.book5); 
       read(is, x.book6);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6, class _DistBook7>
  class basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                          _DistBook3, _DistBook4, _DistBook5, 
                          _DistBook6, _DistBook7, _NullDistBook, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;
    typedef _DistBook3 distbook3_type;
    typedef _DistBook4 distbook4_type;
    typedef _DistBook5 distbook5_type;
    typedef _DistBook6 distbook6_type;
    typedef _DistBook7 distbook7_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2; 
       book3 += x.book3; book4 += x.book4; book5 += x.book5; 
       book6 += x.book6; book7 += x.book7;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       book3.assadd_weighted(x.book3);
       book4.assadd_weighted(x.book4);
       book5.assadd_weighted(x.book5);
       book6.assadd_weighted(x.book6);
       book7.assadd_weighted(x.book7);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
    distbook3_type book3;
    distbook4_type book4;
    distbook5_type book5;
    distbook6_type book6;
    distbook7_type book7;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6, class _DistBook7>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _DistBook3, _DistBook4, _DistBook5, 
                   _DistBook6, _DistBook7, _NullDistBook, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook3).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook4).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook5).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook6).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook7).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2; is>>x.book3; is>>x.book4; 
       is>>x.book5; is>>x.book6; is>>x.book7;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2); 
       read(is, x.book3); read(is, x.book4); read(is, x.book5); 
       read(is, x.book6); read(is, x.book7);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6, class _DistBook7, class _DistBook8>
  class basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                          _DistBook3, _DistBook4, _DistBook5, 
                          _DistBook6, _DistBook7, _DistBook8, 
                          _NullDistBook>
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;
    typedef _DistBook3 distbook3_type;
    typedef _DistBook4 distbook4_type;
    typedef _DistBook5 distbook5_type;
    typedef _DistBook6 distbook6_type;
    typedef _DistBook7 distbook7_type;
    typedef _DistBook8 distbook8_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2; 
       book3 += x.book3; book4 += x.book4; book5 += x.book5; 
       book6 += x.book6; book7 += x.book7; book8 += x.book8;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       book3.assadd_weighted(x.book3);
       book4.assadd_weighted(x.book4);
       book5.assadd_weighted(x.book5);
       book6.assadd_weighted(x.book6);
       book7.assadd_weighted(x.book7);
       book8.assadd_weighted(x.book8);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
    distbook3_type book3;
    distbook4_type book4;
    distbook5_type book5;
    distbook6_type book6;
    distbook7_type book7;
    distbook8_type book8;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6, class _DistBook7, class _DistBook8>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _DistBook3, _DistBook4, _DistBook5, 
                   _DistBook6, _DistBook7, _DistBook8, 
                   _NullDistBook>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook3).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook4).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook5).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook6).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook7).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook8).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2; is>>x.book3; is>>x.book4; 
       is>>x.book5; is>>x.book6; is>>x.book7; is>>x.book8;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2); 
       read(is, x.book3); read(is, x.book4); read(is, x.book5); 
       read(is, x.book6); read(is, x.book7); read(is, x.book8);
    }

    return nsaved;
  }



  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6, class _DistBook7, class _DistBook8, 
           class _DistBook9>
  class basic_user_result
  {
  public:
    //   public types
    typedef _DistBook0 distbook0_type;
    typedef _DistBook1 distbook1_type;
    typedef _DistBook2 distbook2_type;
    typedef _DistBook3 distbook3_type;
    typedef _DistBook4 distbook4_type;
    typedef _DistBook5 distbook5_type;
    typedef _DistBook6 distbook6_type;
    typedef _DistBook7 distbook7_type;
    typedef _DistBook8 distbook8_type;
    typedef _DistBook9 distbook9_type;

    //   computed assignments
    basic_user_result& operator+=(const basic_user_result& x)
    {
       book0 += x.book0; book1 += x.book1; book2 += x.book2; 
       book3 += x.book3; book4 += x.book4; book5 += x.book5; 
       book6 += x.book6; book7 += x.book7; book8 += x.book8; 
       book9 += x.book9;
       return *this;
    }

    basic_user_result& assadd_weighted(const basic_user_result& x)
    {
       book0.assadd_weighted(x.book0);
       book1.assadd_weighted(x.book1);
       book2.assadd_weighted(x.book2);
       book3.assadd_weighted(x.book3);
       book4.assadd_weighted(x.book4);
       book5.assadd_weighted(x.book5);
       book6.assadd_weighted(x.book6);
       book7.assadd_weighted(x.book7);
       book8.assadd_weighted(x.book8);
       book9.assadd_weighted(x.book9);
       return *this;
    }

    //   public acces to the distbook objects
    distbook0_type book0;
    distbook1_type book1;
    distbook2_type book2;
    distbook3_type book3;
    distbook4_type book4;
    distbook5_type book5;
    distbook6_type book6;
    distbook7_type book7;
    distbook8_type book8;
    distbook9_type book9;
  };


  template<class _DistBook0, class _DistBook1, class _DistBook2, 
           class _DistBook3, class _DistBook4, class _DistBook5, 
           class _DistBook6, class _DistBook7, class _DistBook8, 
           class _DistBook9>
  unsigned long read_basic_user_output(std::istream& is,
    basic_user_result<_DistBook0, _DistBook1, _DistBook2, 
                   _DistBook3, _DistBook4, _DistBook5, 
                   _DistBook6, _DistBook7, _DistBook8, 
                   _DistBook9>& x)
  {
    bool txt, flag = false;
    unsigned long int nsaved;
    std::basic_string<char> distbook;
    extern void __helper_basic_user_set_warning();

    is>>txt;
    is>>nsaved;

    is>>distbook; if(typeid(_DistBook0).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook1).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook2).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook3).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook4).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook5).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook6).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook7).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook8).name() != distbook) flag = true;
    is>>distbook; if(typeid(_DistBook9).name() != distbook) flag = true;

    if(flag) __helper_basic_user_set_warning();

    if(txt) {
       is>>x.book0; is>>x.book1; is>>x.book2; is>>x.book3; is>>x.book4; 
       is>>x.book5; is>>x.book6; is>>x.book7; is>>x.book8; is>>x.book9;
    } else {
       read(is, x.book0); read(is, x.book1); read(is, x.book2); 
       read(is, x.book3); read(is, x.book4); read(is, x.book5); 
       read(is, x.book6); read(is, x.book7); read(is, x.book8); 
       read(is, x.book9);
    }

    return nsaved;
  }

}  // namespace nlo

#endif
