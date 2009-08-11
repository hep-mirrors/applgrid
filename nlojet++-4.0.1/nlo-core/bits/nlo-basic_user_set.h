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

//  This header file was generated automatically by mk-basic_user_set-header.
#ifndef __NLO_NLO_BASIC_USER_SET_H__
#define __NLO_NLO_BASIC_USER_SET_H__ 1

//   nlojet++ includes
#include <bits/nlo-basic_user.h>


namespace nlo {

  class _NullUser{};

  template<class _User0=_NullUser, class _User1=_NullUser, 
           class _User2=_NullUser, class _User3=_NullUser, 
           class _User4=_NullUser, class _User5=_NullUser, 
           class _User6=_NullUser, class _User7=_NullUser, 
           class _User8=_NullUser, class _User9=_NullUser>
  class basic_user_set;


  template<class _User0>
  class basic_user_set<_User0, _NullUser, _NullUser, _NullUser, _NullUser, 
                       _NullUser, _NullUser, _NullUser, _NullUser, _NullUser>
    : public _User0
  {
  public:
    //   public types
    typedef _User0 user0_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
    }
  };



  template<class _User0, class _User1>
  class basic_user_set<_User0, _User1, _NullUser, _NullUser, _NullUser, 
                       _NullUser, _NullUser, _NullUser, _NullUser, _NullUser>
    : public _User0, public _User1
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2>
  class basic_user_set<_User0, _User1, _User2, _NullUser, _NullUser, 
                       _NullUser, _NullUser, _NullUser, _NullUser, _NullUser>
    : public _User0, public _User1, public _User2
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2, class _User3>
  class basic_user_set<_User0, _User1, _User2, _User3, _NullUser, 
                       _NullUser, _NullUser, _NullUser, _NullUser, _NullUser>
    : public _User0, public _User1, public _User2, public _User3
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;
    typedef _User3 user3_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
      this ->  _User3::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
      this -> _User3::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
      this -> _User3::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
      this -> _User3::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2, class _User3, 
           class _User4>
  class basic_user_set<_User0, _User1, _User2, _User3, _User4, 
                       _NullUser, _NullUser, _NullUser, _NullUser, _NullUser>
    : public _User0, public _User1, public _User2, public _User3, 
      public _User4
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;
    typedef _User3 user3_type;
    typedef _User4 user4_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
      this ->  _User3::operations_at_the_end_of_event();
      this ->  _User4::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
      this -> _User3::write_typeinfo(os);
      this -> _User4::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
      this -> _User3::write_result_txt(os);
      this -> _User4::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
      this -> _User3::write_result_bin(os);
      this -> _User4::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2, class _User3, 
           class _User4, class _User5>
  class basic_user_set<_User0, _User1, _User2, _User3, _User4, 
                       _User5, _NullUser, _NullUser, _NullUser, _NullUser>
    : public _User0, public _User1, public _User2, public _User3, 
      public _User4, public _User5
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;
    typedef _User3 user3_type;
    typedef _User4 user4_type;
    typedef _User5 user5_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
      this ->  _User3::operations_at_the_end_of_event();
      this ->  _User4::operations_at_the_end_of_event();
      this ->  _User5::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
      this -> _User3::write_typeinfo(os);
      this -> _User4::write_typeinfo(os);
      this -> _User5::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
      this -> _User3::write_result_txt(os);
      this -> _User4::write_result_txt(os);
      this -> _User5::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
      this -> _User3::write_result_bin(os);
      this -> _User4::write_result_bin(os);
      this -> _User5::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2, class _User3, 
           class _User4, class _User5, class _User6>
  class basic_user_set<_User0, _User1, _User2, _User3, _User4, 
                       _User5, _User6, _NullUser, _NullUser, _NullUser>
    : public _User0, public _User1, public _User2, public _User3, 
      public _User4, public _User5, public _User6
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;
    typedef _User3 user3_type;
    typedef _User4 user4_type;
    typedef _User5 user5_type;
    typedef _User6 user6_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
      this ->  _User3::operations_at_the_end_of_event();
      this ->  _User4::operations_at_the_end_of_event();
      this ->  _User5::operations_at_the_end_of_event();
      this ->  _User6::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
      this -> _User3::write_typeinfo(os);
      this -> _User4::write_typeinfo(os);
      this -> _User5::write_typeinfo(os);
      this -> _User6::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
      this -> _User3::write_result_txt(os);
      this -> _User4::write_result_txt(os);
      this -> _User5::write_result_txt(os);
      this -> _User6::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
      this -> _User3::write_result_bin(os);
      this -> _User4::write_result_bin(os);
      this -> _User5::write_result_bin(os);
      this -> _User6::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2, class _User3, 
           class _User4, class _User5, class _User6, class _User7>
  class basic_user_set<_User0, _User1, _User2, _User3, _User4, 
                       _User5, _User6, _User7, _NullUser, _NullUser>
    : public _User0, public _User1, public _User2, public _User3, 
      public _User4, public _User5, public _User6, public _User7
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;
    typedef _User3 user3_type;
    typedef _User4 user4_type;
    typedef _User5 user5_type;
    typedef _User6 user6_type;
    typedef _User7 user7_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
      this ->  _User3::operations_at_the_end_of_event();
      this ->  _User4::operations_at_the_end_of_event();
      this ->  _User5::operations_at_the_end_of_event();
      this ->  _User6::operations_at_the_end_of_event();
      this ->  _User7::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
      this -> _User3::write_typeinfo(os);
      this -> _User4::write_typeinfo(os);
      this -> _User5::write_typeinfo(os);
      this -> _User6::write_typeinfo(os);
      this -> _User7::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
      this -> _User3::write_result_txt(os);
      this -> _User4::write_result_txt(os);
      this -> _User5::write_result_txt(os);
      this -> _User6::write_result_txt(os);
      this -> _User7::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
      this -> _User3::write_result_bin(os);
      this -> _User4::write_result_bin(os);
      this -> _User5::write_result_bin(os);
      this -> _User6::write_result_bin(os);
      this -> _User7::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2, class _User3, 
           class _User4, class _User5, class _User6, class _User7, 
           class _User8>
  class basic_user_set<_User0, _User1, _User2, _User3, _User4, 
                       _User5, _User6, _User7, _User8, _NullUser>
    : public _User0, public _User1, public _User2, public _User3, 
      public _User4, public _User5, public _User6, public _User7, 
      public _User8
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;
    typedef _User3 user3_type;
    typedef _User4 user4_type;
    typedef _User5 user5_type;
    typedef _User6 user6_type;
    typedef _User7 user7_type;
    typedef _User8 user8_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
      this ->  _User3::operations_at_the_end_of_event();
      this ->  _User4::operations_at_the_end_of_event();
      this ->  _User5::operations_at_the_end_of_event();
      this ->  _User6::operations_at_the_end_of_event();
      this ->  _User7::operations_at_the_end_of_event();
      this ->  _User8::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
      this -> _User3::write_typeinfo(os);
      this -> _User4::write_typeinfo(os);
      this -> _User5::write_typeinfo(os);
      this -> _User6::write_typeinfo(os);
      this -> _User7::write_typeinfo(os);
      this -> _User8::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
      this -> _User3::write_result_txt(os);
      this -> _User4::write_result_txt(os);
      this -> _User5::write_result_txt(os);
      this -> _User6::write_result_txt(os);
      this -> _User7::write_result_txt(os);
      this -> _User8::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
      this -> _User3::write_result_bin(os);
      this -> _User4::write_result_bin(os);
      this -> _User5::write_result_bin(os);
      this -> _User6::write_result_bin(os);
      this -> _User7::write_result_bin(os);
      this -> _User8::write_result_bin(os);
    }
  };



  template<class _User0, class _User1, class _User2, class _User3, 
           class _User4, class _User5, class _User6, class _User7, 
           class _User8, class _User9>
  class basic_user_set
    : public _User0, public _User1, public _User2, public _User3, 
      public _User4, public _User5, public _User6, public _User7, 
      public _User8, public _User9
  {
  public:
    //   public types
    typedef _User0 user0_type;
    typedef _User1 user1_type;
    typedef _User2 user2_type;
    typedef _User3 user3_type;
    typedef _User4 user4_type;
    typedef _User5 user5_type;
    typedef _User6 user6_type;
    typedef _User7 user7_type;
    typedef _User8 user8_type;
    typedef _User9 user9_type;

  private:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      this ->  _User0::operations_at_the_end_of_event();
      this ->  _User1::operations_at_the_end_of_event();
      this ->  _User2::operations_at_the_end_of_event();
      this ->  _User3::operations_at_the_end_of_event();
      this ->  _User4::operations_at_the_end_of_event();
      this ->  _User5::operations_at_the_end_of_event();
      this ->  _User6::operations_at_the_end_of_event();
      this ->  _User7::operations_at_the_end_of_event();
      this ->  _User8::operations_at_the_end_of_event();
      this ->  _User9::operations_at_the_end_of_event();
    }

    //   write some information about the type of the sample
    void write_typeinfo(std::ostream& os) {
      this -> _User0::write_typeinfo(os);
      this -> _User1::write_typeinfo(os);
      this -> _User2::write_typeinfo(os);
      this -> _User3::write_typeinfo(os);
      this -> _User4::write_typeinfo(os);
      this -> _User5::write_typeinfo(os);
      this -> _User6::write_typeinfo(os);
      this -> _User7::write_typeinfo(os);
      this -> _User8::write_typeinfo(os);
      this -> _User9::write_typeinfo(os);
    }

    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      this -> _User0::write_result_txt(os);
      this -> _User1::write_result_txt(os);
      this -> _User2::write_result_txt(os);
      this -> _User3::write_result_txt(os);
      this -> _User4::write_result_txt(os);
      this -> _User5::write_result_txt(os);
      this -> _User6::write_result_txt(os);
      this -> _User7::write_result_txt(os);
      this -> _User8::write_result_txt(os);
      this -> _User9::write_result_txt(os);
    }

    void write_result_bin(std::ostream& os) {
      this -> _User0::write_result_bin(os);
      this -> _User1::write_result_bin(os);
      this -> _User2::write_result_bin(os);
      this -> _User3::write_result_bin(os);
      this -> _User4::write_result_bin(os);
      this -> _User5::write_result_bin(os);
      this -> _User6::write_result_bin(os);
      this -> _User7::write_result_bin(os);
      this -> _User8::write_result_bin(os);
      this -> _User9::write_result_bin(os);
    }
  };


}  // namespace nlo

#endif
