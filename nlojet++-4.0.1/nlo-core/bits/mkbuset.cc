#include <iostream>
#include <cstdlib>
#ifndef SIZE
#define SIZE 10
#endif


using namespace std;


void line(const char *str)
{
  cout<<str<<endl;
}


int main(int argc, const char **argv) 
{
  int size = SIZE;

  if(argc > 1) size =  atoi(argv[1]);

  line("//  Copyright (C) 2004 Zoltan Nagy");
  line("//");
  line("//  This program is free software; you can redistribute it and/or modify");
  line("//  it under the terms of the GNU General Public License as published by");
  line("//  the Free Software Foundation; either version 2 of the License, or");
  line("//  (at your option) any later version.");
  line("//");
  line("//  This program is distributed in the hope that it will be useful,");
  line("//  but WITHOUT ANY WARRANTY; without even the implied warranty of");
  line("//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the");
  line("//  GNU General Public License for more details.");
  line("//");
  line("//  You should have received a copy of the GNU General Public License");
  line("//  along with this program; if not, write to the Free Software");
  line("//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA");
  line(""); 
  line("//  This header file was generated automatically by mk-basic_user_set-header.");
  line("#ifndef __NLO_NLO_BASIC_USER_SET_H__");
  line("#define __NLO_NLO_BASIC_USER_SET_H__ 1");
  line(""); 
  line("//   nlojet++ includes");
  line("#include <bits/nlo-basic_user.h>");
  line("");
  line("");
  line("namespace nlo {");
  line("");
  line("  class _NullUser{};");
  line("");

  cout<<"  template<";
  for(int i = 0; i < size-1; i++) {
    cout<<"class _User"<<i<<"=_NullUser, ";
    if(i%2 == 1)cout<<endl<<"           ";
  }
  cout<<"class _User"<<size-1<<"=_NullUser>\n";    
  line("  class basic_user_set;");
  line("");
  line("");

  //   MAIN LOOP
  for(int l = 0; l < size; l++) {
    cout<<"  template<";
    for(int i = 0; i < l; i++) {
      cout<<"class _User"<<i<<", ";
      if((i+1)%4 == 0)cout<<endl<<"           ";
    }
    cout<<"class _User"<<l<<">\n";    


    cout<<"  class basic_user_set";
    if(l < size-1) {
      cout<<"<";
      for(int i = 0; i < size-1; i++) {
	if(i <= l) cout<<"_User"<<i<<", ";
	else cout<<"_NullUser, ";
	if((i+1)%5 == 0)cout<<endl<<"                       ";
      }
      if(size-1==l) cout<<"_User"<<size-1<<">\n";
      else cout<<"_NullUser>";
    }
    
    cout<<"\n    : ";
    for(int i = 0; i < l; i++) {
      cout<<"public _User"<<i<<", ";
      if((i+1)%4 == 0)cout<<endl<<"      ";
    }
    cout<<"public _User"<<l<<"\n";    
   
    line("  {");
    line("  public:");
    line("    //   public types");

    for(int i = 0; i <= l; i++)
      cout<<"    typedef _User"<<i<<" user"<<i<<"_type;\n";

    line("");
    line("  private:");
    line("    //   do any operetaions at the end of the event");
    line("    void operations_at_the_end_of_event() {");
    for(int i = 0; i <= l; i++)
      cout<<"      this ->  _User"<<i<<"::operations_at_the_end_of_event();\n";
    line("    }");
    
    line("");
    line("    //   write some information about the type of the sample");
    line("    void write_typeinfo(std::ostream& os) {");
    for(int i = 0; i <= l; i++)
      cout<<"      this -> _User"<<i<<"::write_typeinfo(os);\n";
    line("    }");
    
    line("");
    line("    //   save the results in txt and binary mode");
    line("    void write_result_txt(std::ostream& os) {");
    for(int i = 0; i <= l; i++)
      cout<<"      this -> _User"<<i<<"::write_result_txt(os);\n";
    line("    }");
    
    line("");
    line("    void write_result_bin(std::ostream& os) {");
    for(int i = 0; i <= l; i++)
      cout<<"      this -> _User"<<i<<"::write_result_bin(os);\n";
    line("    }"); 
    line("  };");
    
    if(l < size-1) {
      line("");
      line("");
      line("");
    }
  }  //  MAIN LOOP
  
  line("");
  line("");
  line("}  // namespace nlo");
  line("");
  line("#endif");


  return 0;
}



    
