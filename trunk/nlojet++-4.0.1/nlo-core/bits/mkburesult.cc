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

  if(argc > 1) size = atoi(argv[1]);

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
  line("//  This header file was generated automatically by mk-basic_user_result-header.");
  line("#ifndef __NLO_NLO_BASIC_USER_RESULT_H__");
  line("#define __NLO_NLO_BASIC_USER_RESULT_H__ 1");
  line(""); 
  line("//   stdc++ includes");
  line("#include <string>");
  line("#include <iostream>");
  line("");
  line("");
  line("namespace nlo {");
  line("");
  line("  class _NullDistBook{};");
  line("");

  cout<<"  template<";
  for(int i = 0; i < size-1; i++) {
    cout<<"class _DistBook"<<i<<"=_NullDistBook, ";
    if(i%2 == 1)cout<<endl<<"           ";
  }
  cout<<"class _DistBook"<<size-1<<"=_NullDistBook>\n";    
  line("  class basic_user_result;");
  line("");
  line("");

  //   MAIN LOOP
  for(int l = 0; l < size; l++) {


    cout<<"  template<";
    for(int i = 0; i < l; i++) {
      cout<<"class _DistBook"<<i<<", ";
      if((i+1)%3 == 0)cout<<endl<<"           ";
    }
    cout<<"class _DistBook"<<l<<">\n";    


    cout<<"  class basic_user_result";
    if(l < size-1) {
      cout<<"<";
      for(int i = 0; i < size-1; i++) {
	if(i <= l) cout<<"_DistBook"<<i<<", ";
	else cout<<"_NullDistBook, ";
	if((i+1)%3 == 0)cout<<endl<<"                          ";
      }
      if(size-1==l) cout<<"_DistBook"<<size-1<<">\n";
      else cout<<"_NullDistBook>";
    }
    
    line("\n  {");
    line("  public:");
    line("    //   public types");
    for(int i = 0; i <= l; i++)
      cout<<"    typedef _DistBook"<<i<<" distbook"<<i<<"_type;\n";

    line("");
    line("    //   computed assignments");
    line("    basic_user_result& operator+=(const basic_user_result& x)"); 
    line("    {");
    cout<<"       ";
    for(int i = 0; i < l; i++) {
      cout<<"book"<<i<<" += x.book"<<i<<"; ";
      if((i+1)%3 == 0)cout<<endl<<"       ";
    }
    cout<<"book"<<l<<" += x.book"<<l<<";\n";
    line("       return *this;"); 
    line("    }\n");


    line("    basic_user_result& assadd_weighted(const basic_user_result& x)"); 
    line("    {");
    for(int i = 0; i <= l; i++)
      cout<<"       book"<<i<<".assadd_weighted(x.book"<<i<<");\n";
    line("       return *this;"); 
    line("    }\n");
   
    
    line("    //   public acces to the distbook objects");
    for(int i = 0; i <= l; i++)
      cout<<"    distbook"<<i<<"_type book"<<i<<";\n";   
    line("  };");
    line("");
    line("");

    cout<<"  template<";
    for(int i = 0; i < l; i++) {
      cout<<"class _DistBook"<<i<<", ";
      if((i+1)%3 == 0)cout<<endl<<"           ";
    }
    cout<<"class _DistBook"<<l<<">\n";    

    line("  unsigned long read_basic_user_output(std::istream& is,"); 
    cout<<"    basic_user_result<";
    for(int i = 0; i < size-1; i++) {
      if(i <= l) cout<<"_DistBook"<<i<<", ";
      else cout<<"_NullDistBook, ";
      if((i+1)%3 == 0)cout<<endl<<"                   ";
    }
    if(size-1==l) cout<<"_DistBook"<<size-1<<">& x)\n";
    else cout<<"_NullDistBook>& x)\n";

  line("  {"); 
  line("    bool txt, flag = false;"); 
  line("    unsigned long int nsaved;"); 
  line("    std::basic_string<char> distbook;"); 
  line("    extern void __helper_basic_user_set_warning();");     
  line("");     
  line("    is>>txt;"); 
  line("    is>>nsaved;"); 
  line("");     
  
  
  for(int i = 0; i <= l; i++) {
    cout<<"    is>>distbook;";
    cout<<" if(typeid(_DistBook"<<i<<").name() != distbook) flag = true;\n"; 
  }

  line("");        
  line("    if(flag) __helper_basic_user_set_warning();");     
  line("");        
  line("    if(txt) {");
  cout<<"       ";
  for(int i = 0; i < l; i++) {
    cout<<"is>>x.book"<<i<<"; ";
    if((i+1)%5 == 0)cout<<endl<<"       ";
  }
  cout<<"is>>x.book"<<l<<";\n";
  line("    } else {");
  cout<<"       ";
  for(int i = 0; i < l; i++) {
    cout<<"read(is, x.book"<<i<<"); ";
    if((i+1)%3 == 0)cout<<endl<<"       ";
  }
  cout<<"read(is, x.book"<<l<<");\n";
  line("    }");
  line("");        
  line("    return nsaved;"); 
  line("  }");
 
  line("");
  if(l < size-1) {
    line("");
    line("");
  }
  }  //  MAIN LOOP

  line("}  // namespace nlo");
  line("");
  line("#endif");


  return 0;
}


