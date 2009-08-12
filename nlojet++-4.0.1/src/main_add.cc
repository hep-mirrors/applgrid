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
#include <fstream>
#include <sstream>
#include <string>
#include <list>

#include <getopt.h>
#include <sys/stat.h>

#include "distribution.h"
#include "bits/nlo-basic_user_result.h"
#include "nlo++-module_add.h"

#include "ltdl.h"
#include "nlojet++.h"


//----- used namespaces -----
using namespace std;
using namespace nlo;


//----- some shorthand notation -----
typedef basic_string<char> string_type;
typedef list<string_type> string_list;

typedef void (*inputfunc_type)(unsigned int&, unsigned int&, unsigned int&);
typedef unsigned long (*module_type)(bool, const string_list&, const string_type&);


//----- generate the list of the file names ----
void file_list(unsigned int nj, const char *proc, const char *contr, 
	       const string_list& names, string_list& filelist)
{
  //----- generate the template file name ----
  struct stat  buf;
  stringstream strs;
  strs<<"@NAME@-"<<proc<<"-"<<contr<<"-"<<nj<<"jet"<<'\0';

  //----- initialize the file list -----
  filelist.clear();
    
  //---- chechk the files and fill the list -----
  string_type file;
  string_list::const_iterator in = names.begin(), end = names.end();
    
  while(in != end) {
    //----- generate the file name with the current name -----
    file = strs.str(); file.replace(file.find("@NAME@"), 6, *(in++));
    
    //----- check the file -----
    if(stat(file.c_str(), &buf) != -1) filelist.push_back(file);
  }
}


//----- print out the help -----
void main_add_help()
{
  cout
    <<"Usage: nlojet --add [-d dir] [-r dir] name1 [name2] [name3] ...\n\n"
    <<"       nlojet --add -h|--help \n\n"
    
    <<"Options:\n" 
    <<"  -h, --help     print this message\n"
    <<"  -d directory   directory where search the data files (default: ./output)\n"
    <<"  -r directory   save the result to this directory (default: ./result)\n"
    <<endl;
}



int main_add(int argc, char **argv)
{
  const option long_opt[] = {
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  //----- defoult parameters -----
  bool weighted = false;
  int njet = -1, from=-1, to = -2;
  const char * proc = 0, *contr = 0, *user = 0;
  string_type dir("./output/"), res("./result/");
  
  //----- start the parsing -----
  while(1) {
    int opt_idx = 0;      
    int c = getopt_long(argc, argv, "hwd:r:f:t:j:p:c:u:", long_opt, &opt_idx);
    if(c == -1) break;
    
    switch(c) {
    case 'h': main_add_help(); exit(0); break;
    case 'w': weighted = true; break;
    case 'd': dir = optarg; break;
    case 'r': res = optarg; break;
    case 'f': from = atoi(optarg); break;
    case 't': to = atoi(optarg); break;
    case 'j': njet = atoi(optarg); break;
    case 'p': proc = optarg; break;
    case 'c': contr = optarg; break;
    case 'u': user = optarg; break;
    case '?': break;
    default: throw; break;
    }
  }
  
  if(dir[dir.size()-1] != '/') dir += '/';
  if(res[res.size()-1] != '/') res += '/';
  
  //----- generate the names -----
  string_list names;
  
  int idx = optind;
  while(idx < argc) {
    stringstream strs;
    strs<<dir<<argv[idx++];
    names.push_back(strs.str());
  }
  
  for(int i = from; i <= to; i++) {
    stringstream strs; 
    strs<<dir<<i;
    names.push_back(strs.str());
  }
  
  //----- default adding module -----
  typedef distbook<double, void>   _Distbook0D;
  typedef distbook<double, double> _Distbook1D;
  typedef distbook<double, histpoint1d> _Distbook1H;
  typedef distbook<double, histpoint2d> _Distbook2H;

  module_type add_module = main_module_add<basic_user_result<_Distbook0D, _Distbook1D, _Distbook1H, _Distbook2H> >;

  //----- using external module instead of the default -----
  lt_dlhandle handle = NULL;
  
  if(user) {
    //----- open the user defined file -----
    if(lt_dlinit() != 0) throw lt_dlerror();
    
    if(!(handle = lt_dlopen(user))) {
      cerr<<"can't open the module "<<user<<"!\n"<<lt_dlerror()<<endl;
      
      lt_dlclose(handle); 
      lt_dlexit();
      handle = NULL;
      goto no_external_module;
    }
    
    const lt_dlinfo *info = lt_dlgetinfo(handle);
    if (!info) {
      cerr<<"can't get module info: "<<lt_dlerror()<<endl;

      lt_dlclose(handle);
      lt_dlexit();
      handle = NULL;
      goto no_external_module;
    }
    
    if(info->name) cout<<"module name: "<<info->name<<"\n";
    else cout<<"module is not a libtool module\n";
    cout<<"module filename: "<<info->filename<<"\n";
    cout<<"module reference count: "<<info->ref_count<<endl;
    
    //----- try to find the symbol table -----
    cout<<"\nTrying to find the symbol table :\n"
	<<"user_defined_functions[]        : ";
    nlo_symbol *symtable = (nlo_symbol *) lt_dlsym(handle, "user_defined_functions");
    
    if(symtable) {
      cout<<"     OK\n"<<endl;

      //----- look for adding module -----
      void *extmod = find_symbol(symtable, "main_module_add");
      if(extmod) add_module = (module_type) extmod;

      //----- set the process id -----
      proc = (const char *) find_symbol(symtable, "procindex");
      
      //----- set the number of the jets ------
      unsigned int nu = 2U, nd = 3U, nj =  njet >= 0 ? njet : 0;
      (*((inputfunc_type) find_symbol(symtable, "inputfunc")))(nj, nu, nd);
      njet = nj;
    } else {
      cout<<"     FAIL\n"<<endl;
      
      lt_dlclose(handle);
      lt_dlexit();
      handle = NULL;
    }
  }
  
 no_external_module:
  
  for(unsigned int ip = 0; proctbl[ip].name; ip++)
    if(proc == 0 || strcmp(proctbl[ip].name, proc) == 0) {
      string_type rdir_p = res + proctbl[ip].name;
      
      for(unsigned int nj = 0; proctbl[nj].njet[nj] != -1; nj++)
		if(proctbl[nj].njet && ((int) nj == njet || njet == -1)) {
		  bool first = true;

		  std::stringstream strs; strs<<"/"<<nj<<"-jet/";
		  string_type rdir_j = rdir_p + strs.str();
	    
		  for(unsigned int ic = 0; contbl[ic]; ic++)
			if(contr == 0 || strcmp(contbl[ic], contr) == 0) {
			  //----- make the file list -----
			  string_list filelist;
			  file_list(nj, proctbl[ip].name, contbl[ic], names, filelist);
	      
			  if(filelist.size() > 0) {
				//----- collect the results write out to the disk ---
				string_type rdir = rdir_j + contbl[ic]; 
		
				make_dir(rdir.c_str());
				unsigned long nsaved = (*add_module)(weighted, filelist, rdir);
		
				if(nsaved > 0UL) {
				  if(first) {
					cout<<proctbl[ip].desc<<"  "<<nj<<"-jet : "<<endl;
					first = false;
				  }
		  
				  cout<<"   "<<contbl[ic]<<" : "<<nsaved<<endl;
				}
			  }
			}
		}
    }
  
  if(handle) {
    lt_dlclose(handle);
    lt_dlexit();
  }
  
  return 0;
}
