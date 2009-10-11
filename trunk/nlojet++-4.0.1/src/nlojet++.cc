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
#include <unistd.h>
#include <sys/wait.h>
#include <cstring>
#include "nlojet++.h"


using namespace std;

extern void main_module_epa(const main_input&);
extern void main_module_dis(const main_input&);
extern void main_module_hhc(const main_input&);
extern void main_module_hhc2ph(const main_input&);
extern void main_module_photo(const main_input&);

//------ process table -----
const process_table proctbl[] = {
  {"epa",      "e+e- annihilation",					       {0, 0, 0, 1, 1, 1,-1}, main_module_epa},
  {"dis",      "deeply inelastic scatering",		       {0, 0, 1, 1, 1,-1}   , main_module_dis},
  {"hhc",      "hadron-hadron collision", 		           {0, 1, 1, 1, 1,-1}   , main_module_hhc},
  {"hhc2ph",   "hadron-hadron collision with two photons", {0, 1,-1}            , main_module_hhc2ph},
  {"photodir", "photoproduction (direct photon)", 	       {0, 1, 1, 1, 1,-1}   , main_module_photo},
  {"photores", "photoproduction (resolved photon)",        {0, 1, 1, 1, 1,-1}   , main_module_hhc},
  {0,0,{-1}, 0}
};

//----- contribution types ------
const char *contbl[] = {"born", "nlo", "full", 0};


int make_dir(const char *dir)
{
  int status;
  pid_t cpid = fork();
  switch (cpid) {
  case -1:    //  Can't fork.
    return -1;
  case 0: 
    execlp("mkdir", "mkdir", "-p", dir, (char *) 0); 
    _exit(1);
  default:
    while(wait(&status) != cpid);
    if(status) return -1;  //  mkdir failed
    return 0;
  }
}


void * find_symbol(const nlo_symbol *handle, const char *symbol)
{
  while(handle->name) {
    if(std::strcmp(symbol, handle->name) == 0)
      return handle -> address;
    ++handle;
  }
  
  return 0;
}
