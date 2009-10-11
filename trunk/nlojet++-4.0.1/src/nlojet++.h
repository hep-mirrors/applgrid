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
#ifndef __NLO_NLOJETXX_H__
#define __NLO_NLOJETXX_H__ 1


#include <climits>
#include <string>

struct main_input;

//----- process list -----
struct process_table {
  const char * name;
  const char * desc;
  const int njet[20];
  void (*module_calc)(const main_input&);
};

extern const process_table proctbl[];
  

//---- table of the contributions -----
extern const char *contbl[];

//----- make directory function -----
int make_dir(const char *);

  
//   find the user defined functions
struct nlo_symbol {
  const char *name;
  void *address;
};
  
void * find_symbol(const nlo_symbol *, const char *);


//      class for storing the input parameters
struct main_input {
  //  type of the contribution (born, nlo, scale)
  int contr;
  
  //  number of the event
  unsigned long nevent;
  
  //  save after nsave points
  unsigned long nsave;

  //  nonphysical cut parameter in the dipole terms
  double alpha;

  //  Monte Carlo helicity sum
  bool mchel;
  
  //  Time rate between the n+1 and n parton integral
  unsigned int time;

  //  true if the output is formated output
  bool txtout;

  //  user defined functions
  nlo_symbol *symtable;
  
  //  file names 
  mutable std::basic_string<char> outfile;

  //  default constructor
  main_input() 
    : contr(-1), nevent(ULONG_MAX), nsave(10000UL),
      alpha(0.1), mchel(true), time(1U), txtout(false) {}
  
  //----- find a symbol in the symbol table ----
  void * find_symbol(const char *symbol) const {
    return ::find_symbol(symtable, symbol);
  }
};

#endif
