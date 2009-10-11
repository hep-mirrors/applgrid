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
#include <cstdlib>
#include <iostream>
#include <getopt.h>


//----- the used namespaces -----
using namespace std;


extern void main_calc_help();
extern void main_add_help();

extern int main_add(int, char **);
extern int main_calc(int, char **);



int main(int argc, char **argv) 
{
  option longopt[] = {
    {"add", no_argument, 0, 0},
    {"calculate", no_argument, 0, 0},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  try {  
    while(1) {
      int optidx = 0;      
      int c = getopt_long(argc, argv, "h", longopt, &optidx);
      if(c == -1) break;
      
      switch(c) {
    case 'h': 
      main_calc_help();
      main_add_help();
      exit(0);
      break;
      case 0:
	switch(optidx) {
	case 0: main_add(argc, argv); break;
	case 1: main_calc(argc, argv); break;
      }
	break;
      case '?': break;
      default: throw; break;
      }
    }
  } catch(const char *message) {
    cerr<<"Error : "<<message<<endl;
    return 1;
  } catch(...) {
    cerr<<"Unexpected error!"<<endl;
    return 1;
  }

  return 0;
}
