//
//   @file    combine.cxx         
//            simple code to add together many grids 
//
//   @author M.Sutton
// 
//   Copyright (C) 2013 M.Sutton (sutt@cern.ch)    
//
//   $Id: combine.cxx, v0.0   Sat 13 Jul 2013 09:54:51 CEST sutt $


#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>

#include "appl_grid/appl_grid.h"
#include "amconfig.h"


int usage(std::ostream& s, int argc, char** argv) { 
  s << "Usage: " << argv[0] << " [OPTIONS] -o output_grid.root  input_grid.root [input_grid1.root ... input_gridN.root]\n\n";
  s << "  APPLgrid \'" << argv[0] << "\' adds " << PACKAGE_STRING << " grid files together into a single grid\n\n"; 
  s << "Configuration: \n";
  s << "    -o filename   \tname of output grid (filename required)\n\n";
  s << "Options: \n";
  s << "    -s, --scale value\tscale output grid by value, \n";
  s << "        --verbose \tdisplay grid documentation during add\n";     
  s << "    -v, --version \tdisplays the APPLgrid version\n";
  s << "    -h, --help    \tdisplay this help\n";
  s << "\nSee " << PACKAGE_URL << " for more details\n"; 
  s << "\nReport bugs to <" << PACKAGE_BUGREPORT << ">";
  s << std::endl;
  return 0;
}



int main(int argc, char** argv) { 

  /// check correct number ofg parameters
  if ( argc<2 ) return usage( std::cerr, argc, argv );


  std::string output_grid = "";

  /// read in the parameters 
  std::vector<std::string> grids;
  for ( int i=1 ; i<argc ; i++ ) { 
    if ( std::string(argv[i])=="-h" || std::string(argv[i])=="--help" )    return usage( std::cout, argc, argv ); 
    if ( std::string(argv[i])=="-v" || std::string(argv[i])=="--version" ) {
      std::cout << argv[0] << " APPLgrid version " << PACKAGE_VERSION << std::endl; 
      return 0;
    }
  }

  double d = 1;

  bool verbose = false; 

  for ( int i=1 ; i<argc ; i++ ) { 
    if      ( std::string(argv[i])=="--verbose" ) verbose = true;
    else if ( std::string(argv[i])=="-o" ) { 
      ++i;
      if ( i<argc ) output_grid = argv[i];
      else  return usage( std::cerr, argc, argv );
    }
    else if ( std::string(argv[i])=="-s" || std::string(argv[i])=="--scale" ) { 
      ++i;
      if ( i<argc ) d = std::atof(argv[i]);
      else  return usage( std::cerr, argc, argv );
    }
    else { 
      grids.push_back(argv[i]);
    }
  }

  if ( grids.size()<1 ) return usage(std::cerr, argc, argv);

  if ( output_grid=="" ) return usage(std::cerr, argc, argv);

  /// now add the grids together

  
  appl::grid  g( grids[0] );

  g.untrim();

  if ( verbose ) std::cout << g.getDocumentation() << std::endl;

  for ( unsigned i=1 ; i<grids.size() ; i++ ) { 
    appl::grid  _g( grids[i] );
    _g.untrim();
    if ( verbose ) std::cout << _g.getDocumentation() << std::endl;
    g += _g;
  }

  if ( d!=1 ) g *= d;

  std::cout << "writing " << output_grid << std::endl;
  g.Write(output_grid);

  std::cout << argv[0] << ": added " << grids.size() << " grids" << std::endl; 
  std::cout << argv[0] << ": output to  " << output_grid << std::endl; 

  return 0;
}
