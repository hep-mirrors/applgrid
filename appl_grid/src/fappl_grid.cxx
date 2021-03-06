//   
//   fappl_grid.cxx        
//     
//      fortran callable wrapper functions for the c++  
//      appl grid project.
//                   
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: fappl.cxx, v1.0   Wed May 21 14:31:36 CEST 2008 sutt

#include <map>
#include <iostream>
using std::cout;
using std::endl;

#include "appl_grid/appl_grid.h"
#include "appl_grid/fastnlo.h"



/// externally defined alpha_s and pdf routines for fortran 
/// callable convolution wrapper
extern "C" double fnalphas_(const double& Q); 
extern "C" void   fnpdf_(const double& x, const double& Q, double* f);

/// create a grid
extern "C" void bookgrid_(int& id, const int& Nobs, const double* binlims);

/// delete a grid
extern "C" void releasegrid_(int& id);

/// delete all grids
extern "C" void releasegrids_();

/// read a grid from a file
extern "C" void readgrid_(int& id, const char* s);

/// write to a file 
extern "C" void writegrid_(int& id, const char* s);

/// add an entry 
extern "C" void fillgrid_(int& id, 
			  const int& ix1, const int& ix2, const int& iQ, 
			  const int& iobs, 
			  const double* w,
			  const int& iorder );  

/// redefine the grid dimensions
extern "C" void redefine_(int& id, 
			  const int& iobs, const int& iorder, 
			  const int& NQ2, const double& Q2min, const double& Q2max, 
			  const int& Nx,  const double&  xmin, const double&  xmax); 

/// get number of observable bins for a grid 
extern "C" int getnbins_(int& id);

/// do the convolution!! hooray!!
extern "C" void convolute_(int& id, double* data);
extern "C" void convoluteorder_(int& id, int& nloops, double* data);

extern "C" void convolutewrap_(int& id, double* data, 
			       void (*pdf)(const double& , const double&, double* ),
			       double (*alphas)(const double& ) );


extern "C" void fullconvolutewrap_(int& id, double* data, 
				   void (*pdf)(const double& , const double&, double* ),
				   double (*alphas)(const double& ),
				   int nloops,
				   double rscale, double fscale  );

/// print a grid
extern "C" void printgrid_(int& id);

/// print all grids
extern "C" void printgrids_();

/// print the grid documentation
extern "C" void printgriddoc_(int& id);

/// create grids from fastnlo
extern "C" void readfastnlogrids_(  int* ids, const char* s );


static int idcounter = 0;
static std::map<int,appl::grid*> _grid;


/// grid map management

extern "C" void ngrids_(int& n) { n=_grid.size(); }

extern "C" void gridids_(int* ids) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.begin();
  for ( int i=0 ; gitr!=_grid.end() ; gitr++, i++ ) ids[i] = gitr->first;
}



void bookgrid_(int& id, const int& Nobs, const double* binlims) 
{
  id = idcounter++;

  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);

  if ( gitr==_grid.end() ) {
    cout << "bookgrid_() creating grid with id " << id << endl; 
    _grid.insert(  std::map<int,appl::grid*>::value_type( id, new appl::grid( Nobs, binlims,
									      2,    10, 1000, 1,
									      12,  1e-5, 1, 3, 
									      "nlojet", 1, 3, "f3") ) ) ;									 
    //  _grid->symmetrise(true);
  }
  else throw appl::grid::exception( std::cerr << "grid with id " << id << " already exists" << std::endl );  

}


void readgrid_(int& id, const char* s) {
  id = idcounter++;
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr==_grid.end() ) { 
    _grid.insert(  std::map<int,appl::grid*>::value_type( id, new appl::grid(s) ) );
  }
  else throw appl::grid::exception( std::cerr << "grid with id " << id << " already exists" << std::endl );  
}


  
void printgrid_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    std::cout << "grid id " << id << "\n" << *gitr->second << std::endl;
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void printgrids_() { 
  std::map<int,appl::grid*>::iterator gitr = _grid.begin();
  for ( ; gitr!=_grid.end() ; gitr++ ) { 
    std::cout << "grid id " << gitr->first << "\n" << *gitr->second << std::endl;
  }
}

  
void printgriddoc_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    std::cout << "grid id " << id << "\n" << gitr->second->getDocumentation() << std::endl;
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void releasegrid_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() )     { 
    delete gitr->second; 
    _grid.erase(gitr);
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void releasegrids_() { 
  std::map<int,appl::grid*>::iterator gitr = _grid.begin();
  for ( ; gitr!=_grid.end() ; gitr++ ) { 
    delete gitr->second; 
    _grid.erase(gitr);
  }
}


void redefine_(int& id, 
	       const int& iobs, const int& iorder, 
	       const int& NQ2, const double& Q2min, const double& Q2max, 
	       const int& Nx,  const double&  xmin, const double&  xmax) 
{
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) {
    gitr->second->redefine(iobs, iorder, 
			   NQ2, Q2min, Q2max, 
			   Nx,   xmin,  xmax); 
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
  
} 



int getnbins_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) return gitr->second->Nobs();
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}



void convolute_(int& id, double* data) { 
  convolutewrap_(id, data, fnpdf_, fnalphas_); 
}


void convoluteorder_(int& id, int& nloops, double* data) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    vector<double> v = g->vconvolute(fnpdf_, fnalphas_, nloops);
    for ( unsigned i=0 ; i<v.size() ; i++ ) data[i] = v[i];      
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}

void convolutewrap_(int& id, double* data, 
		       void (*pdf)(const double& , const double&, double* ),  
		       double (*alphas)(const double& ) ) {  
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    vector<double> v = g->vconvolute( pdf, alphas);
    for ( unsigned i=0 ; i<v.size() ; i++ ) data[i] = v[i];      
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}



void fullconvolutewrap_(int& id, double* data, 
			void (*pdf)(const double& , const double&, double* ),  
			double (*alphas)(const double& ),
			int nloops,
			double rscale, double fscale  ) {  
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    vector<double> v = g->vconvolute( pdf, alphas, nloops, rscale, fscale);
    for ( unsigned i=0 ; i<v.size() ; i++ ) data[i] = v[i];      
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void writegrid_(int& id, const char* s) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    cout << "writegrid_() writing " << s << "\tid " << id << endl;
    appl::grid* g = gitr->second;
    g->trim();
    //   g->print();
    g->Write(s);
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}



void fillgrid_(int& id, 
	       const int& ix1, const int& ix2, const int& iQ,  
	       const int& iobs, 
	       const double* w, 
	       const int& iorder ) { 
  //  cout << "ix " << ix1 << " " << ix2 << "  iQ" << iQ << " " << iobs << "  " << iorder << endl;  
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    gitr->second->fill_index(ix1, ix2, iQ, iobs, w, iorder);
  }  
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
 
}


void readfastnlogrids_( int* ids, const char* s ) { 

  /// create the fastnlo grids
  fastnlo f(s);

  /// don't want the grids managed by the fastnlo object, 
  /// manage them in fortran with the map
  f.manageGrids(false);

  ///copy to the fortran accessible grid map
  std::vector<appl::grid*> grids = f.grids();

  //  std::cout << "hooray!" << std::endl;
  
  for ( unsigned i=0 ; i<grids.size() ; i++ ) { 
    int id = idcounter++;
    std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
    if ( gitr==_grid.end() )  { 
      _grid.insert(  std::map<int,appl::grid*>::value_type( id, grids[i] ) );
      // std::cout << grids[i]->getDocumentation() << std::endl;
    }
    else throw appl::grid::exception( std::cerr << "grid with id " << id << " already exists" << std::endl );
    ids[i] = id;
  }  

}

  

