// emacs: this is -*- c++ -*-
//
//   tsparse2d.h        
//
//   very basic 2d sparse matrix class, needs to be improved, 
//   and extended to 3d.
//
//   NB: be VERY CAREFUL about the use of the (i,j) index operator, 
//       since it checks to see whether the given element exists,
//       and if not, it creates it, so even viewing it such as...
//
//          tsparse2d<...> d(10,10); cout << a(2,3) << endl;                    
//
//       will create the element if it doesn't exist *unless* you 
//       are looking at a const tsparse2d<...>
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: tsparse2d.h, v1.0   Fri Nov  9 00:17:31 CET 2007 sutt


#ifndef __TSPARSE2D_H
#define __TSPARSE2D_H

#include <iostream>
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

#include <iomanip>
using std::setw;
using std::setprecision;

#include <cmath>
using std::fabs;
using std::pow;
using std::log10;



#include "appl_grid/tsparse1d.h"


template<class T> 
class tsparse2d : public tsparse_base {

public:

  tsparse2d(int nx, int ny) : tsparse_base(nx), m_Ny(ny), m_v(NULL) {  
    m_v = new tsparse1d<T>*[m_Nx];
    for ( int i=0 ; i<Nx() ; i++ ) m_v[i] = new tsparse1d<T>(m_Ny);
  }

  tsparse2d(int nx, int ny, int lx, int ly, int _size=1, int _sizey=1) : tsparse_base(nx, lx, _size), m_Ny(ny), m_v(NULL) {  
    if ( _size ) m_v = new tsparse1d<T>*[_size];
    for ( int i=0 ; i<_size ; i++ ) m_v[i] = new tsparse1d<T>(m_Ny,ly, _sizey);
  }


  
  ~tsparse2d() { 
    if ( m_v ) { 
      for ( int i=0 ; i<m_ux-m_lx+1 ; i++ )  if ( m_v[i] ) delete m_v[i];
      delete[] m_v;
    }
  }
  
  



  void trim() { 
    
    // trim each column down 
    for ( int i=0 ; i<m_ux-m_lx+1 ; i++ ) if ( m_v[i] ) m_v[i]->trim();

    int xmin = m_lx;
    int xmax = m_ux;
    
     // now find new limits
    for ( ; xmin<xmax+1 && (*m_v[xmin-m_lx])==0 ; xmin++ )  delete m_v[xmin-m_lx]; 
    for ( ; xmin<xmax   && (*m_v[xmax-m_lx])==0 ; xmax-- )  delete m_v[xmax-m_lx];

    // don't need to change anything
    if ( xmin==m_lx && xmax==m_ux ) return;
 
     // trimmed matrix is empty
    if ( xmax<xmin ) { 
      m_lx = xmin;  
      m_ux = xmax;
      delete[] m_v; 
      m_v = NULL;
      return; 
    }

    // copy to trimmed matrix
    tsparse1d<T>** m_vnew = new tsparse1d<T>*[xmax-xmin+1]; 
    tsparse1d<T>** mnp = m_vnew;
    tsparse1d<T>** mvp = m_v+xmin-m_lx;
    for ( int j=xmin ; j<=xmax ; j++ )  (*mnp++) = (*mvp++);
    delete[] m_v;

    // update the sparse data and corresponding limits 
    m_lx = xmin;
    m_ux = xmax;
    m_v  = m_vnew;
  }
  

  const tsparse1d<T>* operator[](int i) const {
    //    range_check(i);
    if ( i<m_lx || i>m_ux ) return NULL;
    return m_v[i-m_lx];
  } 

  void untrim() { 
    grow(0);
    grow(Nx()-1);
    tsparse1d<T>** mvp = m_v;
    for ( int i=m_lx ; i<=m_ux ; i++ ) (*mvp++)->untrim();
  }


  
  // NB: for some reason, even though this is a const function, 
  //     it allows non const member functions to be called on
  //     non const members so we cast the member to const to enforce
  //     constness.
  T operator()(int i, int j) const {
    //    range_check(i);
    if ( i<m_lx || i>m_ux ) return 0;
    return (*((const tsparse1d<T>**)m_v)[i-m_lx])(j);
  } 


  T& operator()(int i, int j) {
    //  range_check(i);
    grow(i);
    if ( m_v[i-m_lx] ) return (*m_v[i-m_lx])(j);
  } 


  int Ny() const { return m_Ny; } 


  bool trimmed(int i, int j) const { 
    if ( i<m_lx || i>m_ux ) return false;
    else                    return m_v[i-m_lx]->trimmed(j);
  } 

  tsparse1d<T>* trimmed(int i) const { 
    if ( i<m_lx || i>m_ux ) return NULL;
    else                    return m_v[i-m_lx];
  } 


  int size() const { 
    int N=0;
    for ( int i=0 ; i<m_ux-m_lx+1 ; i++ )  if ( m_v[i] ) N += m_v[i]->size();
    return N; //+3*sizeof(int;
  }
 

  int ymin() {
    int minx = m_Ny;
    tsparse1d<T>** mvp = m_v;
    for ( int i=m_lx ; i<=m_ux ; i++ ) { 
      int _minx=(*mvp++)->xmin();
      if ( _minx<=minx ) minx=_minx;
    }
    return minx;
  }


  int ymax() {
    int maxx = -1;
    tsparse1d<T>** mvp = m_v;
    for ( int i=m_lx ; i<=m_ux ; i++ ) { 
      int _minx=(*mvp)->xmin(); 
      int _maxx=(*mvp++)->xmax(); 
      if ( _minx<=_maxx &&_maxx>maxx ) maxx=_maxx;
    }
    if ( maxx==-1 ) maxx=m_Ny-1;
    return maxx;
  }


  void print() const { 
    
    //   cout << "\t\t+";
    //   for ( int j=0 ; j<Ny() ; j++ ) cout << "-";
    //   cout << "+\n" << endl; 
   
    cout << "\t\t+";
    for ( int j=0 ; j<Ny() ; j++ ) cout << "--";
    cout << "+\n";

    for ( int i=0 ; i<Nx() ; i++ ) {
      cout << setw(9) << trimmed(i) << "\t";
      
      cout << "|";

      for ( int j=0 ; j<Ny() ; j++ ) {
#if 1
      	if ( trimmed(i,j) ) { 
	  //	  cout << setprecision(2) << setw(4) << mant((*this)(i,j)) << "\t";
	  cout << "oo";
	}
      	else  if ( i==j ) cout << ". "; 
	else              cout << "  "; 
#else
      	if ( trimmed(i,j) ) { 
	  cout << setprecision(1) << setw(1) << mant(fabs((*this)(i,j)));
	}
      	else  if ( i==j ) cout << "  "; 
	else              cout << "  "; 
#endif
      }
      //      cout << "|\n";
      cout << "|\n";
    }

#if 0
    ////////////////////////////////////////////////////////////////////
    
    cout << "\nPAW null 0 " << Nx() <<  "0 " << Ny() << endl;

    for ( int i=0 ; i<Nx() ; i++ ) {
      cout << setw(9) << trimmed(i) << "\t";
      
      for ( int j=0 ; j<Ny() ; j++ ) {

      	if ( trimmed(i,j) ) { 
	  //	  cout << setprecision(2) << setw(4) << mant((*this)(i,j)) << "\t";
	  cout << "\nPAW box " << i-1 << " " << i << " " << j-1 << " " << j << endl;
	}
	//      	else  if ( i==j ) cout << ". "; 
	//	else              cout << "  "; 
      }
      //      cout << "|\n";
      // cout << "|\n";
    }
    //////////////////////////////////////////////////////////////
    cout << "\nPAW wait" << endl;

#endif

    cout << "\t\t+";
    for ( int j=0 ; j<Ny() ; j++ ) cout << "--";
    cout << "+\n";
   
    //   cout << "\t\t+";
    //   for ( int j=0 ; j<Ny() ; j++ ) cout << "-";
    //   cout << "+\n" << endl; 
    
  }

  
  tsparse1d<T>**  v() { return m_v; }


  tsparse2d& operator*=(const double& d) { 
    for ( int i=0 ; i<Nx() ; i++ ) {     
      if ( tsparse_base::trimmed(i) ) { 
	if ( m_v[i-m_lx] ) (*m_v[i-m_lx])*=d;
      }
    }
    return *this;
  }
  

  tsparse2d& operator+=(const tsparse2d& t) {
    if ( Nx()!=t.Nx() || Ny()!=t.Ny() ) throw out_of_range("bin mismatch");
    for ( int i=0 ; i<Nx() ; i++ ) {
      for ( int j=0 ; i<Ny() ; j++ ) (*this)(i,j) += t(i,j);
    }
    return *this;
  }

private:


 
  void grow(int i) { 
     
    // nothing needs to be done; 
    if ( i>=m_lx && i<=m_ux ) return;

    // is it empty? add a single row
    if ( m_lx>m_ux ) { 
      m_v    = new tsparse1d<T>*[1];
      (*m_v) = new tsparse1d<T>(m_Ny, m_Ny, 0);
      m_lx = m_ux = i;
      return;
    }

    tsparse1d<T>** newv = new tsparse1d<T>*[ ( i<m_lx ? m_ux-i+1 : i-m_lx+1 ) ];
    tsparse1d<T>** newp = newv;
    tsparse1d<T>** mvp  = m_v;

    int lx = m_lx;
    int ux = m_ux;

    for ( ; m_lx>i ; m_lx-- )  (*newp++) = new tsparse1d<T>(m_Ny, m_Ny, 0);
    for ( ; lx<=ux ; lx++ )    (*newp++) = (*mvp++);
    for ( ; m_ux<i ; m_ux++ )  (*newp++) = new tsparse1d<T>(m_Ny, m_Ny, 0);

    delete[] m_v;
    m_v = newv;
  }


private:
    
  int m_Ny;

  tsparse1d<T>**  m_v;

};


// stream IO template
template<class T>ostream& operator<<(ostream& s, const tsparse2d<T>& sp) { 
  for ( int i=0 ; i<sp.Nx() ; i++ ) { 
    for ( int j=0 ; j<sp.Ny() ; j++ ) { 
      s << sp(i,j) << "\t"; 
    }
    s << "\n";
  }
  return s;
}


#endif  // __TSPARSE2D_H 










