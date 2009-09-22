// emacs: this is -*- c++ -*-
//
//   appl_igrid.h        
//
//  grid class header - all the functions needed to create and 
//  fill the grid from an NLO calculation program, based on the 
//  class from D.Clements and T.Carli.
//
//  See contribution from T.Carli et al from the HERA-LHC 
//  workshop - working group 3 
//  
//  Copyright (C) 2007 Mark Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: appl_igrid.h, v2.00  Fri Nov  2 05:31:03 GMT 2007 sutt $

// Fixme: this needs to be tidied up. eg there are too many different, 
//        and too many version of, accessors for x/y, Q2/tau etc there 
//        should be only  one set, for x and Q2 *or* y and tau, but 
//        not both. In fact maybe they should be for x and Q2, since y 
//        and tau should perhaps be purely an internal grid issue of no 
//        concern for the user.


#ifndef __APPL_IGRID_H
#define __APPL_IGRID_H

#include <iostream>
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <string>
using std::string;

#include <cmath>
using std::abs;
using std::fabs;
using std::log;
using std::exp;
using std::sqrt;
using std::pow;

#include "appl_grid/Directory.h"
#include "appl_grid/SparseMatrix3d.h"
#include "appl_grid/appl_pdf.h"


namespace appl {



class igrid {

private:

  // grid error exception
  class exception { 
  public:
    exception(const string& s) { cerr << s << endl; }; 
    exception(ostream& s)      { cerr << s << endl; }; 
  };


  // structure to store the x<->y transform pairs for 
  struct transform_vec {
    transform_vec() : mfx(NULL), mfy(NULL) { } 
    transform_vec(double (*__fx)(double), double (*__fy)(double)) : mfx(__fx), mfy(__fy) { } 
    double (*mfx)(double);
    double (*mfy)(double);
  };
  

public:

  igrid();

  igrid(int NQ2,   double Q2min=10000.0, double Q2max=25000000.0,  int Q2order=5,  
	int Nx=50, double xmin=1e-5,     double xmax=0.9,          int xorder=5, 
	string transform="f", int Nproc=6);

  igrid(const igrid& g);

  // read grid from stored file
  igrid(TFile& f, const string& s);

  ~igrid();

  // optimise the grid dinemsions  
  void optimise() { optimise(m_Ntau, m_Ny1, m_Ny2); }
  void optimise(int NQ2, int Nx1, int Nx2);  
 
  // trim unfilled elements
  void trim() {
    for ( int i=0 ; i<m_Nproc ; i++ )  m_weight[i]->trim();    
  }

  // formatted print
  void print() const {
    header(cout);
    for ( int i=0 ; i<m_Nproc ; i++ ) { 
      cout << "sub process " << i << endl; 
      m_weight[i]->print();
    }
  }

  // return the number of words used for storage
  int size() const {
    int _size = 0;
    for ( int i=0 ; i<m_Nproc ; i++ ) _size += m_weight[i]->size();
    return _size;
  }

  // inflate unfilled elements
  void untrim() { for ( int i=0 ; i<m_Nproc ; i++ ) m_weight[i]->untrim(); }

  // write to the current root directory
  void write(const string& name);
  
  // update grid with one set of event weights
  void fill(const double x1, const double x2, const double Q2, const double* weight);

  // fast filling with no interpolation for optimisation
  void fill_phasespace(const double x1, const double x2, const double Q2, const double* weight);

  void fill_index(const int ix1, const int ix2, const int iQ2, const double* weight);

  // get the sparse structure for easier access  
  const SparseMatrix3d* weightgrid(int ip) { return m_weight[ip]; }
  SparseMatrix3d**      weightgrid()       { return m_weight; } 


  // this section stores the available x<->y transforms.
  // they are static functions and the function pairs are stored in a map with a 
  // (const string) tag - the tag can be saved to a root file to uniquely identify 
  // the transform pair.
  // additional user defined transform pairs can be added to the map using the static 
  // add_transform method before creating the grid

  // transform 
  double (*fy)(double x);
  double (*fx)(double y);

  string transform() const { return m_transform; } 

  // potential static transforms and 
  static map<const string, transform_vec> init_fmap() { 
    map<const string, transform_vec> fmap; 
    fmap["f"]  = transform_vec( _fx , _fy  );
    fmap["f0"] = transform_vec( _fx0, _fy0 );
    fmap["f1"] = transform_vec( _fx1, _fy1 );
    fmap["f2"] = transform_vec( _fx2, _fy2 );
    fmap["f3"] = transform_vec( _fx3, _fy3 );
    return fmap;
  }

  static void add_transform(const string transform, 
			    double (*__fx)(double), double (*__fy)(double)) { 
    if ( m_fmap.find(transform)!=m_fmap.end() ) { 
      throw exception("igrid::add_fmap() transform "+transform+" already in map");
    }
    m_fmap[transform] = transform_vec( __fx, __fy );
  }

  // define all these so that ymin=fy(xmin) rather than ymin=fy(xmax)
  static double _fy(double x) { return log(1/x-1); } 
  static double _fx(double y) { return 1/(1+exp(y)); } 

  static double _fy0(double x) { return -log(x); } 
  static double _fx0(double y) { return exp(-y); } 

  static double _fy1(double x) { return sqrt(-log(x)); } 
  static double _fx1(double y) { return exp(-y*y); } 

  static double _fy2(double x) { return -log(x)+m_transvar*(1-x); }
  static double _fx2(double y) {
    // use Newton-Raphson: y = ln(1/x)
    // solve   y - yp - a(1 - exp(-yp)) = 0
    // deriv:  - 1 -a exp(-yp)

    if ( m_transvar==0 )  return exp(-y); 
        
    const double eps  = 1e-12;  // our accuracy goal
    const int    imax = 100;    // for safety (avoid infinite loops)
    
    double yp = y;
    double x, delta, deriv;
    for ( int iter=imax ; iter-- ; ) {
      x = exp(-yp);
      delta = y - yp - m_transvar*(1-x);
      if ( std::fabs(delta)<eps ) return x; // we have found good solution
      deriv = -1 - m_transvar*x;
      yp  -= delta/deriv;
    }
    // exceeded maximum iterations 
    cerr << "_fx2() iteration limit reached" << endl; 
    return exp(-yp); 
  }
  
  static double _fy3(double x) { return sqrt(-log10(x)); }
  static double _fx3(double y) { return pow(10,-y*y); } 


  // pdf weight function
  // double _fun(double x) { return pow(x,0.65)*pow((1-0.99*x),-3.3); } 
  // double _fun(double x) { double n=(1-0.99*x); return sqrt(x)/(n*n*n); } 
  // double _fun(double x) { return 1; } 
  
  // this is significantly quicker than pow(x,1.5)*pow(1-0.99*x,3) 
  static double _fun(double x)      { double n=(1-0.99*x); return sqrt(x*x*x)/(n*n*n); } 
  static double weightfun(double x) { return _fun(x); }
  

  // using log(log(Q2/mLambda)) or just log(log(Q2)) makes 
  // little difference for the LHC range
  static double ftau(double Q2) { return log(log(Q2/0.0625)); }
  static double fQ2(double tau) { return 0.0625*exp(exp(tau)); }

  
  // grid value accessors
  double gety1(int iy)    const { return m_weight[0]->yaxis()[iy]; }   
  double gety2(int iy)    const { return m_weight[0]->zaxis()[iy]; }   
  //  double gety(int iy)     const { return gety1(iy); } 

  double gettau(int itau) const { return m_weight[0]->xaxis()[itau]; } 

  // number of subprocesses

  int SubProcesses() const { return m_Nproc; } 

  // kinematic variable accessors
  // y (x) 
  int    Ny1()      const { return m_weight[0]->yaxis().N(); }    
  double y1min()    const { return m_weight[0]->yaxis().min(); }  
  double y1max()    const { return m_weight[0]->yaxis().max(); }  
  double deltay1()  const { return m_weight[0]->yaxis().delta(); }

  int    Ny2()      const { return m_weight[0]->zaxis().N(); }    
  double y2min()    const { return m_weight[0]->zaxis().min(); }  
  double y2max()    const { return m_weight[0]->zaxis().max(); }  
  double deltay2()  const { return m_weight[0]->zaxis().delta(); }

  //  int    Ny()      const { return Ny1(); } 
  //  double ymin()    const { return y1min(); }  
  //  double ymax()    const { return y2min(); }  
  //  double deltay()  const { return deltay1(); }  

  //  int yorder(int i)   { return m_yorder=i; }  
  int yorder()  const { return m_yorder; }  

  // tau (Q2)
  int    Ntau()     const { return m_weight[0]->xaxis().N(); } 
  double taumin()   const { return m_weight[0]->xaxis().min(); } 
  double taumax()   const { return m_weight[0]->xaxis().max(); } 
  double deltatau() const { return m_weight[0]->xaxis().delta(); } 
  
  //  int tauorder(int i)  { return m_tauorder=i; }  
  int tauorder() const { return m_tauorder; }  


  // maybe these are redundant and should be removed
  int    getNQ2()     const { return m_weight[0]->xaxis().N(); }
  double getQ2min()   const { return fQ2(taumin()); }
  double getQ2max()   const { return fQ2(taumax()); }
  
  int    getNx1()     const { return m_weight[0]->yaxis().N(); }
  double getx1min()   const { return fx(y1max()); }
  double getx1max()   const { return fx(y1min()); }

  int    getNx2()     const { return m_weight[0]->zaxis().N(); }
  double getx2min()   const { return fx(y2max()); }
  double getx2max()   const { return fx(y2min()); }

  //  int    getNx()     const { return getNx1(); } 
  //  double getxmin()   const { return getx2max(); }
  //  double getxmax()   const { return getx1max(); }

  
  static double transformvar()         { return m_transvar; }
  static double transformvar(double v) { return m_transvar=v; }

  bool   symmetrise(bool t=true)   { return m_symmetrise=t; }
  bool   isSymmetric() const       { return m_symmetrise; }

  bool   isOptimised() const       { return m_optimised; }
  bool   setOptimised(bool t=true) { return m_optimised=t; } 



  static bool reweight(bool t=true)    { return m_reweight=t; }


  // setup the pdf grid for calculating the pdf using 
  // interpolation - needed if you actually want 
  // the interpolated values for the pdf's or if you want 
  // the grid to calculate the cross section for you
  void setuppdf(void (*pdf)(const double& , const double&, double* ),
		double (*alphas)(const double&),
		int nloop=0, 
		double rscale_factor=1,
		double fscale_factor=1,
		double beam_scale=1 );

  // get the interpolated pdf's
  void pdfinterp(double x1, double Q2, double* f);

  double convolute(void   (*pdf)(const double& , const double&, double* ), 
		   appl_pdf* genpdf, 
		   double (*alphas)(const double& ), 
		   int     lo_order=0,  
		   int     nloop=0, 
		   double  rscale_factor=1,
		   double  fscale_factor=1,
		   double Escale=1 );


  double convolute_subproc(int subproc, 
			   void   (*pdf)(const double& , const double&, double* ), 
			   appl_pdf* genpdf, 
			   double (*alphas)(const double& ), 
			   int     lo_order=0,  
			   int     nloop=0, 
			   double  rscale_factor=1,
			   double  fscale_factor=1,
			   double Escale=1 );
  

  // some useful algebraic operators
  igrid& operator=(const igrid& g); 
  
  igrid& operator*=(const double& d) { 
    for ( int ip=0 ; ip<m_Nproc ; ip++ ) if ( m_weight[ip] ) (*m_weight[ip]) *= d; 
    return *this;
  } 

  // should really check all the limits and *everything* is the same
  igrid& operator+=(const igrid& g) { 
    for ( int ip=0 ; ip<m_Nproc ; ip++ ) {
      if ( m_weight[ip] && g.m_weight[ip] ) (*m_weight[ip]) += (*g.m_weight[ip]);
    }
    return *this;
  } 

  // ouput header
  ostream& header(ostream& s) const;


private:

  // internal common construct for the different types of constructor
  void construct();

  // cleanup
  void deleteweights();
  void deletepdftable();

  // interpolation section - inline and static internals for calculation of the 
  // interpolation for storing on the grid nodes 

  // x (y) interpolation formula
  int fk1(double x) const {
    double y = fy(x);
    // make sure we are in the range covered by our binning
    if( y<y1min() || y>y1max() ) {
      if ( y<y1min() ) cerr <<"\tWarning: x1 out of range: x=" << x << "\t(y=" << y << ")\tBelow Delx=" << x-fx(y1min());
      else             cerr <<"\tWarning: x1 out of range: x=" << x << "\t(y=" << y << ")\tAbove Delx=" << x-fx(y1min());
      cerr << " ( " <<  fx(y1max())  << " - " <<  fx(y1min())  << " )"  
	   << "\ty=" << y << "\tDely=" << y-y1min() << " ( " << y1min() << " - " <<  y1max()  << " )" << endl;
      //     cerr << "\t" << m_weight[0]->yaxis() << "\n\t" 
      //          << m_weight[0]->yaxis().transform(fx) << endl; 
    }
    int k = (int)((y-y1min())/deltay1() - (m_yorder>>1)); // fast integer divide by 2
    if ( k<0 ) k=0;  
    // shift interpolation end nodes to enforce range
    if(k+m_yorder>=Ny1())  k=Ny1()-1-m_yorder;    
    return k;
  }

  int fk2(double x) const {
    double y = fy(x);
    // make sure we are in the range covered by our binning
    if( y<y2min() || y>y2max() ) {
      if ( y<y2min() ) cerr <<"\tWarning: x2 out of range: x=" << x << "\t(y=" << y << ")\tBelow Delx=" << x-fx(y2min());
      else             cerr <<"\tWarning: x2 out of range: x=" << x << "\t(y=" << y << ")\tAbove Delx=" << x-fx(y2min());
      cerr << " ( " <<  fx(y2max())  << " - " <<  fx(y2min())  << " )"  
	   << "\ty=" << y << "\tdely=" << y-y2min() << " ( " << y2min() << " - " <<  y2max()  << " )" << endl;
      //     cerr << "\t" << m_weight[0]->yaxis() << "\n\t" << m_weight[0]->yaxis().transform(fx) << endl; 
    }
    int k = (int)((y-y2min())/deltay2() - (m_yorder>>1)); // fast integer divide by 2
    if ( k<0 ) k=0;  
    // shift interpolation end nodes to enforce range
    if(k+m_yorder>=Ny2())  k=Ny2()-1-m_yorder;    
    return k;
  }

  
  // Q2 (tau) interpolation formula 
  int fkappa(double Q2) const {
    double tau = ftau(Q2);
    // make sure we are in the range covered by our binning
    if( tau<taumin() || tau>taumax() ) {
      cerr << "\tWarning: Q2 out of range Q2=" << Q2 
	   << "\t ( " << fQ2(taumin()) << " - " << fQ2(taumax()) << " )" << endl;
      //      cerr << "\t" << m_weight[0]->xaxis() << "\n\t" << m_weight[0]->xaxis().transform(fQ2) << endl; 
    }
    int kappa = (int)((tau-taumin())/deltatau() - (m_tauorder>>1)); // fast integer divide by 2
    // shift interpolation end nodes to enforce range
    if(kappa+m_tauorder>=Ntau()) kappa=Ntau()-1-m_tauorder;  
    if(kappa<0) kappa=0;
    return kappa;
  }
  
  //  int _fk(double x)      const { return fk(x);  }
  //  int _fkappa(double Q2) const { return fkappa(Q2); }

  // fast -1^i function
  static int pow1(int i) { return ( 1&i ? -1 : 1 ); } 

  // fast factorial
  static double fac(int i) {
    int j;
    static int ntop = 4;
    static double f[34] = { 1, 1, 2, 6, 24 }; 
    if ( i<0 )  {  cerr << "igrid::fac() negative input"  << endl; return 0;  }
    if ( i>33 ) {  cerr << "igrid::fac() input too large" << endl; return 0;  }
    while ( ntop<i ) {
      j=ntop++;
      f[ntop]=f[j]*ntop;
    } 
    return f[i];
  }

  // although not the best way to calculate interpolation coefficients,
  // it may be the best for our use, where the "y" values of the nodes 
  // are not yet defined at the time of evaluation.
  static double fI(int i, int n, double u) {
    if ( n==0 && i==0 )     return 1.0;
    if ( fabs(u-i)<1e-8 ) return 1.0;
    //    if ( fabs(u-n)<u*1e-8 ) return 1.0;
    double product = pow1(n-i) / ( fac(i)*fac(n-i)*(u-i) );
    for( int z=0 ; z<=n ; z++ )  product *= (u-z);
    return product;
  }
  

private:

  // ranges of interest

  // x <-> y parameters
  int    m_Ny1; 
  double m_y1min; 
  double m_y1max;
  double m_deltay1;

  int    m_Ny2; 
  double m_y2min; 
  double m_y2max;
  double m_deltay2;

  int    m_yorder; 

  // tau <-> Q2 parameters
  int    m_Ntau; 
  double m_taumin; 
  double m_taumax;
  double m_deltatau;
  int    m_tauorder; 

  int    m_Nproc;   // number of subprocesses

  // grid state information

  // useful transform information for storage in root file 
  string                   m_transform;
  static map<const string, transform_vec> m_fmap;

  static double m_transvar;    // transform function parameter
  static bool   m_reweight;    // reweight the pdf?
  
  bool   m_symmetrise;   // symmetrise the grid or not 
  bool   m_optimised;    // optimised?

  // the actual weight grids
  SparseMatrix3d**   m_weight;

  // pdf value table for convolution 
  // (NB: doesn't need to be a class variable)
  double*** m_fg1; 
  double*** m_fg2; 

  // values from the splitting functions
  double*** m_fsplit1; 
  double*** m_fsplit2; 

  // alpha_s table
  double*   m_alphas;

};

};

std::ostream& operator<<(std::ostream& s, const appl::igrid& mygrid);


#endif // __APPL_IGRID_H 
