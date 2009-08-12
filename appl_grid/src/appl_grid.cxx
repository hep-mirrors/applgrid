
//  appl_grid.cxx        

//   grid class - all the functions needed to create and 
//   fill the grid from an NLO calculation program. 
//  
//  Copyright (C) 2007 Mark Sutton (sutt@hep.ucl.ac.uk)    

// $Id: appl_grid.cxx, v1.00 2007/10/16 17:01:39 sutt

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <iomanip>
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

#include <cmath>
using std::abs;
using std::fabs;

#include "appl_grid/appl_pdf.h"
#include "appl_grid/appl_timer.h"
#include "appl_grid/Directory.h"
#include "appl_grid/appl_grid.h"
#include "appl_grid/TFileString.h"
using appl::grid;

#include "TFile.h"
#include "TObjString.h"
#include "TVectorT.h"




const string grid::m_version = "version-1.0";


grid::grid(int NQ2, double Q2min, double Q2max, int Q2order, 
	   int Nx,  double xmin,  double xmax,  int xorder,
	   int Nobs,  double obsmin, double obsmax, 
	   string genpdfname,
	   int leading_order, int nloops, 
	   string transform ) :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_symmetrise(false), 
  m_transform(transform), m_genpdfname(genpdfname) {
  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("reference","Bin-Info for Observable", Nobs, obsmin, obsmax);
  //  m_genpdf = genpdf_map.find(m_genpdfname)->second;
  m_genpdf = appl_pdf::getpdf(m_genpdfname);
  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform); 
}




grid::grid(int Nobs, const double* obsbins, 
	   int NQ2,  double Q2min, double Q2max, int Q2order, 
	   int Nx,   double xmin, double xmax,   int xorder, 
	   string genpdfname, 
	   int leading_order, int nloops, 
	   string transform ) :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_symmetrise(false),
  m_transform(transform), m_genpdfname(genpdfname) {
  
  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("reference","Bin-Info for Observable", Nobs, obsbins);
  //   m_genpdf = genpdf_map.find(m_genpdfname)->second;  
  m_genpdf = appl_pdf::getpdf(m_genpdfname);
  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform );
}




grid::grid(const vector<double> obs, 
	   int NQ2, double Q2min, double Q2max, int Q2order,
	   int Nx,  double xmin,  double xmax,  int xorder, 
	   string genpdfname,
	   int leading_order, int nloops, 
	   string transform )  :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_symmetrise(false),  
  m_transform(transform), m_genpdfname(genpdfname) { 

  if ( obs.size()==0 ) { 
    cerr << "grid::not enough bins in observable" << endl;
    exit(0);
  } 
  
  double* obsbins = new double[obs.size()];  
  for ( int i=0 ; i<obs.size() ; i++ ) obsbins[i] = obs[i];
  int Nobs = obs.size()-1;

  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("reference","Bin-Info for Observable", Nobs, obsbins);
  //  m_genpdf = genpdf_map.find(m_genpdfname)->second;
  m_genpdf = appl_pdf::getpdf(m_genpdfname);
  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform); 
}


grid::grid(const string& filename, const string& dirname)  :
  m_leading_order(0), m_order(0),
  m_optimised(false), m_trimmed(false), m_transform(""), m_symmetrise(false)  {
  
  cout << "grid() reading grid from file " << filename << endl;
  
  TFile gridfile(filename.c_str());
  
  //  gDirectory->cd(dirname.c_str());

  //  cout << "pwd=" << gDirectory->GetName() << endl;

  //  gDirectory->cd(dirname.c_str());

  //  cout << "pwd=" << gDirectory->GetName() << endl;

  //  Directory d(dirname);
  //  d.push();


  // using the title of a TH1D because I don't know 
  // how else to save a string in a root file 
  // TH1D* _transform = (TH1D*)gridfile.Get("Transform");
  // m_transform = string(_transform->GetTitle());
  // delete _transform;

  //  TH1D* _genpdfname = (TH1D*)gridfile.Get("GenPDF");
  //  m_genpdfname = string(_genpdfname->GetTitle());
  //  delete _genpdfname;

  // get the name of the transform pair and the
  // generalised pdf

  TFileString _tags = *(TFileString*)gridfile.Get((dirname+"/Tags").c_str());  
  // TFileString _tags = *(TFileString*)gridfile.Get("Tags");  
  m_transform  = _tags[0];
  m_genpdfname = _tags[1];

  string _version = _tags[2];

  // check it has the correct version
  if ( _version != m_version ) { 
    throw exception(cerr << "incorrect version " << _version << " expected " << m_version ); 
  }

  //  cout << "Tags=" << _tags << endl;

  m_genpdf = appl_pdf::getpdf(m_genpdfname);

#if 0
  if ( genpdf_map.find(m_genpdfname)==genpdf_map.end() )  { 
    cerr << "grid::grid() generalised pdf " << m_genpdfname << " not in map" << endl;
    exit(0);
  }
  else { 
    m_genpdf = genpdf_map.find(m_genpdfname)->second; 
  }
#endif

  //  cout << "grid::grid() read transform " << m_transform << " from file" << endl;

  // read state information
  // hmmm, have to use TVectorT<double> since TVector<int> 
  // apparently has no constructor (???)
  TVectorT<double>* setup=(TVectorT<double>*)gridfile.Get((dirname+"/State").c_str());
  // TVectorT<double>* setup=(TVectorT<double>*)gridfile.Get("State");
  m_run        = int((*setup)(0)+0.5);
  m_optimised  = ( (*setup)(1)!=0 ? true : false );
  m_symmetrise = ( (*setup)(2)!=0 ? true : false );  

  m_leading_order = int((*setup)(3)+0.5);  
  m_order         = int((*setup)(4)+0.5);  

  //  cout << "grid::grid()::m_symmetrise=" << m_symmetrise << endl; 

  delete setup;

  //  cout << "grid::grid() read setup" << endl;

  // Read observable bins information
  //  gridfile.GetObject("obs_bins", m_obs_bins );
  m_obs_bins = (TH1D*)gridfile.Get((dirname+"/reference").c_str());
  m_obs_bins->SetDirectory(0);

  //  cout << "grid::grid() read obs bins" << endl;

  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    //  cout << "grid::grid() iorder=" << iorder << endl;
    m_grids[iorder] = new igrid*[Nobs()];  
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {
      char name[128];  sprintf(name, (dirname+"/weight[alpha-%d][%03d]").c_str(), iorder, iobs);
      // cout << "grid::grid() reading " << name << "\tiobs=" << iobs << endl;

      m_grids[iorder][iobs] = new igrid(gridfile, name);

      //    _size += m_grids[iorder][iobs]->size();
      //      cout << "grid::grid() done" << endl;
    }
  }

  //  d.pop();

  //  cout << "grid::grid() read from file" << endl;
}


grid::grid(const grid& g) : 
  m_obs_bins(new TH1D(*g.m_obs_bins)), 
  m_leading_order(g.m_leading_order), m_order(g.m_order), 
  m_run(g.m_run), m_optimised(g.m_optimised), m_trimmed(g.m_trimmed), m_symmetrise(g.m_symmetrise) {
 
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    m_grids[iorder] = new igrid*[Nobs()];
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  m_grids[iorder][iobs] = new igrid(*g.m_grids[iorder][iobs]);
  }
} 


// Initialize histogram that saves the correspondence obsvalue<->obsbin

// constructor common internals 
void grid::construct(int Nobs, 
		     int NQ2,  double Q2min, double Q2max, int Q2order,
		     int Nx,   double xmin,  double xmax,  int xorder, 
		     int order, 
		     string transform ) { 

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    m_grids[iorder] = new igrid*[Nobs];
    for ( int iobs=0 ; iobs<Nobs ; iobs++ ) {
      m_grids[iorder][iobs] = new igrid(NQ2, Q2min, Q2max, Q2order, 
					Nx, xmin, xmax, xorder, 
					transform, m_genpdf->Nproc());
    }
  }

}



grid::~grid() {
  cout << "~grid() deleting grid" << endl;
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {  
    if( m_grids[iorder] ) { 
      for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  delete m_grids[iorder][iobs];
      delete[] m_grids[iorder];
      m_grids[iorder] = NULL;
    }
  }
  if(m_obs_bins) delete m_obs_bins;
  m_obs_bins=NULL;
}




// algebraic operators

grid& grid::operator=(const grid& g) { 
  // clear out the old...
  delete m_obs_bins;
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  delete m_grids[iorder][iobs];
    delete m_grids[iorder];
  }
  
  // copy the new
  m_obs_bins = new TH1D(*g.m_obs_bins);
  // copy the state
  m_leading_order = g.m_leading_order;
  m_order         = g.m_order;

  m_run       = g.m_run;
  m_optimised = g.m_optimised;
  m_trimmed   = g.m_trimmed;

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    m_grids[iorder] = new igrid*[Nobs()];
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  m_grids[iorder][iobs] = new igrid(*g.m_grids[iorder][iobs]);
  }
} 
  

grid& grid::operator*=(const double& d) { 
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) (*m_grids[iorder][iobs])*=d; 
  }
  return *this;
}


grid& grid::operator+=(const grid& g) {
  m_run      += g.m_run;
  m_optimised = g.m_optimised;
  m_trimmed   = g.m_trimmed;
  if ( Nobs()!=g.Nobs() )   throw exception("grid::operator+ Nobs bin mismatch");
  if ( m_order!=g.m_order ) throw exception("grid::operator+ different order grids");
  if ( m_leading_order!=g.m_leading_order ) throw exception("grid::operator+ different order processes in grids");
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) (*m_grids[iorder][iobs]) += (*g.m_grids[iorder][iobs]); 
  }
  return *this;
}


// fill the appropriate igrid with these weights
void grid::fill(const double x1, const double x2, const double Q2, 
		const double obs, 
		const double* weight, const int iorder)  {  
  int iobs = m_obs_bins->FindBin(obs)-1;
  if ( iobs<0 || iobs>=Nobs() ) {
    //    cerr << "grid::fill() iobs out of range " << iobs << "\tobs=" << obs << endl;
    //    cerr << "obs=" << obs << "\tobsmin=" << obsmin() << "\tobsmax=" << obsmax() << endl;
    return;
  }
  
  //  cout << "iobs=" << iobs << "\tobs=" << obs;
  //  for ( int i=0 ; i<subProcesses() ; i++ ) cout << "\t" << weight[i];
  //  cout << endl;

  //  cout << "\tiobs=" << iobs << endl;
  if ( m_symmetrise && x2<x1 )  m_grids[iorder][iobs]->fill(x2, x1, Q2, weight);
  else                          m_grids[iorder][iobs]->fill(x1, x2, Q2, weight);
}


// fast fill pre-optimisation don't perform the interpolation and so on
void grid::fill_phasespace(const double x1, const double x2, const double Q2, 
			   const double obs, 
			   const double* weight, const int iorder) {
  int iobs = m_obs_bins->FindBin(obs)-1;
  if ( iobs<0 || iobs>=Nobs() ) {
    //  cerr << "grid::fill() iobs out of range " << iobs << "\tobs=" << obs << endl;
    //  cerr << "obs=" << obs << "\tobsmin=" << obsmin() << "\tobsmax=" << obsmax() << endl;
    return;
  }
  if ( m_symmetrise && x2<x1 )  m_grids[iorder][iobs]->fill_phasespace(x2, x1, Q2, weight);
  else                          m_grids[iorder][iobs]->fill_phasespace(x1, x2, Q2, weight);
}




// fast fill pre-optimisation don't perform the interpolation and so on
void grid::fill_index(const int ix1, const int ix2, const int iQ2, 
		      const int iobs, 
		      const double* weight, const int iorder) {

  if ( iobs<0 || iobs>=Nobs() ) {
    //  cerr << "grid::fill() iobs out of range " << iobs << "\tobs=" << obs << endl;
    //  cerr << "obs=" << obs << "\tobsmin=" << obsmin() << "\tobsmax=" << obsmax() << endl;
    return;
  }
  if ( m_symmetrise && ix2<ix1 )  m_grids[iorder][iobs]->fill_index(ix2, ix1, iQ2, weight);
  else                            m_grids[iorder][iobs]->fill_index(ix1, ix2, iQ2, weight);
}



void grid::trim() {
  m_trimmed = true;
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) m_grids[iorder][iobs]->trim(); 
  }
}

void grid::untrim() {
  m_trimmed = false;
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) m_grids[iorder][iobs]->untrim(); 
  }
}

void grid::print() const {
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {     
      cout << iobs << "\t" 
	   << setprecision(5) << setw(6) << getReference()->GetBinLowEdge(iobs+1) << "\t- " 
	   << setprecision(5) << setw(6) << getReference()->GetBinLowEdge(iobs+2) << "\t"; 
      m_grids[iorder][iobs]->print();       
    }
  }
}

void grid::setuppdf(void (*pdf)(const double&, const double&, double* ) )  {  }
void grid::pdfinterp(double x, double Q2, double* f) {  }


// dump to file
void grid::Write(const string& filename, const string& dirname) {
  
  string _filename(filename);

  if ( FILE* f=fopen(_filename.c_str(), "r") ) { 
    fclose(f);
    _filename += "-save";
    string cmd = "mv " + filename + " " + _filename;
    int i = system(cmd.c_str());
  } 

  //  cout << "grid::Write() writing to file " << _filename << endl;
  //  TFile rootfile(_filename.c_str(),"recreate");

  cout << "grid::Write() writing to file " << filename << endl;
  TFile rootfile(filename.c_str(),"recreate");


  cout << "pwd=" << gDirectory->GetName() << endl;

  Directory d(dirname);
  d.push();
  
  cout << "pwd=" << gDirectory->GetName() << endl;

  // write the name of the transform pair and the
  // generalised pdf
  TFileString _tags("Tags");
  _tags.add(m_transform);
  _tags.add(m_genpdfname);
  _tags.add(m_version);
  _tags.Write();

  //  TH1D* _transform = new TH1D("Transform", m_transform.c_str(), 1, 0, 1);  
  //  _transform->Write();

  //  TH1D* _genpdfname = new TH1D("GenPDF", m_genpdfname.c_str(), 1, 0, 1);
  //  _genpdfname->Write();


  // state information
  TVectorT<double>* setup=new TVectorT<double>(5); 
  (*setup)(0) = m_run;
  (*setup)(1) = ( m_optimised  ? 1 : 0 );
  (*setup)(2) = ( m_symmetrise ? 1 : 0 );
  (*setup)(3) =   m_leading_order ;
  (*setup)(4) =   m_order ;
  setup->Write("State");
  
  //  int _size     = 0;
  //  int trim_size = 0;

  // internal grids
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {
      char name[128];  sprintf(name, "weight[alpha-%d][%03d]", iorder, iobs);
      // cout << "writing grid " << name << endl;
      //   _size += m_grids[iorder][iobs]->size();
      m_grids[iorder][iobs]->write(name);
      //   trim_size += m_grids[iorder][iobs]->size();
    }
  }
  //  cout <<"grid::Write() size(untrimmed)=" << _size 
  //     << "\tsize(trimmed)="              << trim_size << endl;
  //  d.pop();
  m_obs_bins->Write();
  rootfile.Close();
  d.pop();
}


// takes pdf as the pdf lib wrapper for the pdf set for the convolution.
// type specifies which sort of partons should be included:

TH1D* grid::convolute(void (*pdf)(const double& , const double&, double* ), 
		      double (*alphas)(const double& ), 
		      int     nloops, 
		      double  rscale_factor,
		      double  fscale_factor,
		      void (*splitting)(const double& , const double&, double* ) )
{ 
  
  struct timeval _ctimer = appl_timer_start();

  TH1D* h = new TH1D(*m_obs_bins);
  h->SetName("xsec");

  string label;

  if ( nloops>=m_order ) { 
    cerr << "too many loops for grid nloops=" << nloops << "\tgrid=" << m_order << endl;   
    return h;
  } 

  for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  

    double dsigma = 0;
   
    if ( nloops==0 ) {
      label = "lo      ";
      // leading order cross section
      dsigma = m_grids[0][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order, 0);
    }
    else if ( nloops==1 ) { 
      label = "nlo     ";
      // next to leading order cross section
      // cout << "convolute() nloop=1" << endl;
      // leading order contribution and scale dependent born dependent terms
      double dsigma_lo  = m_grids[0][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order, 1, rscale_factor, fscale_factor, splitting);
      // cout << "dsigma_lo=" << dsigma_lo << endl;
      // next to leading order contribution
      //      double dsigma_nlo = m_grids[1][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order+1, 0);
      // GPS: the NLO piece must use the same rscale_factor and fscale_factor as
      //      the LO piece -- that's the convention that defines how NLO calculations
      //      are done.
      double dsigma_nlo = m_grids[1][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order+1, 0, rscale_factor, fscale_factor);
      // cout << "dsigma_nlo=" << dsigma_nlo << endl;
      dsigma = dsigma_lo + dsigma_nlo;
    }
    else if ( nloops==-1 ) {
      label = "nlo only";
      // nlo contribution only (only strict nlo contributions) 
      dsigma = m_grids[1][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order+1, 0, rscale_factor, fscale_factor);
    } 
    else if ( nloops==2 ) {
      // FIXME: not implemented completely yet 
      return h;
      label = "nnlo    ";
      // next to next to leading order contribution 
      // NB: NO scale dependendent parts so only  muR=muF=mu
      double dsigma_lo  = m_grids[0][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order, 1);
      // next to leading order contribution      
      double dsigma_nlo = m_grids[1][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order+1, 0);
      // next to next to leading order contribution
      double dsigma_nnlo = m_grids[2][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order+2, 0);
      dsigma = dsigma_lo + dsigma_nlo + dsigma_nnlo;
    }


    double deltaobs = h->GetBinLowEdge(iobs+2)-h->GetBinLowEdge(iobs+1);
    h->SetBinContent(iobs+1, dsigma/(deltaobs));
    h->SetBinError(iobs+1, 0);

    //    cout << "dsigma[" << iobs << "]=" << dsigma/deltaobs << endl;
  }  // iobs   

  double _ctime = appl_timer_stop(_ctimer);

  cout << "grid::convolute() " << label << " convolution time=" << _ctime << " ms" << endl;

  return h;
}







TH1D* grid::convolute_subproc(int subproc,
			      void (*pdf)(const double& , const double&, double* ), 
			      double (*alphas)(const double& ), 
			      int     nloops, 
			      double  rscale_factor,
			      double  fscale_factor,
			      void (*splitting)(const double& , const double&, double* ) )
{ 
  
  struct timeval _ctimer = appl_timer_start();

  //  genpdf = genpdf_map[m_genpdfname];
    
  TH1D* h = new TH1D(*m_obs_bins);
  h->SetName("xsec");

  string label;

  int lo_order = m_leading_order;

  if ( nloops>=m_order ) { 
    cerr << "too many loops for grid nloops=" << nloops << "\tgrid=" << m_order << endl;   
    return h;
  } 
  
  for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  

    double dsigma = 0;
   
    if ( nloops==0 ) {
      label = "lo      ";
      //      cout << "convolute() nloop=0" << iobs << endl;
      // leading order cross section
      dsigma = m_grids[0][iobs]->convolute_subproc(subproc, pdf, m_genpdf, alphas, lo_order, 0);
    }
    else if ( nloops==1 ) { 
      label = "nlo     ";
      // cout << "convolute() nloop=1" << endl;
      // next to leading order cross section
      // leading order contribution and scale dependent born dependent terms
      double dsigma_lo  = m_grids[0][iobs]->convolute_subproc(subproc, pdf, m_genpdf, alphas, lo_order,   1, rscale_factor, fscale_factor, splitting);
      // cout << "dsigma_lo=" << dsigma_lo << endl;
      // next to leading order contribution
      double dsigma_nlo = m_grids[1][iobs]->convolute_subproc(subproc, pdf, m_genpdf, alphas, lo_order+1, 0, rscale_factor, fscale_factor, splitting );
      // cout << "dsigma_nlo=" << dsigma_nlo << endl;
      dsigma = dsigma_lo + dsigma_nlo;
    }
    else if ( nloops==-1 ) { 
      label = "nlo only";
      // nlo contribution only (only strict nlo contributions) 
      dsigma = m_grids[1][iobs]->convolute_subproc(subproc, pdf, m_genpdf, alphas, lo_order+1, 0);
    }
    else if ( nloops==2 ) { 
      // FIXME: not implemented completely yet 
      return h;
      label = "nnlo    ";
      // next to next to leading order contribution 
      // NB: NO scale dependendent parts, so only muR=muF=mu
      double dsigma_lo  = m_grids[0][iobs]->convolute_subproc(subproc, pdf, m_genpdf, alphas, lo_order, 0);
      double dsigma_nlo = m_grids[1][iobs]->convolute_subproc(subproc, pdf, m_genpdf, alphas, lo_order+1, 0);
      double dsigma_nnlo = m_grids[2][iobs]->convolute_subproc(subproc, pdf, m_genpdf, alphas, lo_order+2, 0);
      dsigma = dsigma_lo + dsigma_nlo + dsigma_nnlo;
    }


    double deltaobs = h->GetBinLowEdge(iobs+2)-h->GetBinLowEdge(iobs+1);
    h->SetBinContent(iobs+1, dsigma/(deltaobs));
    h->SetBinError(iobs+1, 0);

    //    cout << "dsigma[" << iobs << "]=" << dsigma/deltaobs << endl;

    
    //    cout << "obs bin " << iobs 
    //         << "\t" <<  h->GetBinLowEdge(iobs+1) << " - " << h->GetBinLowEdge(iobs+2)
    //	       << "\tdsigma=" << dsigma << endl;
  }  // iobs   

  double _ctime = appl_timer_stop(_ctimer);

  cout << "grid::convolute_subproc(" << subproc << ") " << label << " convolution time=" << _ctime << " ms" << endl;

  return h;
}





void grid::optimise() {
  m_optimised = true;
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  { 
      cout << "grid::optimise() bin " << iobs << "\t";
      m_grids[iorder][iobs]->optimise();
    }
  }
}

void grid::optimise(int NQ2, int Nx) {  optimise(NQ2, Nx, Nx);  }

void grid::optimise(int NQ2, int Nx1, int Nx2) {
  m_optimised = true;
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  { 
      cout << "grid::optimise() bin " << iobs << "\t";
      m_grids[iorder][iobs]->optimise(NQ2, Nx1, Nx2);
    }
  }
}




// redefine the limits by hand
void grid::redefine(int iobs, int iorder,
		    int NQ2, double Q2min, double Q2max, 
		    int Nx,  double  xmin, double  xmax ) 
{ 
  
  if  ( iorder>=m_order ) { 
    cerr << "grid does not extend to this order" << endl;
    return;
  }
  
  if ( iobs<0 || iobs>=Nobs() ) { 
    cerr << "observable bin out of range" << endl;
    return;
  }
  
  if ( iorder==0 ) { 
    cout << "grid::redefine() iobs=" << iobs 
	 << "NQ2="  << NQ2 << "\tQmin=" << sqrt(Q2min) << "\tQmax=" << sqrt(Q2max) 
	 << "\tNx=" << Nx  << "\txmin=" <<       xmin  << "\txmax=" <<       xmax << endl; 
  }
  
  igrid* oldgrid =  m_grids[iorder][iobs];
  
  //  m_grids[iorder][iobs]->redefine(NQ2, Q2min, Q2max, Nx, xmin, xmax);

  m_grids[iorder][iobs] = new igrid(NQ2, Q2min, Q2max, oldgrid->tauorder(),
				    Nx,  xmin,  xmax,  oldgrid->yorder(), 
				    oldgrid->transform(), m_genpdf->Nproc());

  delete oldgrid;
}
  





ostream& operator<<(ostream& s, const appl::grid& g) {
  s << "==================================================" << endl;
  s << "appl::grid version " << g.version() << "\t(" << g.subProcesses() << " initial states, " << g.Nobs() << " observable bins)" << endl;
  if ( g.isOptimised() ) s << "Optimised grid" << endl;
  if ( g.isSymmetric() ) s << "Symmetrised in x1, x2" << endl;
  else                   s << "Unsymmetrised in x1, x2" << endl;
  s << "leading order of processes  "  << g.leadingOrder() << endl;
  s << "number of loops for grid    " << g.nloops() << endl;   
  s << "x->y coordinate transform:  "  << g.getTransform() << endl;
  s << "genpdf in use: " << g.getGenpdf() << endl;
  s << "--------------------------------------------------" << endl;
  s << "Observable binning: [ " << g.Nobs() 
    << " bins : " << g.obsmin() << ",  " << g.obsmax() << " ]" << endl;

  //  for( int iorder=0 ; iorder<1 ; iorder++ ) {
  for( int iobs=0 ; iobs<g.Nobs() ; iobs++ ) {
    s << iobs << "\t" 
      << setprecision(5) << setw(5) << g.getReference()->GetBinLowEdge(iobs+1) << "\t- " 
      << setprecision(5) << setw(5) << g.getReference()->GetBinLowEdge(iobs+2) << "\t"; 
    s << "   " << *(g.weightgrid(0,iobs)) << endl;
  }
  //  }
  s << endl;
  
  return s;
}

