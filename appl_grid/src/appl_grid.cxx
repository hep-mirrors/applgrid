
//  appl_grid.cxx        

//   grid class - all the functions needed to create and 
//   fill the grid from an NLO calculation program. 
//  
//  Copyright (C) 2007 Mark Sutton (sutt@hep.ucl.ac.uk)    

// $Id: appl_grid.cxx, v1.00 2007/10/16 17:01:39 sutt $

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

#include <iostream>
#include <iomanip>
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

#include <cmath>
// using std::abs;
// using std::fabs;

#include "appl_grid/appl_pdf.h"
#include "appl_grid/appl_timer.h"
#include "appl_grid/Directory.h"
#include "appl_grid/appl_grid.h"
#include "appl_grid/TFileString.h"
#include "appl_grid/TFileVector.h"
using appl::grid;



#include "TFile.h"
#include "TObjString.h"
#include "TVectorT.h"



/// this is a compatability flag for persistent versions 
/// of the grid
/// NB: ONLY change this if the persistent class
///     changes in a non-backwards compatible way.

const string grid::m_version = "version-3.0";

/// check if we have hoppet included 
#include "amconfig.h"

#include "appl_grid/hoppet_init.h"

#ifdef HAVE_HOPPET

#include "hoppet_v1.h"

// include hoppet splitting function code

static hoppet_init* hoppet = NULL;

void Splitting(const double& x, const double& Q, double* xf) {
  static const int nLoops    = 1;
  static const int nFlavours = 5;
  hoppetevalsplit_( x, Q, nLoops, nFlavours, xf); 
  return;
}

#else

void Splitting(const double& x, const double& Q, double* xf) {
  throw grid::exception( std::cerr << "hoppet library not included - cannot call splitting function"  ); 
  return; // technically, don't need this - should throw an exception
}

#endif



grid::grid(int NQ2, double Q2min, double Q2max, int Q2order, 
	   int Nx,  double xmin,  double xmax,  int xorder,
	   int Nobs,  double obsmin, double obsmax, 
	   string genpdfname,
	   int leading_order, int nloops, 
	   string transform ) :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_normalised(false), m_symmetrise(false), 
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),
  m_documentation("") {
  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", Nobs, obsmin, obsmax);
  m_obs_bins->SetDirectory(0);

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform); 
}




grid::grid(int Nobs, const double* obsbins, 
	   int NQ2,  double Q2min, double Q2max, int Q2order, 
	   int Nx,   double xmin, double xmax,   int xorder, 
	   string genpdfname, 
	   int leading_order, int nloops, 
	   string transform ) :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false),  m_normalised(false), m_symmetrise(false),
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),
  m_documentation("") {
  
  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", Nobs, obsbins);
  m_obs_bins->SetDirectory(0);

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform );
}




grid::grid(const vector<double> obs, 
	   int NQ2, double Q2min, double Q2max, int Q2order,
	   int Nx,  double xmin,  double xmax,  int xorder, 
	   string genpdfname,
	   int leading_order, int nloops, 
	   string transform )  :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_normalised(false), m_symmetrise(false),  
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),
  m_documentation("") {
  
  if ( obs.size()==0 ) { 
    cerr << "grid::not enough bins in observable" << endl;
    exit(0);
  } 
  
  double* obsbins = new double[obs.size()];  
  for ( unsigned i=0 ; i<obs.size() ; i++ ) obsbins[i] = obs[i];
  int Nobs = obs.size()-1;

  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", Nobs, obsbins);
  m_obs_bins->SetDirectory(0);
  delete[] obsbins;

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform); 
}



grid::grid(const vector<double> obs, 
	   string genpdfname,
	   int leading_order, int nloops, 
	   string transform )  :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_normalised(false), m_symmetrise(false),  
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),  
  m_documentation("")
{ 

  if ( obs.size()==0 ) { 
    cerr << "grid::not enough bins in observable" << endl;
    exit(0);
  } 
  
  double* obsbins = new double[obs.size()];  
  for ( unsigned i=0 ; i<obs.size() ; i++ ) obsbins[i] = obs[i];
  int Nobs = obs.size()-1;

  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", Nobs, obsbins);
  m_obs_bins->SetDirectory(0);
  delete[] obsbins;

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) m_grids[iorder] = new igrid*[Nobs];

}



grid::grid(const string& filename, const string& dirname)  :
  m_leading_order(0),  m_order(0),
  m_optimised(false),  m_trimmed(false), 
  m_normalised(false),
  m_symmetrise(false), m_transform(""), 
  m_applyCorrections(false),
  m_documentation("") 
{

  struct stat stfileinfo;
  if ( stat(filename.c_str(),&stfileinfo) )   {    
    throw exception(std::cerr << "grid::grid() cannot open file " << filename << std::endl ); 
  }

  std::cout << "grid() reading grid from file " << filename << std::endl;
  
  TFile* gridfilep = TFile::Open(filename.c_str());
  
  if (gridfilep->IsZombie()) {
    throw exception(std::cerr << "grid::grid() cannot open file: zombie " << filename << std::endl ); 
    delete gridfilep;
  }

  // TFile gridfile(filename.c_str());
  
  //  gDirectory->cd(dirname.c_str());

  //  cout << "pwd=" << gDirectory->GetName() << endl;

  //  gDirectory->cd(dirname.c_str());

  //  cout << "pwd=" << gDirectory->GetName() << endl;

  //  Directory d(dirname);
  //  d.push();

  TFileString _tags = *(TFileString*)gridfilep->Get((dirname+"/Tags").c_str());  
  // TFileString _tags = *(TFileString*)gridfile.Get("Tags");  
  m_transform  = _tags[0];
  m_genpdfname = _tags[1];

  string _version = _tags[2];

  if ( _tags.size()>3 ) m_documentation = _tags[3];

  // check it has the correct version
  // if ( _version > m_version ) { 
  //      throw exception(cerr << "incorrect version " << _version << " expected " << m_version ); 
  // }
  //  m_version = _version;
  
  std::cout << "appl::grid " << m_version << "\t" << m_documentation << std::endl; 
  
  //  cout << "Tags=" << _tags << endl;

  //  cout << "grid::grid() read transform " << m_transform << " from file" << endl;

  // read state information
  // hmmm, have to use TVectorT<double> since TVector<int> 
  // apparently has no constructor (???)
  TVectorT<double>* setup=(TVectorT<double>*)gridfilep->Get((dirname+"/State").c_str());
 
  m_run        = (*setup)(0);
  m_optimised  = ( (*setup)(1)!=0 ? true : false );
  m_symmetrise = ( (*setup)(2)!=0 ? true : false );  

  m_leading_order = int((*setup)(3)+0.5);  
  m_order         = int((*setup)(4)+0.5);  

  if ( setup->GetNoElements()>5 ) m_cmsScale = (*setup)(5);
  else                            m_cmsScale = 0;
 
  if ( setup->GetNoElements()>6 ) m_normalised = ( (*setup)(6)!=0 ? true : false );
  else                            m_normalised = true;

  if ( setup->GetNoElements()>7 ) m_applyCorrections = ( (*setup)(7)!=0 ? true : false );
  else                            m_applyCorrections = false;

 
  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  delete setup;

  //  cout << "grid::grid() read setup" << endl;

  // Read observable bins information
  //  gridfile.GetObject("obs_bins", m_obs_bins );
  m_obs_bins = (TH1D*)gridfilep->Get((dirname+"/reference").c_str());
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Scale(run());
  m_obs_bins->SetName("referenceInternal");

  //  cout << "grid::grid() read obs bins" << endl;

  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    //  cout << "grid::grid() iorder=" << iorder << endl;
    m_grids[iorder] = new igrid*[Nobs()];  
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {
      char name[128];  sprintf(name, (dirname+"/weight[alpha-%d][%03d]").c_str(), iorder, iobs);
      // cout << "grid::grid() reading " << name << "\tiobs=" << iobs << endl;

      m_grids[iorder][iobs] = new igrid(*gridfilep, name);

      //    _size += m_grids[iorder][iobs]->size();
      //      cout << "grid::grid() done" << endl;
    }
  }

  //  d.pop();

  /// bin-by-bin correction labels                                       
  TFileString* correctionLabels = (TFileString*)gridfilep->Get((dirname+"/CorrectionLabels").c_str());  
  if ( correctionLabels ) { 
    for ( unsigned i=0 ; i<correctionLabels->size() ; i++ ) {
      m_correctionLabels.push_back( (*correctionLabels)[i] ); // copy the correction label
    }
  }

  /// bin-by-bin correction values
  TFileVector* corrections = (TFileVector*)gridfilep->Get((dirname+"/Corrections").c_str());  
  if ( corrections ) { 
    for ( unsigned i=0 ; i<corrections->size() ; i++ ) {
      m_corrections.push_back( (*corrections)[i] ); // copy the correction histograms
    }
  }

  //  cout << "grid::grid() read from file Nobs = " << Nobs() << endl;

  //  std::cout << "read grid" << std::endl;

  delete gridfilep;
}


grid::grid(const grid& g) : 
  m_obs_bins(new TH1D(*g.m_obs_bins)), 
  m_leading_order(g.m_leading_order), m_order(g.m_order), 
  m_run(g.m_run), m_optimised(g.m_optimised), m_trimmed(g.m_trimmed), 
  m_normalised(true),
  m_symmetrise(g.m_symmetrise),
  m_transform(g.m_transform),
  m_genpdfname(g.m_genpdfname), 
  m_cmsScale(g.m_cmsScale),
  m_applyCorrections(g.m_applyCorrections), 
  m_documentation(g.m_documentation)
{
  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

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

  //  std::cout << "appl::grid::construct() m_order " << m_order << "\tNobs " << Nobs << std::endl; 

  for ( int iorder=0 ; iorder<3 ; iorder++ ) m_grids[iorder] = 0;

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    m_grids[iorder] = new igrid*[Nobs];
    for ( int iobs=0 ; iobs<Nobs ; iobs++ ) {
      m_grids[iorder][iobs] = new igrid(NQ2, Q2min, Q2max, Q2order, 
					Nx, xmin, xmax, xorder, 
					transform, m_genpdf[iorder]->Nproc());
    }
  }
  //  std::cout << "appl::grid::construct() return" << std::endl; 
}




  // number of subprocesses 
int grid::subProcesses(int i) const { 
  if ( i<0 || i>=m_order ) throw exception( std::cerr << "grid::subProcess(int i) " << i << " out or range [0-" << m_order-1 << "]" << std::endl );
  return m_grids[i][0]->SubProcesses();     
}  




// add a single grid
void grid::add_igrid(int bin, int order, igrid* g) { 

  if ( !(order>=0 && order<m_order) ) { 
    std::cerr << "grid::add_igrid() order out of range " << order << std::endl; 
    return;
  } 

  if ( !(bin>=0 && bin<Nobs() ) ) {
    std::cerr << "grid::add_igrid() observable bin out of range " << bin << std::endl; 
    return;
  }

  m_grids[order][bin] = g;

  if ( g->transform()!=m_transform ) { 
    std::cerr << "grid::add_igrid() transform " << m_transform 
	      << " does not match that from added igrid, " << g->transform() << std::endl;
    m_transform = g->transform();
  }

}



grid::~grid() {
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {  
    if( m_grids[iorder] ) { 
      for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) { 
	delete m_grids[iorder][iobs];
      }
      delete[] m_grids[iorder];
      m_grids[iorder] = 0;
    }
  }
  if(m_obs_bins) delete m_obs_bins;
  m_obs_bins=0;

#ifdef HAVE_HOPPET
  if ( hoppet ) delete hoppet; 
  hoppet=0; 
#endif

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
  return *this;
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


/// check grids match
bool grid::operator==(const grid& g) const {
  
  bool match = true; 

  if ( Nobs()!=g.Nobs() )    match = false;
  if ( m_order!=g.m_order )  match = false;
  if ( m_leading_order!=g.m_leading_order ) match = false;
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) match &= ( (*m_grids[iorder][iobs]) == (*g.m_grids[iorder][iobs]) ); 
  }
  return match;
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
  //  for ( int i=0 ; i<subProcesses(iorder) ; i++ ) cout << "\t" << weight[i];
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




 /// get the required pdf combinations from those registered   
void grid::findgenpdf( std::string s ) { 
    std::vector<std::string> names = parse( s, ":" );
    if ( names.size()==unsigned(m_order) ) for ( int i=0 ; i<3 ; i++ ) m_genpdf[i] = appl_pdf::getpdf( names[i] );
    else  if ( names.size()==1 ) m_genpdf[0] = m_genpdf[1] = m_genpdf[2] = appl_pdf::getpdf( names[0] );
    else  { 
      throw exception( std::cerr << "requested " << m_order << " pdf combination but given " << names.size() << std::endl );
    }
}


void grid::addpdf( std::string s ) {

    /// parse names, if they are contain .dat, then create the new generic pdfs
    /// they will be added to the pdf map automatically 
    std::vector<std::string> names = parse( s, ":" );

    unsigned imax = unsigned(m_order); 

    /// check to see whether we have a different pdf for each order
    if ( names.size()!=imax ) { 
      if ( names.size()==1 ) imax = 1;
      else { 
	throw exception( std::cerr << "requested " << m_order << " pdf combination but given " << names.size() << std::endl );
      }
    }

    /// loop through all the required pdfs checking whether they exist already,
    /// if not (from thrown exception) then create it, otherwise, don't need to 
    /// do anything 
    for ( unsigned i=0 ; i<imax ; i++ ) { 
      if ( names[i].find(".dat")!=std::string::npos ) { 
	try {
	  appl_pdf::getpdf(names[i]);
	}
	catch ( appl_pdf::exception e ) { 
	  std::cout << "creating new generic_pdf " << names[i] << std::endl;
	  new generic_pdf(names[i]);
	}
      }
    } 

}



void grid::setuppdf(void (*pdf)(const double&, const double&, double* ) )  {  }
// void grid::pdfinterp(double x, double Q2, double* f) {  }

// set the rewight flag of the internal grids
bool grid::reweight(bool t) { 
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {     
      m_grids[iorder][iobs]->reweight(t);       
    }
  }
  return t;
}

// dump to file
void grid::Write(const string& filename, const string& dirname) { 
 
  string _filename(filename);

  if ( FILE* f=fopen(_filename.c_str(), "r") ) { 
    fclose(f);
    _filename += "-save";
    string cmd = "mv " + filename + " " + _filename;
    //    int i = 
    system(cmd.c_str());
  } 

  //  cout << "grid::Write() writing to file " << _filename << endl;
  //  TFile rootfile(_filename.c_str(),"recreate");

  //  cout << "grid::Write() writing to file " << filename << endl;
  TFile rootfile(filename.c_str(),"recreate");

  //  cout << "pwd=" << gDirectory->GetName() << endl;

  Directory d(dirname);
  d.push();
  
  //  cout << "pwd=" << gDirectory->GetName() << endl;

  // write the name of the transform pair and the
  // generalised pdf
  TFileString _tags("Tags");
  _tags.add(m_transform);
  _tags.add(m_genpdfname);
  _tags.add(m_version);
  if ( m_documentation!="" ) _tags.add(m_documentation);
  _tags.Write();

  //  TH1D* _transform = new TH1D("Transform", m_transform.c_str(), 1, 0, 1);  
  //  _transform->Write();

  //  TH1D* _genpdfname = new TH1D("GenPDF", m_genpdfname.c_str(), 1, 0, 1);
  //  _genpdfname->Write();


  //  cout << "state vector=" << endl;

  // state information
  TVectorT<double>* setup=new TVectorT<double>(10); // add a few extra just in case 
  (*setup)(0) = m_run;
  (*setup)(1) = ( m_optimised  ? 1 : 0 );
  (*setup)(2) = ( m_symmetrise ? 1 : 0 );
  (*setup)(3) =   m_leading_order ;
  (*setup)(4) =   m_order ;
  (*setup)(5) =   m_cmsScale ;
  (*setup)(6) = ( m_normalised ? 1 : 0 );
  (*setup)(7) = ( m_applyCorrections ? 1 : 0 );
  setup->Write("State");
  
  //  int _size     = 0;
  //  int trim_size = 0;

  //  cout << "grids Nobs = " << Nobs() << endl;

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

  //  cout << "reference" << endl;
  
  TH1D* reference = (TH1D*)m_obs_bins->Clone("reference");
  
  if ( !getNormalised() )
    if ( run() ) reference->Scale(1/double(run()));
  
  reference->Write();
  delete reference;

  //  cout << "corrections" << endl;

  /// correction histograms

  if ( m_corrections.size()>0 ) {

    /// Fixme: should add the labels to the actual corrections rather than save separately
    /// write labels
    TFileVector* corrections = new TFileVector("Corrections");
    for ( unsigned i=0 ; i<m_corrections.size() ; i++ )  corrections->add( m_corrections[i] );    
    corrections->Write("Corrections");

    /// write actual corrections
    TFileString correctionLabels("CorrectionLabels");
    for ( unsigned i=0 ; i<m_correctionLabels.size() ; i++ )  correctionLabels.add( m_correctionLabels[i] );
    correctionLabels.Write("CorrectionLabels");

  }

  //  std::cout << "close file" << std::endl;

  rootfile.Close();
  d.pop();

  //  std::cout << "written" << std::endl;
}



// takes pdf as the pdf lib wrapper for the pdf set for the convolution.
// type specifies which sort of partons should be included:

std::vector<double> grid::vconvolute(void (*pdf)(const double& , const double&, double* ), 
				     double (*alphas)(const double& ), 
				     int     nloops, 
				     double  rscale_factor,
				     double  fscale_factor,
				     double Escale )
{ 
  
  //  struct timeval _ctimer = appl_timer_start();
  
  double Escale2 = 1;
 
  if ( Escale!=1 ) Escale2 = Escale*Escale;
  
  std::vector<double> hvec;

  double invNruns = 1;
  if ( (!m_normalised) && run() ) invNruns /= double(run());

  //  std::cout << "grid::run() " << run() << std::endl; 

  //  TH1D* h = new TH1D(*m_obs_bins);
  //  h->SetName("xsec");

#ifdef HAVE_HOPPET
  // check if we need to use the splitting function, and if so see if we 
  // need to initialise it again, and do so if required
  if ( fscale_factor!=1 ) {
    if ( hoppet == NULL ) hoppet = new hoppet_init();
    bool newpdf = hoppet->compareCache(pdf);
    if ( newpdf ) hoppet->fillCache( pdf );
  }
#endif

  string label;

  if ( nloops>=m_order ) { 
    cerr << "too many loops for grid nloops=" << nloops << "\tgrid=" << m_order << endl;   
    return hvec;
  } 

  for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  

    double dsigma = 0;
   
    if ( nloops==0 ) {
      label = "lo      ";
      // leading order cross section
      dsigma = m_grids[0][iobs]->convolute(pdf, m_genpdf[0], alphas, m_leading_order, 0, 1, 1, Escale);
    }
    else if ( nloops==1 ) { 
      label = "nlo     ";
      // next to leading order cross section
      // cout << "convolute() nloop=1" << endl;
      // leading order contribution and scale dependent born dependent terms
      double dsigma_lo  = m_grids[0][iobs]->convolute(pdf, m_genpdf[0], alphas, m_leading_order, 1, rscale_factor, fscale_factor, Escale);
      // cout << "dsigma_lo=" << dsigma_lo << endl;
      // next to leading order contribution
      //      double dsigma_nlo = m_grids[1][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order+1, 0);
      // GPS: the NLO piece must use the same rscale_factor and fscale_factor as
      //      the LO piece -- that's the convention that defines how NLO calculations
      //      are done.
      double dsigma_nlo = m_grids[1][iobs]->convolute(pdf, m_genpdf[1], alphas, m_leading_order+1, 0, rscale_factor, fscale_factor, Escale);
      // cout << "dsigma_nlo=" << dsigma_nlo << endl;
      dsigma = dsigma_lo + dsigma_nlo;
    }
    else if ( nloops==-1 ) {
      label = "nlo only";
      // nlo contribution only (only strict nlo contributions) 
      dsigma = m_grids[1][iobs]->convolute(pdf, m_genpdf[1], alphas, m_leading_order+1, 0, rscale_factor, fscale_factor, Escale);
    } 
    else if ( nloops==2 ) {
      // FIXME: not implemented completely yet 
      //      return hvec;
      label = "nnlo    ";
      // next to next to leading order contribution 
      // NB: NO scale dependendent parts so only  muR=muF=mu
      double dsigma_lo  = m_grids[0][iobs]->convolute(pdf, m_genpdf[0], alphas, m_leading_order, 0);
      // next to leading order contribution      
      double dsigma_nlo = m_grids[1][iobs]->convolute(pdf, m_genpdf[1], alphas, m_leading_order+1, 0);
      // next to next to leading order contribution
      double dsigma_nnlo = m_grids[2][iobs]->convolute(pdf, m_genpdf[2], alphas, m_leading_order+2, 0);
      dsigma = dsigma_lo + dsigma_nlo + dsigma_nnlo;
    }
    else if ( nloops==-2 ) {
      label = "nnlo only";
      // next to next to leading order contribution
      dsigma = m_grids[2][iobs]->convolute(pdf, m_genpdf[2], alphas, m_leading_order+2, 0);
    }

    //   double deltaobs = h->GetBinLowEdge(iobs+2)-h->GetBinLowEdge(iobs+1);
    //   h->SetBinContent(iobs+1, dsigma/(deltaobs));
    //   h->SetBinError(iobs+1, 0);

    double deltaobs = m_obs_bins->GetBinLowEdge(iobs+2)-m_obs_bins->GetBinLowEdge(iobs+1);
    
    hvec.push_back( invNruns*Escale2*dsigma/deltaobs );
    // hvec.push_back( Escale2*dsigma/deltaobs );

    //    cout << "dsigma[" << iobs << "]=" << dsigma/deltaobs << endl;
  }  // iobs   

  //  double _ctime = appl_timer_stop(_ctimer);
  //  cout << "grid::convolute() " << label << " convolution time=" << _ctime << " ms" << endl;

  if ( getApplyCorrections() ) applyCorrections(hvec);

  return hvec;
}



#if 0
double grid::vconvolute(void (*pdf)(const double& , const double&, double* ), 
			double (*alphas)(const double& ), 
			int     nloops, 
			double  rscale_factor,
			double  fscale_factor,
			double Escale )
{ 
  
  
}
#endif


std::vector<double> grid::vconvolute_subproc(int subproc,
					     void (*pdf)(const double& , const double&, double* ), 
					     double (*alphas)(const double& ), 
					     int     nloops, 
					     double  rscale_factor, double Escale )
{ 
  
  //  struct timeval _ctimer = appl_timer_start();

  double Escale2 = 1;
 
  if ( Escale!=1 ) Escale2 = Escale*Escale;
  
  double invNruns = 1;
  if ( (!m_normalised) && run() ) invNruns /= double(run());

#ifdef HAVE_HOPPET
  //  factorisation scale variation is disabled for the subprocess
  //  convolution
  //  // check if we need to use the splitting function, and if so see if we 
  //  // need to initialise it again, and do so if required
  //  if ( fscale_factor!=1 ) {
  //    if ( hoppet == NULL ) hoppet = new hoppet_init();
  //    bool newpdf = hoppet->compareCache(pdf);
  //    //   if ( newpdf ) hoppet->fillCache( pdf );
  //  }
#endif

  //  std::cout << "grid::run() " << run() << std::endl; 

  //  genpdf = genpdf_map[m_genpdfname];
    
  //  TH1D* h = new TH1D(*m_obs_bins);
  //  h->SetName("xsec");

  std::vector<double> hvec;

  string label;

  int lo_order = m_leading_order;
  if ( nloops>=m_order ) { 
    cerr << "too many loops for grid nloops=" << nloops << "\tgrid=" << m_order << endl;   
    return hvec;
  } 

  for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  

    double dsigma = 0;
   
    if ( nloops==0 ) {
      label = "lo      ";
      //      cout << "convolute() nloop=0" << iobs << endl;
      // leading order cross section
      dsigma = m_grids[0][iobs]->convolute_subproc(subproc, pdf, m_genpdf[0], alphas, lo_order, 0, 1, 1, Escale);
    }
    else if ( nloops==1 ) { 
      label = "nlo     ";
      // next to leading order cross section
      // leading and next to order contributions and scale dependent born dependent terms
      double dsigma_lo  = m_grids[0][iobs]->convolute_subproc(subproc, pdf, m_genpdf[0], alphas, lo_order,   1, rscale_factor, 1,  Escale );
      double dsigma_nlo = m_grids[1][iobs]->convolute_subproc(subproc, pdf, m_genpdf[1], alphas, lo_order+1, 0, rscale_factor, 1,  Escale );
      dsigma = dsigma_lo + dsigma_nlo;
    }
    else if ( nloops==-1 ) { 
      label = "nlo only";
      // nlo contribution only (only strict nlo contributions)
      dsigma = m_grids[1][iobs]->convolute_subproc(subproc, pdf, m_genpdf[1], alphas, lo_order+1, 0, rscale_factor, 1, Escale );
    }
    else if ( nloops==2 ) { 
      // FIXME: not implemented completely yet 
      return hvec;
      label = "nnlo    ";
      // next to next to leading order contribution 
      // NB: NO scale dependendent parts, so only muR=muF=mu
      double dsigma_lo   = m_grids[0][iobs]->convolute_subproc(subproc, pdf, m_genpdf[0], alphas, lo_order,   0);
      double dsigma_nlo  = m_grids[1][iobs]->convolute_subproc(subproc, pdf, m_genpdf[1], alphas, lo_order+1, 0);
      double dsigma_nnlo = m_grids[2][iobs]->convolute_subproc(subproc, pdf, m_genpdf[2], alphas, lo_order+2, 0);
      dsigma = dsigma_lo + dsigma_nlo + dsigma_nnlo;
    }
    
    
    //   double deltaobs = h->GetBinLowEdge(iobs+2)-h->GetBinLowEdge(iobs+1);
    //   h->SetBinContent(iobs+1, dsigma/(deltaobs));
    //   h->SetBinError(iobs+1, 0);

    double deltaobs = m_obs_bins->GetBinLowEdge(iobs+2)-m_obs_bins->GetBinLowEdge(iobs+1);
    
    hvec.push_back( invNruns*Escale2*dsigma/deltaobs );
    // hvec.push_back( Escale2*dsigma/deltaobs );


    //    cout << "dsigma[" << iobs << "]=" << dsigma/deltaobs << endl;

    
    //    cout << "obs bin " << iobs 
    //         << "\t" <<  h->GetBinLowEdge(iobs+1) << " - " << h->GetBinLowEdge(iobs+2)
    //	       << "\tdsigma=" << dsigma << endl;
  }  // iobs   

  //  double _ctime = appl_timer_stop(_ctimer);
  //  cout << "grid::convolute_subproc(" << subproc << ") " << label << " convolution time=" << _ctime << " ms" << endl;

  if ( getApplyCorrections() ) applyCorrections(hvec);

  return hvec;
}



TH1D* grid::convolute(void (*pdf)(const double& , const double&, double* ), 
		      double (*alphas)(const double& ), 
		      int     nloops, 
		      double  rscale_factor,
		      double  fscale_factor,
		      double Escale ) {

    TH1D* h = new TH1D(*m_obs_bins);
    h->SetName("xsec");
    
    std::vector<double> dvec = vconvolute( pdf, alphas, nloops, rscale_factor, fscale_factor, Escale );
    
    for ( unsigned i=0 ; i<dvec.size() ; i++ ) { 
      h->SetBinContent( i+1, dvec[i] );
      h->SetBinError( i+1, 0 );
    }
    
    return h;

}





TH1D* grid::convolute_subproc(int subproc,
			      void (*pdf)(const double& , const double&, double* ), 
			      double (*alphas)(const double& ), 
			      int     nloops, 
			      double  rscale_factor, double Escale ) {

    TH1D* h = new TH1D(*m_obs_bins);
    h->SetName("xsec");
    
    std::vector<double> dvec = vconvolute_subproc( subproc, pdf, alphas, nloops, rscale_factor, Escale );
    
    for ( unsigned i=0 ; i<dvec.size() ; i++ ) { 
      h->SetBinContent( i+1, dvec[i] );
      h->SetBinError( i+1, 0 );
    }
  
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
				    oldgrid->transform(), m_genpdf[iorder]->Nproc());

  delete oldgrid;
}
  


void grid::setRange(int ilower, int iupper, double xScaleFactor) { 
  if ( ilower>=0 && iupper <Nobs() ) {  
    double lower = getReference()->GetBinLowEdge(ilower+1);
    double upper = getReference()->GetBinLowEdge(iupper+2); 
    setRange( lower, upper, xScaleFactor );
  }
}


void grid::setRange(double lower, double upper, double xScaleFactor) { 
  
  std::cout << "grid::SetRange() " << lower << " " << upper << std::endl; 

  std::vector<bool>   used;
  std::vector<double> limits;
  std::vector<double> contents;
  std::vector<double> errors;

  /// get the occupied bins
  //  int Nbins = 0;
  double last = 0;
  for ( int i=1 ; i<=m_obs_bins->GetNbinsX() ; i++ ) { 
    double bin =  m_obs_bins->GetBinCenter(i);
    if ( bin>lower && bin<upper ) { 
      limits.push_back( m_obs_bins->GetBinLowEdge(i) );
      contents.push_back( m_obs_bins->GetBinContent(i) );
      errors.push_back( m_obs_bins->GetBinError(i) );

      last = m_obs_bins->GetBinLowEdge(i+1);
      used.push_back(true);
    }
    else { 
      used.push_back(false);
    }
  }
  
  /// copy the range of the reference histogram
  if ( limits.size()>0 ) limits.push_back( last );
  else { 
    throw grid::exception( std::cerr << "new range does not include any bins"  ); 
  }

  if ( xScaleFactor!=1 ) { 
    for ( unsigned i=0 ; i<limits.size(); i++ ) limits[i] *= xScaleFactor;
  }

  /// save bins somewhere so can overwrite them 
  TH1D* h = m_obs_bins;
  
  m_obs_bins = new TH1D( h->GetTitle(), h->GetName(), limits.size()-1, &(limits[0]) );
  
  for ( int i=0 ; i<m_obs_bins->GetNbinsX() ; i++ ) { 
    m_obs_bins->SetBinContent( i+1, contents[i] );
    m_obs_bins->SetBinError( i+1, errors[i] );
  }

  /// copy the igrids for the observable bins in the range 

  igrid** grids[3];

  /// save old grids
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) grids[iorder] = m_grids[iorder];
  
  int Nobs = m_obs_bins->GetNbinsX();

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    m_grids[iorder] = new igrid*[Nobs];
    int iobs = 0;
    for ( int igrid=0 ; igrid<h->GetNbinsX() ; igrid++ ) {
      if ( used[igrid] ) m_grids[iorder][iobs++] = grids[iorder][igrid];
      else               delete grids[iorder][igrid];                           
    }
  }

  delete h;

}


/// methods to handle the documentation
void grid::setDocumentation(const std::string& s) { m_documentation = s; }
void grid::addDocumentation(const std::string& s) {   
  if ( m_documentation.size() ) m_documentation += s;
  else                          setDocumentation(s);    
}






/// methods to handle bin-by-bin corrections

/// add a correction as a vector
void grid::addCorrection( std::vector<double>& v, const std::string& label) {
  //  std::cout << "addCorrections(vector) " << v.size() << " " << m_obs_bins->GetNbinsX() << std::endl;
  if ( v.size()==unsigned(m_obs_bins->GetNbinsX()) ) {
    m_corrections.push_back(v);
    m_correctionLabels.push_back(label);
    //  std::cout << "appl::grid::addCorrection(vector) now " << m_corrections.size() << " corrections" << std::endl;
  }
}


/// add a correction by histogram
void grid::addCorrection(TH1D* h, const std::string& label) {
  // std::cout << "addCorrections(TH1D*) " << h->GetNbinsX() << " " << m_obs_bins->GetNbinsX() << std::endl;
  if ( h->GetNbinsX()==m_obs_bins->GetNbinsX() ) {
    for ( int i=1 ; i<=h->GetNbinsX()+1 ; i++ ) { 
      if ( std::fabs(h->GetBinLowEdge(i+1)-m_obs_bins->GetBinLowEdge(i+1))>1e-10 ) { 
	std::cerr << "grid::addCorrection(TH1D* h): bin mismatch, not adding correction" << std::endl;
	return;
      }
    }

    std::vector<double> v(h->GetNbinsX());
    for ( int i=0 ; i<h->GetNbinsX() ; i++ ) v[i] = h->GetBinContent(i+1);
    if ( label=="" ) addCorrection(v, h->GetName());
    else             addCorrection(v, label);
  }
  else { 
    std::cerr << "grid::addCorrection(TH1D* h): bin mismatch, not adding correction" << std::endl;
  }
  
}




/// apply corrections to a vector
void grid::applyCorrections(std::vector<double>& v) {
  //  std::cout << "grid::applyCorrections(vector) " << m_corrections.size() << std::endl;
  for ( unsigned i=0 ; i<m_corrections.size() ; i++ ) { 
    std::vector<double>& correction = m_corrections[i];
    //      TH1D* hc = m_corrections[i];
    for ( unsigned j=0 ; j<v.size() ; j++ ) v[j] *= correction[j];
  }
  //  std::cout << "grid::applyCorrections(vector) done" << std::endl;
}




ostream& operator<<(ostream& s, const appl::grid& g) {
  s << "==================================================" << endl;
  //  s << "appl::grid version " << g.version() << "\t(" << g.subProcesses(0) << " initial states, " << g.Nobs() << " observable bins)" << endl;

  std::string order[3] = {  "-LO, ",  "-NLO, ",  "-NNLO, " };  

  s << "appl::grid version " << g.version() << "\t( "; 
  for ( int i=0 ; i<g.nloops()+1 ; i++ ) s << g.subProcesses(i) << order[i];
  s << "initial states, " << g.Nobs() << " observable bins )" << endl;
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
