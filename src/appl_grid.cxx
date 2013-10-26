
//  appl_grid.cxx        

//   grid class - all the functions needed to create and 
//   fill the grid from an NLO calculation program. 
//  
//  Copyright (C) 2007 Mark Sutton (sutt@hep.ucl.ac.uk)    

// $Id: appl_grid.cxx, v1.00 2007/10/16 17:01:39 sutt $

#include <cstdlib>
#include <stdio.h>
#include <sys/stat.h>

#include <iostream>
#include <iomanip>
#include <cmath>

#include "appl_grid/appl_pdf.h"
#include "appl_grid/appl_timer.h"
#include "appl_grid/Directory.h"
#include "appl_grid/appl_grid.h"

#include "appl_grid/generic_pdf.h"
#include "appl_grid/lumi_pdf.h"

#include "appl_igrid.h"
#include "Cache.h"


#include "TFileString.h"
#include "TFileVector.h"

#include "TFile.h"
#include "TObjString.h"
#include "TVectorT.h"


#include "amconfig.h"



/// this is a compatability flag for persistent versions 
/// of the grid
/// NB: ONLY change this if the persistent class
///     changes in a non-backwards compatible way.

const std::string appl::grid::m_version = "version-3.2";

std::string appl::grid::appl_version() const { return PACKAGE_VERSION; }

#include "hoppet_init.h"

#ifdef HAVE_HOPPET

#include "hoppet_v1.h"

// include hoppet splitting function code

static hoppet_init* hoppet = 0;

void Splitting(const double& x, const double& Q, double* xf) {
  static const int nLoops    = 1;
  static const int nFlavours = 5;
  hoppetevalsplit_( x, Q, nLoops, nFlavours, xf); 
  return;
}

#else

void Splitting(const double& x, const double& Q, double* xf) {
  throw appl::grid::exception( std::cerr << "hoppet library not included - cannot call splitting function"  ); 
  return; // technically, don't need this - should throw an exception
}

#endif




/// helper function
static bool contains(const std::string& s, const std::string& reg ) { 
  return s.find(reg)!=std::string::npos;
}

/// make sure pdf std::map is initialised
// bool pdf_ready = appl::appl_pdf::create_map(); 


appl::grid::grid(int NQ2, double Q2min, double Q2max, int Q2order, 
		 int Nx,  double xmin,  double xmax,  int xorder,
		 int Nobs,  double obsmin, double obsmax, 
		 std::string genpdfname,
		 int leading_order, int nloops, 
		 std::string transform ) :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_normalised(false), m_symmetrise(false), 
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),
  m_documentation(""),
  m_type(STANDARD)
{
  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", Nobs, obsmin, obsmax);
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Sumw2(); /// grrr root is so rubbish - not scaling errors properly

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  //  if      ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  //  else if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  if ( contains(m_genpdfname, ".dat") ||  contains(m_genpdfname, ".config") ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform); 
}




appl::grid::grid(int Nobs, const double* obsbins, 
		 int NQ2,  double Q2min, double Q2max, int Q2order, 
		 int Nx,   double xmin, double xmax,   int xorder, 
		 std::string genpdfname, 
		 int leading_order, int nloops, 
		 std::string transform ) :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false),  m_normalised(false), m_symmetrise(false),
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),
  m_documentation(""),
  m_type(STANDARD)
{
  
  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", Nobs, obsbins);
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Sumw2(); /// grrr root is so rubbish - not scaling errors properly

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  //  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  if ( contains(m_genpdfname, ".dat") ||  contains(m_genpdfname, ".config") ) addpdf(m_genpdfname);

  findgenpdf( m_genpdfname );

  construct(Nobs, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform );
}




appl::grid::grid(const std::vector<double>& obs, 
		 int NQ2, double Q2min, double Q2max, int Q2order,
		 int Nx,  double xmin,  double xmax,  int xorder, 
		 std::string genpdfname,
		 int leading_order, int nloops, 
		 std::string transform )  :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_normalised(false), m_symmetrise(false),  
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),
  m_documentation(""),
  m_type(STANDARD)
{
  
  if ( obs.size()==0 ) { 
    std::cerr << "grid::not enough bins in observable" << std::endl;
    std::exit(0);
  } 
  
  //  double* obsbins = new double[obs.size()];  
  //  for ( unsigned i=0 ; i<obs.size() ; i++ ) obsbins[i] = obs[i];
  //  int Nobs = obs.size()-1;

  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", obs.size()-1, &obs[0] );
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Sumw2(); /// grrr root is so rubbish - not scaling errors properly
  //  delete[] obsbins;

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  //   if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  if ( contains(m_genpdfname, ".dat") ||  contains(m_genpdfname, ".config") ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  construct( obs.size()-1, NQ2, Q2min, Q2max, Q2order, Nx, xmin, xmax, xorder, m_order, m_transform); 
}



appl::grid::grid(const std::vector<double>& obs, 
		 std::string genpdfname,
		 int leading_order, int nloops, 
		 std::string transform )  :
  m_leading_order(leading_order), m_order(nloops+1), 
  m_run(0), m_optimised(false), m_trimmed(false), m_normalised(false), m_symmetrise(false),  
  m_transform(transform), m_genpdfname(genpdfname), m_cmsScale(0),
  m_applyCorrections(false),  
  m_documentation(""),
  m_type(STANDARD)
{ 

  if ( obs.size()==0 ) { 
    std::cerr << "grid::not enough bins in observable" << std::endl;
    std::exit(0);
  } 
  
  //  double* obsbins = new double[obs.size()];  
  //  for ( unsigned i=0 ; i<obs.size() ; i++ ) obsbins[i] = obs[i];
  //  int Nobs = obs.size()-1;

  // Initialize histogram that saves the correspondence obsvalue<->obsbin
  m_obs_bins=new TH1D("referenceInternal","Bin-Info for Observable", obs.size()-1, &obs[0] );
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Sumw2(); /// grrr root is so rubbish - not scaling errors properly
  //  delete[] obsbins;

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  //  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  if ( contains(m_genpdfname, ".dat") ||  contains(m_genpdfname, ".config") ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) m_grids[iorder] = new igrid*[obs.size()-1];

}



appl::grid::grid(const std::string& filename, const std::string& dirname)  :
  m_leading_order(0),  m_order(0),
  m_optimised(false),  m_trimmed(false), 
  m_normalised(false),
  m_symmetrise(false), m_transform(""), 
  m_applyCorrections(false),
  m_documentation(""),
  m_type(STANDARD)
{
  
  struct stat _fileinfo;
  if ( stat(filename.c_str(),&_fileinfo) )   {    
    throw exception(std::cerr << "grid::grid() cannot open file " << filename << std::endl ); 
  }

  std::cout << "appl::grid() reading grid from file " << filename << std::endl;
  
  TFile* gridfilep = TFile::Open(filename.c_str());
  
  if (gridfilep->IsZombie()) {
    throw exception(std::cerr << "grid::grid() cannot open file: zombie " << filename << std::endl ); 
    delete gridfilep;
  }

  // TFile gridfile(filename.c_str());
  
  //  gDirectory->cd(dirname.c_str());

  //  std::cout << "pwd=" << gDirectory->GetName() << std::endl;

  //  gDirectory->cd(dirname.c_str());

  //  std::cout << "pwd=" << gDirectory->GetName() << std::endl;

  //  Directory d(dirname);
  //  d.push();

  TFileString _tags = *(TFileString*)gridfilep->Get((dirname+"/Tags").c_str());  
  // TFileString _tags = *(TFileString*)gridfile.Get("Tags");  
  m_transform  = _tags[0];
  m_genpdfname = _tags[1];

  std::string _version = _tags[2];

  //  std::cout << "tags:: transform: " << m_transform << "\tpdfname: " << m_genpdfname << "\tdoc: " << m_documentation << std::endl;  

  if ( _tags.size()>3 ) m_documentation = _tags[3];

  // check it has the correct version
  // if ( _version > m_version ) { 
  //      throw exception(cerr << "incorrect version " << _version << " expected " << m_version ); 
  // }
  //  m_version = _version;
  
  //  std::cout << "appl::grid " << m_version << "\t" << m_documentation << std::endl; 
  
  //  std::cout << "Tags=" << _tags << std::endl;

  //  std::cout << "grid::grid() read transform " << m_transform << " from file" << std::endl;

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

  //  std::vector<double> _ckmsum;
  std::vector<std::vector<double> > _ckm2;
  std::vector<std::vector<double> > _ckm;




  if ( setup->GetNoElements()>8 && (*setup)(8)!=0 ) {

    std::cout << "grid::grid() read ckm matrices" << std::endl;
    
    /// try 13x13 squared ckm matrix 

    TVectorT<double>* ckm2flat=(TVectorT<double>*)gridfilep->Get((dirname+"/CKM2").c_str());

    if ( ckm2flat ) { 
      if ( ckm2flat->GetNrows()>0 ) { 
	_ckm2 = std::vector<std::vector<double> >(13, std::vector<double>(13) );
      
	for ( int ic=0 ; ic<13 ; ic++ ) { 
	  for ( int id=0 ; id<13 ; id++ ) _ckm2[ic][id] = (*ckm2flat)(ic*13+id); 
	}
      }  
    }

    /// now try usual 3x3 matrix

    TVectorT<double>* ckmflat=(TVectorT<double>*)gridfilep->Get((dirname+"/CKM").c_str());

    if ( ckmflat ) { 
      if ( ckmflat->GetNrows()>0 ) { 
	_ckm = std::vector<std::vector<double> >(3, std::vector<double>(3) );
	
	for ( int ic=0 ; ic<3 ; ic++ ) { 
	  for ( int id=0 ; id<3 ; id++ ) _ckm[ic][id] = (*ckmflat)(ic*3+id); 
	}
      }  
    }

  }

  if ( setup->GetNoElements()>9 ) m_type = (CALCULATION)int( (*setup)(9)+0.5 );
  else                            m_type = STANDARD;

  std::cout << "appl::grid() reading grid calculation type: " << _calculation(m_type) << std::endl;

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf (or lumi_pdf for amcatnlo)
  //  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  
  if ( contains(m_genpdfname, ".config") ) { 
    /// decode the pdf combination if appropriate

    /// again have to use TVectorT<double> because TVectorT<int> has no constructor!!!
    /// I ask you!! what's the point of a template if it doesn't actually instantiate
    /// it's pathetic!
    TVectorT<double>* _combinations = (TVectorT<double>*)gridfilep->Get((dirname+"/Combinations").c_str());

    std::vector<int> combinations(_combinations->GetNoElements());

    for ( unsigned ic=0 ; ic<combinations.size() ; ic++ ) {
      combinations[ic] = int((*_combinations)(ic)); 
    }
    
    addpdf(m_genpdfname, combinations);
  }
  else { 
    /// of just create the generic from the file
    if ( contains(m_genpdfname, ".dat") ) addpdf(m_genpdfname);
  }

  /// returieve the pdf routine 
  findgenpdf( m_genpdfname );

  //  std::cout << "grid::grid() read " << m_genpdfname << " " << m_genpdf[0]->getckmsum().size() << std::endl; 

  // set the ckm matrices 
  if      ( _ckm.size()>0 )  setckm( _ckm );
  else if ( _ckm2.size()>0 ) setckm2( _ckm2 );

  delete setup;

  //  std::cout << "grid::grid() read setup" << std::endl;

  // Read observable bins information
  //  gridfile.GetObject("obs_bins", m_obs_bins );
  m_obs_bins = (TH1D*)gridfilep->Get((dirname+"/reference").c_str());
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Scale(run());
  m_obs_bins->SetName("referenceInternal");


  //  std::cout << "grid::grid() read obs bins" << std::endl;

  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    //  std::cout << "grid::grid() iorder=" << iorder << std::endl;
    m_grids[iorder] = new igrid*[Nobs()];  
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {
      char name[128];  sprintf(name, (dirname+"/weight[alpha-%d][%03d]").c_str(), iorder, iobs);
      // std::cout << "grid::grid() reading " << name << "\tiobs=" << iobs << std::endl;

      m_grids[iorder][iobs] = new igrid(*gridfilep, name);

      //    _size += m_grids[iorder][iobs]->size();
      //      std::cout << "grid::grid() done" << std::endl;
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

  //  std::cout << "grid::grid() read from file Nobs = " << Nobs() << std::endl;

  //  std::cout << "read grid" << std::endl;

  gridfilep->Close();
  delete gridfilep;
}


appl::grid::grid(const grid& g) : 
  m_obs_bins(new TH1D(*g.m_obs_bins)), 
  m_leading_order(g.m_leading_order), m_order(g.m_order), 
  m_run(g.m_run), m_optimised(g.m_optimised), m_trimmed(g.m_trimmed), 
  m_normalised(true),
  m_symmetrise(g.m_symmetrise),
  m_transform(g.m_transform),
  m_genpdfname(g.m_genpdfname), 
  m_cmsScale(g.m_cmsScale),
  m_applyCorrections(g.m_applyCorrections), 
  m_documentation(g.m_documentation),
  m_ckmsum(g.m_ckmsum), /// need a deep copy of the contents
  m_ckm2(g.m_ckm2),     /// need a deep copy of the contents
  m_ckm(g.m_ckm),       /// need a deep copy of the contents
  m_type(g.m_type)
{
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Sumw2();

  /// check to see if we require a generic pdf from a text file, and 
  /// and if so, create the required generic pdf
  //  if ( m_genpdfname.find(".dat")!=std::string::npos ) addpdf(m_genpdfname);
  if ( contains(m_genpdfname, ".dat") ||  contains(m_genpdfname, ".config") ) addpdf(m_genpdfname);
  findgenpdf( m_genpdfname );

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    m_grids[iorder] = new igrid*[Nobs()];
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  m_grids[iorder][iobs] = new igrid(*g.m_grids[iorder][iobs]);
  }
} 


// Initialize histogram that saves the correspondence obsvalue<->obsbin

// constructor common internals 
void appl::grid::construct(int Nobs, 
			   int NQ2,  double Q2min, double Q2max, int Q2order,
			   int Nx,   double xmin,  double xmax,  int xorder, 
			   int order, 
			   std::string transform ) { 
  
  //  std::cout << "appl::grid::construct() m_order " << m_order << "\tNobs " << Nobs << std::endl; 

  if ( m_order!=order ) std::cerr << "appl::grid::construct() order mismatch" << std::endl;

  for ( int iorder=0 ; iorder<m_order ; iorder++ ) m_grids[iorder] = 0;

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
int appl::grid::subProcesses(int i) const { 
  if ( i<0 || i>=m_order ) throw exception( std::cerr << "grid::subProcess(int i) " << i << " out or range [0-" << m_order-1 << "]" << std::endl );
  return m_grids[i][0]->SubProcesses();     
}  


/// access the transform functions for the appl::igrid so that the 
/// igrid can be hidden 
double appl::grid::transformvar()         { return igrid::transformvar(); }
double appl::grid::transformvar(double v) { return igrid::transformvar(v); }



// add a single grid
void appl::grid::add_igrid(int bin, int order, igrid* g) { 

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




appl::grid::~grid() {
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

appl::grid& appl::grid::operator=(const appl::grid& g) { 
  // clear out the old...
  delete m_obs_bins;
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  delete m_grids[iorder][iobs];
    delete m_grids[iorder];
  }
  
  // copy the new
  m_obs_bins = new TH1D(*g.m_obs_bins);
  m_obs_bins->SetDirectory(0);
  m_obs_bins->Sumw2();

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
  

appl::grid& appl::grid::operator*=(const double& d) { 
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) (*m_grids[iorder][iobs])*=d; 
  }
  getReference()->Scale( d );
  return *this;
}



appl::grid& appl::grid::operator+=(const appl::grid& g) {
  m_run      += g.m_run;
  m_optimised = g.m_optimised;
  m_trimmed   = g.m_trimmed;
  if ( Nobs()!=g.Nobs() )   throw exception("grid::operator+ Nobs bin mismatch");
  if ( m_order!=g.m_order ) throw exception("grid::operator+ different order grids");
  if ( m_leading_order!=g.m_leading_order ) throw exception("grid::operator+ different order processes in grids");
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) (*m_grids[iorder][iobs]) += (*g.m_grids[iorder][iobs]); 
  }

  /// grrr use root TH1::Add() even though I don't like it. 
  getReference()->Add( g.getReference() );

  return *this;
}




/// check grids match
bool appl::grid::operator==(const appl::grid& g) const {
  
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
void appl::grid::fill(const double x1, const double x2, const double Q2, 
		      const double obs, 
		      const double* weight, const int iorder)  {  
  int iobs = m_obs_bins->FindBin(obs)-1;
  if ( iobs<0 || iobs>=Nobs() ) {
    //    cerr << "grid::fill() iobs out of range " << iobs << "\tobs=" << obs << std::endl;
    //    cerr << "obs=" << obs << "\tobsmin=" << obsmin() << "\tobsmax=" << obsmax() << std::endl;
    return;
  }
  
  //  std::cout << "iobs=" << iobs << "\tobs=" << obs;
  //  for ( int i=0 ; i<subProcesses(iorder) ; i++ ) std::cout << "\t" << weight[i];
  //  std::cout << std::endl;

  //  std::cout << "\tiobs=" << iobs << std::endl;
  if ( m_symmetrise && x2<x1 )  m_grids[iorder][iobs]->fill(x2, x1, Q2, weight);
  else                          m_grids[iorder][iobs]->fill(x1, x2, Q2, weight);
}


// fast fill pre-optimisation don't perform the interpolation and so on
void appl::grid::fill_phasespace(const double x1, const double x2, const double Q2, 
				 const double obs, 
				 const double* weight, const int iorder) {
  int iobs = m_obs_bins->FindBin(obs)-1;
  if ( iobs<0 || iobs>=Nobs() ) {
    //  cerr << "grid::fill() iobs out of range " << iobs << "\tobs=" << obs << std::endl;
    //  cerr << "obs=" << obs << "\tobsmin=" << obsmin() << "\tobsmax=" << obsmax() << std::endl;
    return;
  }
  if ( m_symmetrise && x2<x1 )  m_grids[iorder][iobs]->fill_phasespace(x2, x1, Q2, weight);
  else                          m_grids[iorder][iobs]->fill_phasespace(x1, x2, Q2, weight);
}




// fast fill pre-optimisation don't perform the interpolation and so on
void appl::grid::fill_index(const int ix1, const int ix2, const int iQ2, 
			    const int iobs, 
			    const double* weight, const int iorder) {

  if ( iobs<0 || iobs>=Nobs() ) {
    //  cerr << "grid::fill() iobs out of range " << iobs << "\tobs=" << obs << std::endl;
    //  cerr << "obs=" << obs << "\tobsmin=" << obsmin() << "\tobsmax=" << obsmax() << std::endl;
    return;
  }
  if ( m_symmetrise && ix2<ix1 )  m_grids[iorder][iobs]->fill_index(ix2, ix1, iQ2, weight);
  else                            m_grids[iorder][iobs]->fill_index(ix1, ix2, iQ2, weight);
}


void appl::grid::trim() {
  m_trimmed = true;
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) m_grids[iorder][iobs]->trim(); 
  }
}

void appl::grid::untrim() {
  m_trimmed = false;
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) m_grids[iorder][iobs]->untrim(); 
  }
}

std::ostream& appl::grid::print(std::ostream& s) const {
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {     
      s << iobs << "\t" 
	<< std::setprecision(5) << std::setw(6) << getReference()->GetBinLowEdge(iobs+1) << "\t- " 
	<< std::setprecision(5) << std::setw(6) << getReference()->GetBinLowEdge(iobs+2) << "\t"; 
      m_grids[iorder][iobs]->print(s);       
    }
  }
  return s;
}




 /// get the required pdf combinations from those registered   
void appl::grid::findgenpdf( std::string s ) { 
    std::vector<std::string> names = parse( s, ":" );
    if ( names.size()==unsigned(m_order) ) for ( int i=0 ; i<m_order ; i++ ) m_genpdf[i] = appl_pdf::getpdf( names[i] );
    else  if ( names.size()==1 )           for ( int i=0 ; i<m_order ; i++ ) m_genpdf[i] = appl_pdf::getpdf( names[0] );
    else  { 
      throw exception( std::cerr << "requested " << m_order << " pdf combination but given " << names.size() << std::endl );
    }
}


void appl::grid::addpdf( const std::string& s, const std::vector<int>& combinations ) {

  //  std::cout << "addpdf() in " << std::endl;

    /// parse names, if they contain .dat, then create the new generic pdfs
    /// they will be added to the pdf std::map automatically 
    std::vector<std::string> names = parse( s, ":" );

    unsigned imax = unsigned(m_order); 

    /// check to see whether we have a different pdf for each order
    if ( names.size()!=imax ) { 
      if ( names.size()==1 ) imax = 1;
      else { 
	throw exception( std::cerr << "requested " << m_order << " pdf combination but given " << names.size() << std::endl );
      }
    }

    //    std::cout << "imax " << imax << std::endl; 

    /// loop through all the required pdfs checking whether they exist already,
    /// if not (from thrown exception) then create it, otherwise, don't need to 
    /// do anything 
    for ( unsigned i=0 ; i<imax ; i++ ) { 

      //      std::cout << "\ti " << i<< std::endl; 

      if ( names[i].find(".dat")!=std::string::npos ) { 
	try {
	  appl_pdf::getpdf(names[i]); // , false);
	}
	catch ( appl_pdf::exception e ) { 
	  std::cout << "creating new generic_pdf " << names[i] << std::endl;
	  new generic_pdf(names[i]);
	}
      }
      else if ( names[i].find(".config")!=std::string::npos ) { 
	try {
	  appl_pdf::getpdf(names[i]); // , false);
	}
	catch ( appl_pdf::exception e ) { 
	  std::cout << "creating new lumi_pdf " << names[i] << std::endl;
	  new lumi_pdf(names[i], combinations);
	  //	  std::cout << "created" << names[i] << std::endl;

	}
      }

    }

}



void appl::grid::setckm2( const std::vector<std::vector<double> >& ckm2 ) { 
  for ( int i=0 ; i<m_order ; i++ ) m_genpdf[i]->setckm2(ckm2);
}


void appl::grid::setckm( const std::vector<std::vector<double> >& ckm ) { 
  for ( int i=0 ; i<m_order ; i++ ) m_genpdf[i]->setckm(ckm);
}


void appl::grid::setckm( const std::vector<double>& ckm ) {
  std::vector<std::vector<double> > _ckm(3, std::vector<double>(3,0) );
  for ( unsigned i=0 ; i<ckm.size() && i<9 ; i++ ) _ckm[i/3][i%3] = ckm[i]; 
  for ( int i=0 ; i<m_order ; i++ ) m_genpdf[i]->setckm(_ckm);
}

void appl::grid::setckm( const double* ckm ) { 
  std::vector<std::vector<double> > _ckm(3, std::vector<double>(3,0) );
  for ( unsigned i=0 ; i<9 ; i++ ) _ckm[i/3][i%3] = ckm[i]; 
  for ( int i=0 ; i<m_order ; i++ ) m_genpdf[i]->setckm(_ckm);
}


const std::vector<std::vector<double> >& appl::grid::getckm()  const { return m_genpdf[0]->getckm(); }  

const std::vector<std::vector<double> >& appl::grid::getckm2() const { return m_genpdf[0]->getckm2(); }  


// void appl::grid::setuppdf(void (*pdf)(const double&, const double&, double* ) )  {  }
// void grid::pdfinterp(double x, double Q2, double* f) {  }



// set the rewight flag of the internal grids
bool appl::grid::reweight(bool t) { 
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) m_grids[iorder][iobs]->reweight(t);       
  }
  return t;
}



// dump to file
void appl::grid::Write(const std::string& filename, const std::string& dirname) { 
 

  struct stat sb;
  
  if ( stat( filename.c_str(), &sb)==0 ) { // && S_ISREG(sb.st_mode ))
    std::string filename_save = filename + "-save";
    if ( std::rename( filename.c_str(), filename_save.c_str() ) ) std::cerr << "could not rename grid file " << filename << std::endl;
  } 


  //  std::cout << "grid::Write() writing to file " << filename << std::endl;
  TFile rootfile(filename.c_str(),"recreate");

  //  std::cout << "pwd=" << gDirectory->GetName() << std::endl;

  Directory d(dirname);
  d.push();
  
  //  std::cout << "pwd=" << gDirectory->GetName() << std::endl;

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


  //  std::cout << "state std::vector=" << std::endl;

  // state information
  TVectorT<double>* setup=new TVectorT<double>(12); // add a few extra just in case 
  (*setup)(0) = m_run;
  (*setup)(1) = ( m_optimised  ? 1 : 0 );
  (*setup)(2) = ( m_symmetrise ? 1 : 0 );
  (*setup)(3) =   m_leading_order ;
  (*setup)(4) =   m_order ;
  (*setup)(5) =   m_cmsScale ;
  (*setup)(6) = ( m_normalised ? 1 : 0 );
  (*setup)(7) = ( m_applyCorrections ? 1 : 0 );

  if ( m_genpdf[0]->getckmsum().size()==0 ) (*setup)(8) = 0;
  else                                      (*setup)(8) = 1;

  (*setup)(9) = (int)m_type;

  setup->Write("State");
  
  if ( (*setup)(8) == 1 ) { 
    
#if 0
    /// no longer write out squared ckm matrix - just use the 3x3
    TVectorT<double>* ckm2flat = new TVectorT<double>(169);
    const std::vector<std::vector<double> >& _ckm2 = m_genpdf[0]->getckm2();

    for ( int ic=0 ; ic<13 ; ic++ ) { 
      for ( int id=0 ; id<13 ; id++ ) (*ckm2flat)(ic*13+id) = _ckm2[ic][id];
    }

    ckm2flat->Write("CKM2");
#endif

    /// no longer write out squared ckm matrix - just use the 3x3
    TVectorT<double>* ckmflat = new TVectorT<double>(9);
    const std::vector<std::vector<double> >& _ckm = m_genpdf[0]->getckm();

    for ( int ic=0 ; ic<3 ; ic++ ) { 
      for ( int id=0 ; id<3 ; id++ ) (*ckmflat)(ic*3+id) = _ckm[ic][id];
    }

    ckmflat->Write("CKM");

  }


  /// encode the pdf combination if appropriate

  if ( contains( m_genpdfname, ".config" ) ) { 
    std::vector<int>   combinations = dynamic_cast<lumi_pdf*>(m_genpdf[0])->serialise();  
    TVectorT<double>* _combinations = new TVectorT<double>(combinations.size()); // add a few extra just in case 
    for ( unsigned ic=0 ; ic<combinations.size() ; ic++ ) { 
      //     std::cout << "write " << ic << "\tcombinations " << combinations[ic] << std::endl;
      /// because root stupidly doesn't have a constructor for TVectorT<int> 
      /// we need to store these integers as doubles - this mean we add (or subtract) 
      /// 0.5 from each value to ensure that the double->int conversion doesn't mess 
      /// up with 0.9999 -> 0 etc    
      if ( combinations[ic]<0 ) (*_combinations)(ic) = double(combinations[ic]-0.5);
      else                      (*_combinations)(ic) = double(combinations[ic]+0.5);
    }
    _combinations->Write("Combinations");
  }

  

  //  int _size     = 0;
  //  int trim_size = 0;

  //  std::cout << "grids Nobs = " << Nobs() << std::endl;

  // internal grids
  for( int iorder=0 ; iorder<m_order ; iorder++ ) {
    for( int iobs=0 ; iobs<Nobs() ; iobs++ ) {
      char name[128];  sprintf(name, "weight[alpha-%d][%03d]", iorder, iobs);
      // std::cout << "writing grid " << name << std::endl;
      //   _size += m_grids[iorder][iobs]->size();
      m_grids[iorder][iobs]->write(name);
      //   trim_size += m_grids[iorder][iobs]->size();
    }
  }
  //  std::cout <<"grid::Write() size(untrimmed)=" << _size 
  //     << "\tsize(trimmed)="              << trim_size << std::endl;
  //  d.pop();

  //  std::cout << "reference" << std::endl;
  
  TH1D* reference = (TH1D*)m_obs_bins->Clone("reference");
  
  if ( !getNormalised() )  if ( run() ) reference->Scale(1/double(run()));
  
  reference->Write();
  delete reference;

  //  std::cout << "corrections" << std::endl;

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

  d.pop();
  rootfile.Close();

  //  std::cout << "written" << std::endl;
}



// takes pdf as the pdf lib wrapper for the pdf set for the convolution.
// type specifies which sort of partons should be included:

std::vector<double> appl::grid::vconvolute(void (*pdf)(const double& , const double&, double* ), 
					   double (*alphas)(const double& ), 
					   int     nloops, 
					   double  rscale_factor,
					   double  fscale_factor,
					   double Escale )
{ 
  return vconvolute( pdf, 0, alphas, nloops, rscale_factor, fscale_factor, Escale );
}

std::vector<double> appl::grid::vconvolute(void (*pdf1)(const double& , const double&, double* ), 
					   void (*pdf2)(const double& , const double&, double* ), 
					   double (*alphas)(const double& ), 
					   int     nloops, 
					   double  rscale_factor,
					   double  fscale_factor,
					   double Escale )
{ 

  NodeCache cache1( pdf1 );
  NodeCache cache2;

  cache1.reset();

  NodeCache* _pdf1 = &cache1;
  NodeCache* _pdf2 = 0;
  
  if ( pdf2!=0 && pdf1!=pdf2 ) { 
    cache2 = NodeCache( pdf2 );
    cache2.reset();
  
    _pdf2    = &cache2;
  }

 


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

    if ( pdf2==0 || pdf1==pdf2 ) { 

      if ( hoppet == 0 ) hoppet = new hoppet_init();
 
      bool newpdf = hoppet->compareCache(pdf1);
      
      if ( newpdf ) hoppet->fillCache( pdf1 );

    }

  }
#endif

  std::string label;
  
  if ( nloops>=m_order ) { 
    std::cerr << "too many loops for grid nloops=" << nloops << "\tgrid=" << m_order << std::endl;   
    return hvec;
  } 
  
  
 
  if ( m_type==STANDARD ) { 

    //    std::cout << "standard convolution" << std::endl;

    for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  

      double dsigma = 0;
     
      if ( nloops==0 ) {
	label = "lo      ";
	// leading order cross section
	dsigma = m_grids[0][iobs]->convolute( _pdf1, _pdf2, m_genpdf[0], alphas, m_leading_order, 0, 1, 1, Escale);
      }
      else if ( nloops==1 ) { 
	label = "nlo     ";
	// next to leading order cross section
	// std::cout << "convolute() nloop=1" << std::endl;
	// leading order contribution and scale dependent born dependent terms
	double dsigma_lo  = m_grids[0][iobs]->convolute( _pdf1, _pdf2, m_genpdf[0], alphas, m_leading_order, 1, rscale_factor, fscale_factor, Escale);
	// std::cout << "dsigma_lo=" << dsigma_lo << std::endl;
	// next to leading order contribution
	//      double dsigma_nlo = m_grids[1][iobs]->convolute(pdf, m_genpdf, alphas, m_leading_order+1, 0);
	// GPS: the NLO piece must use the same rscale_factor and fscale_factor as
	//      the LO piece -- that's the convention that defines how NLO calculations
	//      are done.
	double dsigma_nlo = m_grids[1][iobs]->convolute( _pdf1, _pdf2, m_genpdf[1], alphas, m_leading_order+1, 0, rscale_factor, fscale_factor, Escale);
	// std::cout << "dsigma_nlo=" << dsigma_nlo << std::endl;
	dsigma = dsigma_lo + dsigma_nlo;
      }
      else if ( nloops==-1 ) {
	label = "nlo only";
	// nlo contribution only (only strict nlo contributions) 
	dsigma = m_grids[1][iobs]->convolute( _pdf1, _pdf2, m_genpdf[1], alphas, m_leading_order+1, 0, rscale_factor, fscale_factor, Escale);
      } 
      else if ( nloops==2 ) {
	// FIXME: not implemented completely yet 
	//      return hvec;
	label = "nnlo    ";
	// next to next to leading order contribution 
	// NB: NO scale dependendent parts so only  muR=muF=mu
	double dsigma_lo  = m_grids[0][iobs]->convolute( _pdf1, _pdf2, m_genpdf[0], alphas, m_leading_order, 0);
	// next to leading order contribution      
	double dsigma_nlo = m_grids[1][iobs]->convolute( _pdf1, _pdf2, m_genpdf[1], alphas, m_leading_order+1, 0);
	// next to next to leading order contribution
	double dsigma_nnlo = m_grids[2][iobs]->convolute( _pdf1, _pdf2, m_genpdf[2], alphas, m_leading_order+2, 0);
	dsigma = dsigma_lo + dsigma_nlo + dsigma_nnlo;
      }
      else if ( nloops==-2 ) {
	label = "nnlo only";
	// next to next to leading order contribution
	dsigma = m_grids[2][iobs]->convolute( _pdf1, _pdf2, m_genpdf[2], alphas, m_leading_order+2, 0);
      }
      else { 
	throw grid::exception( std::cerr << "invalid value for nloops " << nloops ); 
      }

      double deltaobs = m_obs_bins->GetBinLowEdge(iobs+2)-m_obs_bins->GetBinLowEdge(iobs+1);      
      hvec.push_back( invNruns*Escale2*dsigma/deltaobs );
    }
    
  }
  else if ( m_type==AMCATNLO ) {  

    //    std::cout << "amc@NLO convolution" << std::endl;
    
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  

      double dsigma = 0; 
      
      if ( nloops==0 ) {
	  /// this is the amcatnlo LO calculation (without FKS shower)
	  label = "lo";
	  /// work out how to call from the igrid - maybe just implement additional 
	  double dsigma_B = m_grids[3][iobs]->amc_convolute( _pdf1, _pdf2, m_genpdf[3], alphas, m_leading_order,   0, 1, 1, Escale );
 
   	  dsigma = dsigma_B;
      }
      else if ( nloops==1 || nloops==-1 ) {
  	  /// this is the amcatnlo NLO calculation (without FKS shower)
	  label = "nlo only"; /// for the time being ...
	  /// work out how to call from the igrid - maybe just implement additional 
	  /// convolution routines and call them here
	  double dsigma_0 = m_grids[0][iobs]->amc_convolute( _pdf1, _pdf2, m_genpdf[0], alphas, m_leading_order+1, 0, 1, 1, Escale );
	  double dsigma_R = m_grids[1][iobs]->amc_convolute( _pdf1, _pdf2, m_genpdf[1], alphas, m_leading_order+1, 0, rscale_factor, 1, Escale );
	  double dsigma_F = m_grids[2][iobs]->amc_convolute( _pdf1, _pdf2, m_genpdf[2], alphas, m_leading_order+1, 0, 1, fscale_factor, Escale );
	  
	  dsigma = dsigma_0 + dsigma_R + dsigma_F;
      
	  if ( nloops==1 ) { 
	    /// this is the amcatnlo NLO calculation (without FKS shower)
	    label = "nlo";
	    /// work out how to call from the igrid - maybe just implement additional 
	    /// convolution routines and call them here
	    double dsigma_B = m_grids[3][iobs]->amc_convolute( _pdf1, _pdf2, m_genpdf[3], alphas, m_leading_order,   0, 1, 1, Escale );
	
	    dsigma += dsigma_B;
	  }

      }
      else { 
	throw grid::exception( std::cerr << "invalid value for nloops " << nloops ); 
      }

      double deltaobs = m_obs_bins->GetBinLowEdge(iobs+2)-m_obs_bins->GetBinLowEdge(iobs+1);      
      hvec.push_back( invNruns*Escale2*dsigma/deltaobs );
    }
  }
  else if ( m_type == SHERPA ) { 
    
    //    std::cout << "sherpa convolution" << std::endl;

    for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  
    
      double dsigma = 0; 
  
      if ( nloops==0 ) {
	label = "lo      ";
	// leading order cross section
	dsigma = m_grids[0][iobs]->convolute( _pdf1, _pdf2, m_genpdf[0], alphas, m_leading_order, 0, 1, 1, Escale);
      }
      else if ( nloops==1 ) { 
	label = "nlo     ";
	// next to leading order cross section
	// leading order contribution and scale dependent born dependent terms

	// will eventually add the other nlo terms ...
	double dsigma_lo  = m_grids[0][iobs]->convolute( _pdf1, _pdf2, m_genpdf[0], alphas, m_leading_order, 0, rscale_factor, fscale_factor, Escale);
	double dsigma_nlo = m_grids[1][iobs]->convolute( _pdf1, _pdf2, m_genpdf[1], alphas, m_leading_order+1, 0, rscale_factor, fscale_factor, Escale);
  
	dsigma = dsigma_lo + dsigma_nlo;
      }
      else if ( nloops==-1 ) { 
	label = "nlo     ";
	// next to leading order cross section
	// leading order contribution and scale dependent born dependent terms

	double dsigma_nlo = m_grids[1][iobs]->convolute( _pdf1, _pdf2, m_genpdf[1], alphas, m_leading_order+1, 0, rscale_factor, fscale_factor, Escale);

	dsigma = dsigma_nlo;
      }


      double deltaobs = m_obs_bins->GetBinLowEdge(iobs+2)-m_obs_bins->GetBinLowEdge(iobs+1);      
      hvec.push_back( invNruns*Escale2*dsigma/deltaobs );
    }

  }

    
  //  double _ctime = appl_timer_stop(_ctimer);
  //  std::cout << "grid::convolute() " << label << " convolution time=" << _ctime << " ms" << std::endl;
  
  if ( getApplyCorrections() ) applyCorrections(hvec);
  else { 
    for ( unsigned i=0 ; i<m_corrections.size() ; i++ ) if ( getApplyCorrection(i) ) applyCorrection(i,hvec);
  }

  cache1.stats();
  if ( cache2.ncalls() ) cache2.stats();
  
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


std::vector<double> appl::grid::vconvolute_subproc(int subproc,
						   void (*pdf)(const double& , const double&, double* ), 
						   double (*alphas)(const double& ), 
						   int     nloops, 
						   double  rscale_factor, double Escale )
{ 

  NodeCache cache1 = NodeCache( pdf );
  NodeCache* _pdf = &cache1;
  
  cache1.reset();

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

  std::string label;

  int lo_order = m_leading_order;
  if ( nloops>=m_order ) { 
    std::cerr << "too many loops for grid nloops=" << nloops << "\tgrid=" << m_order << std::endl;   
    return hvec;
  } 

  for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) {  

    double dsigma = 0;
   
    if ( nloops==0 ) {
      label = "lo      ";
      //      std::cout << "convolute() nloop=0" << iobs << std::endl;
      // leading order cross section
      dsigma = m_grids[0][iobs]->convolute_subproc(subproc, _pdf, 0, m_genpdf[0], alphas, lo_order, 0, 1, 1, Escale);
    }
    else if ( nloops==1 ) { 
      label = "nlo     ";
      // next to leading order cross section
      // leading and next to order contributions and scale dependent born dependent terms
      double dsigma_lo  = m_grids[0][iobs]->convolute_subproc(subproc, _pdf, 0, m_genpdf[0], alphas, lo_order,   1, rscale_factor, 1,  Escale );
      double dsigma_nlo = m_grids[1][iobs]->convolute_subproc(subproc, _pdf, 0, m_genpdf[1], alphas, lo_order+1, 0, rscale_factor, 1,  Escale );
      dsigma = dsigma_lo + dsigma_nlo;
    }
    else if ( nloops==-1 ) { 
      label = "nlo only";
      // nlo contribution only (only strict nlo contributions)
      dsigma = m_grids[1][iobs]->convolute_subproc(subproc, _pdf, 0, m_genpdf[1], alphas, lo_order+1, 0, rscale_factor, 1, Escale );
    }
    else if ( nloops==2 ) { 
      // FIXME: not implemented completely yet 
      return hvec;
      label = "nnlo    ";
      // next to next to leading order contribution 
      // NB: NO scale dependendent parts, so only muR=muF=mu
      double dsigma_lo   = m_grids[0][iobs]->convolute_subproc(subproc, _pdf, 0, m_genpdf[0], alphas, lo_order,   0);
      double dsigma_nlo  = m_grids[1][iobs]->convolute_subproc(subproc, _pdf, 0, m_genpdf[1], alphas, lo_order+1, 0);
      double dsigma_nnlo = m_grids[2][iobs]->convolute_subproc(subproc, _pdf, 0, m_genpdf[2], alphas, lo_order+2, 0);
      dsigma = dsigma_lo + dsigma_nlo + dsigma_nnlo;
    }
    
    
    //   double deltaobs = h->GetBinLowEdge(iobs+2)-h->GetBinLowEdge(iobs+1);
    //   h->SetBinContent(iobs+1, dsigma/(deltaobs));
    //   h->SetBinError(iobs+1, 0);

    double deltaobs = m_obs_bins->GetBinLowEdge(iobs+2)-m_obs_bins->GetBinLowEdge(iobs+1);
    
    hvec.push_back( invNruns*Escale2*dsigma/deltaobs );
    // hvec.push_back( Escale2*dsigma/deltaobs );


    //    std::cout << "dsigma[" << iobs << "]=" << dsigma/deltaobs << std::endl;

    
    //    std::cout << "obs bin " << iobs 
    //         << "\t" <<  h->GetBinLowEdge(iobs+1) << " - " << h->GetBinLowEdge(iobs+2)
    //	       << "\tdsigma=" << dsigma << std::endl;
  }  // iobs   

  //  double _ctime = appl_timer_stop(_ctimer);
  //  std::cout << "grid::convolute_subproc(" << subproc << ") " << label << " convolution time=" << _ctime << " ms" << std::endl;

  if ( getApplyCorrections() ) applyCorrections(hvec);
  else { 
    for ( unsigned i=0 ; i<m_corrections.size() ; i++ ) if ( getApplyCorrection(i) ) applyCorrection(i,hvec);
  }

  cache1.stats();

  return hvec;
}


TH1D* appl::grid::convolute(void (*pdf)(const double& , const double&, double* ), 
			    double (*alphas)(const double& ), 
			    int     nloops, 
			    double  rscale_factor,
			    double  fscale_factor,
			    double Escale )
{
  return convolute( pdf, 0, alphas, nloops, rscale_factor, fscale_factor, Escale );
}



TH1D* appl::grid::convolute(void (*pdf1)(const double& , const double&, double* ), 
			    void (*pdf2)(const double& , const double&, double* ), 
			    double (*alphas)(const double& ), 
			    int     nloops, 
			    double  rscale_factor,
			    double  fscale_factor,
			    double Escale ) {

    TH1D* h = new TH1D(*m_obs_bins);
    h->SetName("xsec");
    
    std::vector<double> dvec = vconvolute( pdf1, pdf2, alphas, nloops, rscale_factor, fscale_factor, Escale );
    
    for ( unsigned i=0 ; i<dvec.size() ; i++ ) { 
      h->SetBinContent( i+1, dvec[i] );
      h->SetBinError( i+1, 0 );
    }
    
    return h;

}





TH1D* appl::grid::convolute_subproc(int subproc,
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





void appl::grid::optimise(bool force) {
  if ( !force && m_optimised ) return;
  m_optimised = true;
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  { 
      std::cout << "grid::optimise() bin " << iobs << "\t";
      m_grids[iorder][iobs]->optimise();
    }
  }
  m_obs_bins->Reset();
}

void appl::grid::optimise(int NQ2, int Nx) {  optimise(NQ2, Nx, Nx);  }

void appl::grid::optimise(int NQ2, int Nx1, int Nx2) {
  m_optimised = true;
  for ( int iorder=0 ; iorder<m_order ; iorder++ ) { 
    for ( int iobs=0 ; iobs<Nobs() ; iobs++ )  { 
      std::cout << "grid::optimise() bin " << iobs << "\t";
      m_grids[iorder][iobs]->optimise(NQ2, Nx1, Nx2);
    }
  }
  m_obs_bins->Reset();
}




// redefine the limits by hand
void appl::grid::redefine(int iobs, int iorder,
			  int NQ2, double Q2min, double Q2max, 
			  int Nx,  double  xmin, double  xmax ) 
{ 
  
  if  ( iorder>=m_order ) { 
    std::cerr << "grid does not extend to this order" << std::endl;
    return;
  }
  
  if ( iobs<0 || iobs>=Nobs() ) { 
    std::cerr << "observable bin out of range" << std::endl;
    return;
  }
  
  if ( iorder==0 ) { 
    std::cout << "grid::redefine() iobs=" << iobs 
	      << "NQ2="  << NQ2 << "\tQmin=" << std::sqrt(Q2min) << "\tQmax=" << std::sqrt(Q2max) 
	      << "\tNx=" << Nx  << "\txmin=" <<            xmin  << "\txmax=" <<            xmax << std::endl; 
  }
  
  igrid* oldgrid =  m_grids[iorder][iobs];
  
  //  m_grids[iorder][iobs]->redefine(NQ2, Q2min, Q2max, Nx, xmin, xmax);

  m_grids[iorder][iobs] = new igrid(NQ2, Q2min, Q2max, oldgrid->tauorder(),
				    Nx,  xmin,  xmax,  oldgrid->yorder(), 
				    oldgrid->transform(), m_genpdf[iorder]->Nproc());

  delete oldgrid;
}
  


void appl::grid::setRange(int ilower, int iupper, double xScaleFactor) { 
  if ( ilower>=0 && iupper <Nobs() ) {  
    double lower = getReference()->GetBinLowEdge(ilower+1);
    double upper = getReference()->GetBinLowEdge(iupper+2); 
    setRange( lower, upper, xScaleFactor );
  }
}


void appl::grid::setRange(double lower, double upper, double xScaleFactor) { 
  
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

  igrid** grids[appl::MAXGRIDS];

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
void appl::grid::setDocumentation(const std::string& s) { m_documentation = s; }
void appl::grid::addDocumentation(const std::string& s) {   
  if ( m_documentation.size() ) m_documentation += s;
  else                          setDocumentation(s);    
}






/// methods to handle bin-by-bin corrections

/// add a correction as a std::vector
void appl::grid::addCorrection( std::vector<double>& v, const std::string& label) {
  //  std::cout << "addCorrections(vector) " << v.size() << " " << m_obs_bins->GetNbinsX() << std::endl;
  if ( v.size()==unsigned(m_obs_bins->GetNbinsX()) ) {
    m_corrections.push_back(v);
    m_correctionLabels.push_back(label);
    m_applyCorrection.push_back(false);
    //  std::cout << "appl::grid::addCorrection(vector) now " << m_corrections.size() << " corrections" << std::endl;
  }
}


/// add a correction by histogram
void appl::grid::addCorrection(TH1D* h, const std::string& label) {
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



// find the number of words used for storage
int appl::grid::size() const { 
    int _size = 0;
    for( int iorder=0 ; iorder<2 ; iorder++ ) {
      for( int iobs=0 ; iobs<Nobs() ; iobs++ ) _size += m_grids[iorder][iobs]->size();
    }
    return _size;
}


/// apply corrections to a std::vector
void appl::grid::applyCorrections(std::vector<double>& v) {
  //  std::cout << "grid::applyCorrections(vector) " << m_corrections.size() << std::endl;
  for ( unsigned i=0 ; i<m_corrections.size() ; i++ ) { 
    std::vector<double>& correction = m_corrections[i];
    //      TH1D* hc = m_corrections[i];
    for ( unsigned j=0 ; j<v.size() ; j++ ) v[j] *= correction[j];
  }
  //  std::cout << "grid::applyCorrections(vector) done" << std::endl;
}


/// apply correction to a std::vector
void appl::grid::applyCorrection(unsigned i, std::vector<double>& v) {
  //  std::cout << "grid::applyCorrections(vector) " << m_corrections.size() << std::endl;
  for ( unsigned j=0 ; j<m_corrections.size() ; j++ ) { 
    if ( j==i ) { 
      std::vector<double>& correction = m_corrections[j];
      //      TH1D* hc = m_corrections[i];
      for ( unsigned k=0 ; k<v.size() ; k++ ) v[k] *= correction[k];
    }
  }
  //  std::cout << "grid::applyCorrections(vector) done" << std::endl;
}




std::ostream& operator<<(std::ostream& s, const appl::grid& g) {
  s << "==================================================" << std::endl;
  //  s << "appl::grid version " << g.version() << "\t(" << g.subProcesses(0) << " initial states, " << g.Nobs() << " observable bins)" << std::endl;

  std::string basis[5] = {  "-LO, ",  "-NLO, ",  "-NNLO, ", "-Xtra0", "-Xtra1" };  
  std::string order[appl::MAXGRIDS];
  for ( int i=0 ; i<appl::MAXGRIDS ; i++ ) { 
    if ( i<5) order[i] = basis[i];
    else      order[i] = "-Unknown";
  }

  s << "appl::grid version " << g.version() << "\t( "; 
  for ( int i=0 ; i<g.nloops()+1 ; i++ ) s << g.subProcesses(i) << order[i];
  s << "initial states, " << g.Nobs() << " observable bins )" << std::endl;
  if ( g.isOptimised() ) s << "Optimised grid" << std::endl;
  if ( g.isSymmetric() ) s << "Symmetrised in x1, x2" << std::endl;
  else                   s << "Unsymmetrised in x1, x2" << std::endl;
  s << "leading order of processes  "  << g.leadingOrder() << std::endl;
  s << "number of loops for grid    " << g.nloops() << std::endl;   
  s << "x->y coordinate transform:  "  << g.getTransform() << std::endl;
  s << "genpdf in use: " << g.getGenpdf() << std::endl;
  s << "--------------------------------------------------" << std::endl;
  s << "Observable binning: [ " << g.Nobs() 
    << " bins : " << g.obsmin() << ",  " << g.obsmax() << " ]" << std::endl;

  //  for( int iorder=0 ; iorder<1 ; iorder++ ) {
  for( int iobs=0 ; iobs<g.Nobs() ; iobs++ ) {
    s << iobs << "\t" 
      << std::setprecision(5) << std::setw(5) << g.getReference()->GetBinLowEdge(iobs+1) << "\t- " 
      << std::setprecision(5) << std::setw(5) << g.getReference()->GetBinLowEdge(iobs+2) << "\t"; 
    s << "   " << *(g.weightgrid(0,iobs)) << std::endl;
  }
  //  }

  s << std::endl;
  
  return s;
}
