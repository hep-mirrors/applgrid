//
//   gridwrap.cxx        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: gridwrap.cxx, v   Sat May  3 10:15:04 BST 2008 sutt


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <iomanip>
using std::setprecision;
using std::setw;

#include <string>
using std::string;

#include <stdlib.h> // exit()
#include <sys/time.h> 


#include "TFile.h"
#include "TH1D.h"
#include "TMatrixT.h"
#include "TVectorT.h"

#include "appl_grid/appl_grid.h"
using appl::grid;
using appl::igrid;

// #include "appl_grid/appl_pdf.h"

extern "C" { 
  
  extern struct {
    int nproc;
  } nproc_;
  
  extern struct  {
    double sqrts;
  } energy_;
  
  
  extern struct {
    bool creategrid;
    int nSubProcess;
  } grid_;

}


static const int mxpart = 12;        // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

static const int _Ngrids = 4;
static const int  Ngrids = _Ngrids;
appl::grid* mygrid[_Ngrids];
static const char* gridFiles[_Ngrids] = 
  {
    "_eta3.root",
    "_pt3.root",
    "_eta4.root",
    "_pt4.root"
  };

static double Observable[_Ngrids] = { 0.0, 0.0, 0.0, 0.0 };   // observable array
int nObsBins[_Ngrids] = {100, 50, 100, 50}; // eta4, pt4 cental eta-bin, pt4 forward eta-bin

static const double eta[] = { 
  -5 , -4.9 , -4.8 , -4.7 , -4.6 , -4.5 , -4.4 , -4.3 , -4.2 , -4.1 , 
  -4 , -3.9 , -3.8 , -3.7 , -3.6 , -3.5 , -3.4 , -3.3 , -3.2 , -3.1 , 
  -3 , -2.9 , -2.8 , -2.7 , -2.6 , -2.5 , -2.4 , -2.3 , -2.2 , -2.1 , 
  -2 , -1.9 , -1.8 , -1.7 , -1.6 , -1.5 , -1.4 , -1.3 , -1.2 , -1.1 , 
  -1 , -0.9 , -0.8 , -0.7 , -0.6 , -0.5 , -0.4 , -0.3 , -0.2 , -0.1 , 
  0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 
  1 , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 , 
  2 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , 2.6 , 2.7 , 2.8 , 2.9 , 
  3 , 3.1 , 3.2 , 3.3 , 3.4 , 3.5 , 3.6 , 3.7 , 3.8 , 3.9 , 
  4 , 4.1 , 4.2 , 4.3 , 4.4 , 4.5 , 4.6 , 4.7 , 4.8 , 4.9 , 5 
};
static const double pt[] = { 
  0 , 2 , 4 , 6 , 8 , 10 , 12 , 14 , 16 , 18 , 20 , 22 , 24 , 26 , 28 , 
  30 , 32 , 34 , 36 , 38 , 40 , 42 , 44 , 46 , 48 , 50 , 52 , 54 , 56 , 58 , 
  60 , 62 , 64 , 66 , 68 , 70 , 72 , 74 , 76 , 78 , 80 , 82 , 84 , 86 , 88 , 
  90 , 92 , 94 , 96 , 98 , 100
};

long unsigned int runs  =  0;
bool isBooked           =  false;
string glabel           =  "";

void getObservable( const double evt[][mxpart] );
int  cuts(int);



extern "C" void book_grid_()  // inital grid booking
{
  if (isBooked) return;
  
  time_t _t;
  time(&_t);
  
  std::cout<<" ***********************************************"<<std::endl;
  std::cout<<" booking the grids " << ctime(&_t) << std::endl;
  
  // binning information for the grid constructor
  int    nXbins  = 30;
  double xLow    = 1.0e-8, xUp = 1.0;
  int    xorder  = 5;
  int    nQ2bins = 3;
  double q2Low   = 6399.99, q2Up = 6400.01;
  int    qorder  = 0;
  // set transform2 value
  double apramval=5.;
  appl::igrid::transformvar(apramval);
  
  // lowest order in alphas	
  int lowest_order = 0;
  // how many loops
  int nloops = 1;
  
  // number of observables and binning for observables  
  const double *obsBins[_Ngrids] = { eta, pt, eta, pt };
  
  
  string pdf_function;
  string labels[ _Ngrids ] = { "_eta3", "_pt3", "_eta4", "_pt4" };
  glabel = "grid-30";
  
  cout << "Process : " << nproc_.nproc;
  if      ( nproc_.nproc == 1 )  
    {
      cout << " W+ production"; 
      pdf_function = "mcfm-wp"; 
      glabel+="-Wplus";
    }  
  else if ( nproc_.nproc == 6 )  
    {
      cout << " W- production"; 
      pdf_function = "mcfm-wm"; 
      glabel+="-Wminus";
    }  
  else if ( nproc_.nproc == 31 ) 
    {
      cout << " Z production"; 
      pdf_function = "mcfm-z"; 
      glabel+="-Z0";
      q2Low = 8280.99, q2Up = 8281.01;
    }  
  else                           
    { 
      cerr << "don't know which process" << endl; 
      exit(-1); 
    } 
  cout << endl;
  
  for(int igrid=0; igrid < Ngrids; igrid++) 
    {
      
      // if the file does not exist, create a new grid
      TFile testFile( (glabel+gridFiles[igrid]).c_str() );
      if ( testFile.IsZombie() ) 
	{ 
	  cout<<"Creating NEW grid... "<<endl;
	  
	  mygrid[igrid] = new appl::grid( nObsBins[igrid], obsBins[igrid],      // obs bins
					  nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
					  nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
					  pdf_function, lowest_order, nloops ); 
	  
	  mygrid[igrid]->setCMSScale( energy_.sqrts );
	  grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  
	  cout << "reference histo name = " 
	       << mygrid[igrid]->getReference()->GetName() << endl;
	  
	  std::cout<<*mygrid[igrid]<<std::endl;  
	}
      else 
	{
	  
	  mygrid[igrid] = new appl::grid(glabel+gridFiles[igrid]); //optimise grid x,Q2 bins
	  grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  mygrid[igrid]->getReference()->Reset();
	  mygrid[igrid]->optimise(nQ2bins, nXbins);
	  
	  std::cout<<*(mygrid[igrid])<<std::endl;  
	}

    }

  runs = 0;
  isBooked = true;
  std::cout<<" ***********************************************"<<std::endl;
}


extern "C" void book_grid__()  // inital grid booking
{
  book_grid_();
}

extern "C" void  fill_grid_( const double evt[][mxpart] )
{
  if (!isBooked) 
    {    
      book_grid_();
      return;
    }

  getObservable( evt );
  
  for(int igrid = 0; igrid < Ngrids; igrid++)
    if(cuts(igrid))
      mygrid[igrid]->fillMCFM( Observable[igrid] );
  
  
  runs++; // counter of number of events (shouldn't this be after cuts)? or is it the number of runs?"
}

extern "C" void  fill_grid__( const double evt[][mxpart] )
{
  fill_grid_( evt );
}

//
// just normalise to bin width
//
void Normalise(TH1D* h) 
{ 
  for ( int ibin=1 ; ibin<=h->GetNbinsX() ; ibin++ ) 
    { 
      double width = h->GetBinLowEdge(ibin+1) - h->GetBinLowEdge(ibin);
      h->SetBinContent( ibin, h->GetBinContent(ibin)/width );
    }
  return;
}

extern "C" void write_grid_(double& xstotal)   // writes out grid after some events
{
  std::cout<<"\t Write out the grid ..."<<std::endl;
  for(int igrid = 0; igrid < Ngrids; igrid++)
    {
      cout << "\tsaving grid N=" << igrid+1 << "\tof total " << Ngrids << endl;
      
      mygrid[igrid]->setNormalised( true );
      mygrid[igrid]->run() = runs;
      
      mygrid[igrid]->trim();
      int trim_size = mygrid[igrid]->size();

      mygrid[igrid]->untrim();
      int untrim_size = mygrid[igrid]->size();
      
      // normalise the reference histogram by bin width
      Normalise( mygrid[igrid]->getReference() );
      
      mygrid[igrid]->Write(glabel+gridFiles[igrid]);
      
      cout << "size(untrimmed)=" << untrim_size 
	   << "\tsize(trimmed)=" << trim_size 
	   << "\tfraction="      << 100.*trim_size/untrim_size << " %" << endl;

      int nsub = mygrid[igrid]->subProcesses();

      delete mygrid[igrid];
      
    }
  
  time_t _t;
  time(&_t);
  
  std::cout<<" ***********************************************"<<std::endl;
  std::cout<<" saved grids " << ctime(&_t) << std::endl;
  std::cout<<" ***************************************************"<<std::endl;
}
 

extern "C" void write_grid__(double& xstotal)   // writes out grid after some events
{
  write_grid_(xstotal);   // writes out grid after some events
}

//
// ----------------------------------------------
//    analysis
// ----------------------------------------------
//

void getObservable(const double evt[][mxpart])
{
  // evt[momentum][particle number-1]
  // momentum[0,1,2,3] = (x,y,z,E)
  //

  // calculate observables
  for(int igrid = 0; igrid < Ngrids; igrid++)Observable[igrid] = 0.0; // initialize
  
  double p3[4] = {evt[3][2],evt[0][2],evt[1][2],evt[2][2]}; // (E,x,y,z)
  double p4[4] = {evt[3][3],evt[0][3],evt[1][3],evt[2][3]};
  
  
  double rapidity3 = 0.0;
  rapidity3 = (p3[0] + p3[3])/(p3[0] - p3[3]);
  (rapidity3 < 1e-13) ? rapidity3 = 100.0 : rapidity3 = 0.5*std::log(rapidity3);
  
  double rapidity4 = 0.0;
  rapidity4 = (p4[0] + p4[3])/(p4[0] - p4[3]);
  (rapidity4 < 1e-13) ? rapidity4 = 100.0 : rapidity4 = 0.5*std::log(rapidity4);
  
  double rapidity34 = 0.0;                      // rapidity of particle (3+4) in event record
  rapidity34  = (p3[0] + p4[0]) + (p3[3] + p4[3]);
  rapidity34 /= (p3[0] + p4[0]) - (p3[3] + p4[3]);
  
  (rapidity34 < 1e-13) ? rapidity34 = 100.0 : rapidity34 = 0.5*std::log(rapidity34);
  
  double pt3 = 0;
  pt3 = std::sqrt( p3[1]*p3[1] + p3[2]*p3[2] );
  
  double pt4 = 0;
  pt4 = std::sqrt( p4[1]*p4[1] + p4[2]*p4[2] );
  
  double pt34 = 0;
  pt34 = std::sqrt( std::pow(p3[1] + p4[1],2) + std::pow(p3[2] + p4[2],2) );
  
  Observable[ 0 ] = rapidity3;
  Observable[ 1 ] = pt3;
  Observable[ 2 ] = rapidity4;
  Observable[ 3 ] = pt4;
  
}

int cuts(int igrid)
{
  int fill;
  switch(igrid)
    {
    case(0):
      fill = 1;
      break;
    case(1):
      fill = 1;
      //      (std::abs(Observable[0]) <= 0.5) ? fill = 1 : fill = 0;
      //TC    (std::abs(Observable[0] <= 0.5)) ? fill = 1 : fill = 0;
      break;
    case(2):
      fill = 1;
      //TC      (std::abs(Observable[0] >= 3.0)) ? fill = 1 : fill = 0;
      //      (std::abs(Observable[0]) >= 3.0) ? fill = 1 : fill = 0;
      break;
    case(3):
      fill = 1;
      break;
    default: 
      std::cerr<<" In gridwrap.cpp::cuts(int). No such process : "<<igrid<<std::endl;
      exit(-1);
    }
  return fill;
}


// extern "C" void  fill_grid_reference_(const double evt[][mxpart], 
// 				      const double &wt, 
// 				      const double &wt2) // fills reference
// {
//   getObservable(evt);
//   for(int igrid = 0; igrid < Ngrids; igrid++)
//     if ( mygrid[igrid]->isOptimised() )
//       if( cuts(igrid) ) 
// 	if ( !( (Observable[igrid] > mygrid[igrid]->obsmax()) || (Observable[igrid] < mygrid[igrid]->obsmin())) )
// 	  mygrid[igrid]->getReference()->Fill(Observable[igrid], wt);
  
//   return;
// }

// extern "C" void  fill_grid_reference__(const double evt[][mxpart], 
// 				       const double &wt, 
// 				       const double &wt2) // fills reference
// { 
//   fill_grid_reference_(evt, wt, wt2);
// }

