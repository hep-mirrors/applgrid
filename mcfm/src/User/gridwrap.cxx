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
  extern struct   {
    double Vsq[11][11],Vsum[11];
  } ckm_;
  
  
  extern struct {
    int nproc;
  } nproc_;
  
  extern struct  {
    int nflav;
  } nflav_;
  
  
  extern struct {
    bool creategrid;
    int nSubProcess;
  } grid_;
}



static const int mxpart = 12;        // mcfm parameter : max number of partons in event record

static const int _Ngrids = 3;
static const int  Ngrids = 3;
appl::grid* mygrid[_Ngrids];
static const char* gridFiles[_Ngrids] = {"weight_eta4.root",
					 "weight_pt4c.root",
					 "weight_pt4f.root"};

static double Observable[_Ngrids] = { 0.0, 0.0, 0.0 };   // observable array


int nObsBins[_Ngrids] = {20, 26, 26}; // eta4, pt4 cental eta-bin, pt4 forward eta-bin


static const double eta[] = { -4.0, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6,
			      -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,
			       1.6,  2.0,  2.4,  2.8,  3.2,  3.6,  4.0 };

static const double pt[] = {  0.,   5.,  10.,  15.,  20.,  25.,  30.,  35., 
			     40.,  45.,  50.,  55.,  60.,  65.,  70.,  75.,  
			     80.,  85.,  90.,  95., 100., 150., 200., 250., 
			    300., 400., 500.};



// static int count111 = 0;
// int nEvent = 100;
// int upCount = 0*3*(nEvent-1), dnCount = 3*(nEvent);

long unsigned int runs = 0;
TH1D ***procReference;


bool file_exists(const string& s) {   
  //  FILE* testfile=fopen(s.c_str(),"r");

  if ( FILE* testfile=fopen(s.c_str(),"r") ) { 
    fclose(testfile);
    return true;
  }
  else return false;
}

// #include "Readcards.h"


string glabel = "";


extern "C" void book_grid_()  // inital grid booking
{
  
  time_t _t;
  time(&_t);

  std::cout<<" ***********************************************"<<std::endl;
  std::cout<<" booking the grids " << ctime(&_t) << std::endl;

  // binning information for the grid constructor
  int    nXbins = 30;
  double xLow = 1.0e-5, xUp = 1.0;

  int    nQ2bins = 3;
  double q2Low = 6399.99, q2Up = 6400.01;

  // number of subrocesses for the process of interest
  // int nSubProcess = 6; // don't need this anymore - the pdf 
  // combination class knows how many


  // number of observables and binning for observables  
  const double *obsBins[_Ngrids] = { eta, pt, pt };


  bool write_ckm = false;

#if 0
  cout << "ckmsum: " << endl;
  for ( int i=0 ; i<11 ; i++ )  cout << "\t" << setw(15) << setprecision(12) << ckm_.Vsum[i];
  cout << "\n" << endl;
    

  cout << "ckmsum: " << endl;
  for ( int i=0 ; i<11 ; i++ ) { 
    for ( int j=0 ; j<11 ; j++ ) { 
      cout << "\t" << setw(15) << setprecision(12) << ckm_.Vsq[i][j];
    }
    cout << endl;
  }
#endif

  if ( write_ckm )   {
    cout << "writing ckm matrices" << endl;
    TFile ckmfile("ckm.root", "recreate");
    
    TMatrixT<double>* mat = new TMatrixT<double>(13,13);
    
    for (int f1 = 0; f1 <= 10; f1++) {
      for (int f2 = 0; f2 <= 10; f2++) (*mat)(f2+1,f1+1) = ckm_.Vsq[f1][f2];
    }
    
    TVectorT<double>* vec = new TVectorT<double>(13);
    for (int f1 = 0; f1 <= 10; f1++)  {
      (*vec)(f1+1) = ckm_.Vsum[f1];
    }
    
    mat->Write("ckm2");
    vec->Write("ckmsum");
    ckmfile.Close();
  }

  // NB don't know what the processes are - no documentation in
  // the mcfm code.

  string pdf_function;

  cout << "Process : " << nproc_.nproc;
  if      ( nproc_.nproc == 1 )  { cout << " W production"; pdf_function = "mcfm-w"; }  
  else if ( nproc_.nproc == 31 ) { cout << " Z production"; pdf_function = "mcfm-z"; }  
  else                           { cerr << "don't know which process" << endl; exit(-1); } 
  cout << endl;

  procReference = new TH1D**[Ngrids];


  //  set parameters from file
  //  ReadCards cards("cards.dat");

  string labels[3] = { "_eta", "_ptc", "_ptf" };

  //  glabel = cards.GetString("Label");
  //  nXbins = cards.GetValue("Nx");
  //  int xorder = cards.GetValue("xorder");

  glabel = "grid-30";
  nXbins = 30;
  int xorder = 5;


  if ( nproc_.nproc ==  1 ) glabel+="-W";
  if ( nproc_.nproc == 31 ) glabel+="-Z";


  for(int igrid=0; igrid < Ngrids; igrid++) {
      
      // if the file does not exist, create a new grid
      if ( ! file_exists(glabel+gridFiles[igrid]) ) { 
	cout<<"Creating NEW grid... "<<endl;
	
	// set transform2 value
	double apramval=5.;
	appl::igrid::transformvar(apramval);

	// lowest order in alphas	
	int lowest_order = 0;
	
	// how many loops
	int nloops = 1;
	
	double Nbins = 103;
	double bins[104];

	double del = ((*obsBins)[nObsBins[igrid]]-(*obsBins)[0])/nObsBins[igrid];

	for ( int ii=0 ; ii<=Nbins ; ii++ ) { 
	  bins[ii] = del*ii+(*obsBins)[0];
	}

	mygrid[igrid] = new appl::grid( nObsBins[igrid], obsBins[igrid], // obs bins
					nQ2bins, q2Low, q2Up, 0,         // Q2 bins and interpolation order
					nXbins,   xLow,  xUp, xorder,    // x bins and interpolation order
					pdf_function, lowest_order, nloops ); 

	//	mygrid[igrid]->symmetrise();
	
#if 0
	mygrid[igrid] = new appl::grid( Nbins, bins, // obs bins
					2, q2Low, q2Up, 1,         // Q2 bins and interpolation order
					12,   xLow,  xUp, 5,         // x bins and interpolation order
					pdf_function, lowest_order, nloops ); 
					//   "mcfm-w", lowest_order, nloops ); 
#endif       

	grid_.nSubProcess = mygrid[igrid]->subProcesses();
	cout << "reference histo name = " 
	     << mygrid[igrid]->getReference()->GetName() << endl;
	
	std::cout<<*mygrid[igrid]<<std::endl;  
      }
      else {
	//  cout<<" grid in " << gridFiles[igrid] << "... "<<endl;
	// user can optimise grid when they want, either before writing, of after 
	// reading, you decide!!
	//  mygrid[igrid] = new appl::grid(nXbins, nQ2bins, (char*)gridFiles[igrid]); //optimierte grid x,Q2 bins
	
	mygrid[igrid] = new appl::grid(glabel+gridFiles[igrid]); //optimise grid x,Q2 bins
	
	cout << "zero reference histogram" << endl;
	// zero reference histogram 
	TH1D* h = mygrid[igrid]->getReference();
	for ( int i=0 ; i<=h->GetNbinsX()+1 ; i++ ) h->SetBinContent(i, 0);
	
	mygrid[igrid]->optimise(nQ2bins, nXbins);
	
	grid_.nSubProcess = mygrid[igrid]->subProcesses();
	
	std::cout<<*(mygrid[igrid])<<std::endl;  
      }
      
      
      // what is this for ???
      double xbins[mygrid[igrid]->Nobs() + 1];
      for(int i = 1; i <= mygrid[igrid]->Nobs() + 1; i++) {
	xbins[i-1] = mygrid[igrid]->getReference()->GetBinLowEdge(i);
      }
      

      //  write out at the end ??? 
      //      gDirectory->cd(mygrid[igrid]->getFileName());
      procReference[igrid] = new TH1D*[mygrid[igrid]->subProcesses()];
      
      for (int isub = 0; isub < mygrid[igrid]->subProcesses(); isub++)
	{
	  char hname[50], htitle[50];
	  sprintf(hname,"procReference_%i",isub);
	  sprintf(htitle,"Contribution from SubProcess # %i",isub);
	  procReference[igrid][isub] = new TH1D(hname, htitle, mygrid[igrid]->Nobs(), xbins);
	  //	  procReference[igrid][isub]->SetDirectory();
	}

      //      gDirectory->cd("..");
  }
  
  runs = 0;
  
  std::cout<<" ***********************************************"<<std::endl;
}


extern "C" void book_grid__()  // inital grid booking
{
  book_grid_();
}



// just normalise to bin width
void Normalise(TH1D* h) { 
  for ( int ibin=1 ; ibin<=h->GetNbinsX() ; ibin++ ) { 
    double width = h->GetBinLowEdge(ibin+1) - h->GetBinLowEdge(ibin);
    
    h->SetBinContent( ibin, h->GetBinContent(ibin)/width );
    // if ( double e = h->GetBinError(ibin) ) h->SetBinError(ibin, e/width );
  }
}



extern "C" void write_grid_(double& xstotal)   // writes out grid after some events
{
  std::cout<<"\t Write out the grid ..."<<std::endl;
  for(int igrid = 0; igrid < Ngrids; igrid++)
    {
      cout << "\tsaving grid N=" << igrid+1 << "\tof total " << Ngrids << endl;
      
      mygrid[igrid]->run() = runs;
      // fix this up later
      //      mygrid[igrid]->setCrossSection(xstotal); 

      mygrid[igrid]->trim();
      int trim_size = mygrid[igrid]->size();

      //      mygrid[igrid]->print();

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

  
      TFile myfile((glabel+gridFiles[igrid]).c_str(),"update");
      {
	// Normalise the subprocesses referenece histograms 
	// by bin width also
	for (int isub = 0; isub < nsub; isub++)
	  {
	    Normalise(procReference[igrid][isub]);

	    procReference[igrid][isub]->Write("",TObject::kOverwrite);
	  }
	
      }
      myfile.Close();

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
  
  
  
  Observable[0] = rapidity4;
  Observable[1] = pt4;
  Observable[2] = pt4;
  
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
      (std::abs(Observable[0]) <= 0.5) ? fill = 1 : fill = 0;
//TC    (std::abs(Observable[0] <= 0.5)) ? fill = 1 : fill = 0;
      break;
    case(2):
//TC      (std::abs(Observable[0] >= 3.0)) ? fill = 1 : fill = 0;
      (std::abs(Observable[0]) >= 3.0) ? fill = 1 : fill = 0;
      break;
    default: 
      std::cerr<<" In gridwrap.cpp::cuts(int). No such process : "<<igrid<<std::endl;
      exit(-1);
    }
  return fill;
}






extern "C" void  fill_grid_(const double &x1, 
			    const double &x2, 
			    const double& fac, 
			    const double evt[][mxpart], const double *wt, 
			    const int &order)
{
  double scale2 = fac*fac;
  getObservable(evt);
  
  for(int igrid = 0; igrid < Ngrids; igrid++)
    {
      if (!(mygrid[igrid]->isOptimised()))
	{
	  //   cout << "not optimised";
	  if(cuts(igrid)) mygrid[igrid]->fill_phasespace(x1, x2, scale2, Observable[igrid], wt, order);
	}
      else
	{
	  //   cout << "optimised";
	  if(cuts(igrid)) mygrid[igrid]->fill(x1, x2, scale2, Observable[igrid], wt, order);
	}
    }
  runs++; // counter of number of events (shouldn't this be after cuts)? or is it the number of runs?"
}


extern "C" void  fill_grid__(const double &x1, 
			     const double &x2, 
			     const double& fac, 
			     const double evt[][mxpart], const double *wt, 
			     const int &order)
{
  fill_grid_( x1, x2, fac, evt, wt, order);
}


extern "C" void  fill_grid_reference_(const double evt[][mxpart], 
				      const double &wt, 
				      const double &wt2) // fills reference
{
  getObservable(evt);
  for(int igrid = 0; igrid < Ngrids; igrid++)
    {
      if(!((Observable[igrid] > mygrid[igrid]->obsmax()) || (Observable[igrid] < mygrid[igrid]->obsmin())) && cuts(igrid)) 
	{
	  //	  mygrid[igrid]->getReference()->Fill(Observable[igrid], wt/mygrid[igrid]->deltaobs(Observable[igrid]));
	  mygrid[igrid]->getReference()->Fill(Observable[igrid], wt);

	  // add a histo of squared of weights
	  //	  mygrid[igrid]->getReference2()->Fill(Observable[igrid],wt2/mygrid[igrid]->deltaobs(Observable[igrid]));
	  //  mygrid[igrid]->getReference2()->Fill(Observable[igrid],wt2);
	}
    }
}

extern "C" void  fill_grid_reference__(const double evt[][mxpart], 
				      const double &wt, 
				      const double &wt2) // fills reference
{ 
  fill_grid_reference_(evt, wt, wt2);
}


extern "C" void  fill_process_reference_(const double evt[][mxpart], 
					 const double *wt) // fills sub process reference
{
  getObservable(evt);
  for(int igrid = 0; igrid < Ngrids; igrid++)
    {
      if(!((Observable[igrid] > mygrid[igrid]->obsmax()) || (Observable[igrid] < mygrid[igrid]->obsmin())) && cuts(igrid)) 
	{
	  for (int isub = 0; isub < mygrid[igrid]->subProcesses(); isub++)
	    {
	      //      double currentWeight = wt[isub]/(mygrid[igrid]->deltaobs(Observable[igrid]));
	      // procReference[igrid][isub]->Fill(Observable[igrid], currentWeight);
	      procReference[igrid][isub]->Fill(Observable[igrid], wt[isub]);
	    }
	}
    }
}

extern "C" void  fill_process_reference__(const double evt[][mxpart], 
					  const double *wt) 
{
  fill_process_reference_( evt, wt );
}
