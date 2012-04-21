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

#include "mcfm_grid.h"


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



bool file_exists(const string& s) {   
  //  FILE* testfile=fopen(s.c_str(),"r");

  if ( FILE* testfile=fopen(s.c_str(),"r") ) { 
    fclose(testfile);
    return true;
  }
  else return false;
}


static const int mxpart = 12;        // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

static const int _Ngrids = 4;
static const int  Ngrids = _Ngrids;
appl::mcfm_grid* mygrid[_Ngrids];
static const char* gridFiles[_Ngrids] = 
  {
    "_eta3.root",
    "_pt3.root",
    "_eta4.root",
    "_pt4.root"
//    "_eta34.root",
//    "_pt34.root",
//    "_mass34.root"
  };

static double Observable[_Ngrids] = { 0.0 };                                // observable array
//int nObsBins[_Ngrids] = { 50, 100, 50, 100, 50, 100, 100 };              // eta4, pt4 cental eta-bin, pt4 forward eta-bin
//int nObsBins[_Ngrids] = { 40, 200, 40, 200 };                             // eta4, pt4 cental eta-bin, pt4 forward eta-bin
//int nObsBins[_Ngrids] = {100, 50, 100, 50};                               // eta4, pt4 cental eta-bin, pt4 forward eta-bin


int nObsBins[_Ngrids] = {100, 50, 100, 50}; // eta4, pt4 cental eta-bin, pt4 forward eta-bin

static const double eta[] = {
  -5.1, -4.9 , -4.8 , -4.7 , -4.6 , -4.5 , -4.4 , -4.3 , -4.2 , -4.1 ,
  -4 , -3.9 , -3.8 , -3.7 , -3.6 , -3.5 , -3.4 , -3.3 , -3.2 , -3.1 ,
  -3 , -2.9 , -2.8 , -2.7 , -2.6 , -2.5 , -2.4 , -2.3 , -2.2 , -2.1 ,
  -2 , -1.9 , -1.8 , -1.7 , -1.6 , -1.5 , -1.4 , -1.3 , -1.2 , -1.1 ,
  -1 , -0.9 , -0.8 , -0.7 , -0.6 , -0.5 , -0.4 , -0.3 , -0.2 , -0.1 ,
  0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 ,
  1 , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
  2 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , 2.6 , 2.7 , 2.8 , 2.9 ,
  3 , 3.1 , 3.2 , 3.3 , 3.4 , 3.5 , 3.6 , 3.7 , 3.8 , 3.9 ,
  4 , 4.1 , 4.2 , 4.3 , 4.4 , 4.5 , 4.6 , 4.7 , 4.8 , 4.9, 5.1
};
static const double pt[] = {
   0 , 2 , 4 , 6 , 8 , 10 , 12 , 14 , 16 , 18 , 20 , 22 , 24 , 26 , 28 ,
  30 , 32 , 34 , 36 , 38 , 40 , 42 , 44 , 46 , 48 , 50 , 52 , 54 , 56 , 58 ,
  60 , 62 , 64 , 66 , 68 , 70 , 72 , 74 , 76 , 78 , 80 , 82 , 84 , 86 , 88 ,
  90, 92, 94, 96, 98,100
};


//static const double eta[] = 
//  { 
//    -4 , -3.8 , -3.6 , -3.4 , -3.2 , -3 , -2.8 , -2.6 , -2.4 , -2.2 , 
//    -2 , -1.8 , -1.6 , -1.4 , -1.2 , -1 , -0.8 , -0.6 , -0.4 , -0.2 , 
//    0 , 0.2 , 0.4 , 0.6 , 0.8 , 1 , 1.2 , 1.4 , 1.6 , 1.8 , 2 , 
//    2.2 , 2.4 , 2.6 , 2.8 , 3 , 3.2 , 3.4 , 3.6 , 3.8 , 4 
//  };
// static const double eta[] = { 
// -5 , -4.8 , -4.6 , -4.4 , -4.2 , -4 , -3.8 , -3.6 , -3.4 , -3.2 , -3 , -2.8 , -2.6 , -2.4 , -2.2 , 
// -2 , -1.8 , -1.6 , -1.4 , -1.2 , -1 , -0.8 , -0.6 , -0.4 , -0.2 , 0 , 0.2 , 0.4 , 0.6 , 0.8 , 
// 1 , 1.2 , 1.4 , 1.6 , 1.8 , 2 , 2.2 , 2.4 , 2.6 , 2.8 , 3 , 3.2 , 3.4 , 3.6 , 3.8 , 4 , 4.2 , 4.4 , 4.6 , 4.8 , 5  
//   -5 , -4.9 , -4.8 , -4.7 , -4.6 , -4.5 , -4.4 , -4.3 , -4.2 , -4.1 , 
//   -4 , -3.9 , -3.8 , -3.7 , -3.6 , -3.5 , -3.4 , -3.3 , -3.2 , -3.1 , 
//   -3 , -2.9 , -2.8 , -2.7 , -2.6 , -2.5 , -2.4 , -2.3 , -2.2 , -2.1 , 
//   -2 , -1.9 , -1.8 , -1.7 , -1.6 , -1.5 , -1.4 , -1.3 , -1.2 , -1.1 , 
//   -1 , -0.9 , -0.8 , -0.7 , -0.6 , -0.5 , -0.4 , -0.3 , -0.2 , -0.1 , 
//   0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 
//   1 , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 , 
//   2 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , 2.6 , 2.7 , 2.8 , 2.9 , 
//   3 , 3.1 , 3.2 , 3.3 , 3.4 , 3.5 , 3.6 , 3.7 , 3.8 , 3.9 , 
//   4 , 4.1 , 4.2 , 4.3 , 4.4 , 4.5 , 4.6 , 4.7 , 4.8 , 4.9 , 5 
// };

//static const double pt[] = {
//  0 , 25 , 50 , 75 , 100 , 125 , 150 , 175 , 200 , 225 , 250 , 275 , 300 , 325 , 350 , 375 , 400 , 425 , 450 , 475 , 
//  500 , 525 , 550 , 575 , 600 , 625 , 650 , 675 , 700 , 725 , 750 , 775 , 800 , 825 , 850 , 875 , 900 , 925 , 950 , 975 , 
//  1000 , 1025 , 1050 , 1075 , 1100 , 1125 , 1150 , 1175 , 1200 , 1225 , 1250 , 1275 , 1300 , 1325 , 1350 , 1375 , 1400 , 1425 , 1450 , 1475 , 
//  1500 , 1525 , 1550 , 1575 , 1600 , 1625 , 1650 , 1675 , 1700 , 1725 , 1750 , 1775 , 1800 , 1825 , 1850 , 1875 , 1900 , 1925 , 1950 , 1975 , 
//  2000 , 2025 , 2050 , 2075 , 2100 , 2125 , 2150 , 2175 , 2200 , 2225 , 2250 , 2275 , 2300 , 2325 , 2350 , 2375 , 2400 , 2425 , 2450 , 2475 , 2500 
//   0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 , 80 , 85 , 90 , 95 , 100 , 
// 105 , 110 , 115 , 120 , 125 , 130 , 135 , 140 , 145 , 150 , 155 , 160 , 165 , 170 , 175 , 180 , 185 , 190 , 195 , 
// 200 , 205 , 210 , 215 , 220 , 225 , 230 , 235 , 240 , 245 , 250 , 255 , 260 , 265 , 270 , 275 , 280 , 285 , 290 , 295 , 
// 300 , 305 , 310 , 315 , 320 , 325 , 330 , 335 , 340 , 345 , 350 , 355 , 360 , 365 , 370 , 375 , 380 , 385 , 390 , 395 , 
// 400 , 405 , 410 , 415 , 420 , 425 , 430 , 435 , 440 , 445 , 450 , 455 , 460 , 465 , 470 , 475 , 480 , 485 , 490 , 495 , 
// 500 , 505 , 510 , 515 , 520 , 525 , 530 , 535 , 540 , 545 , 550 , 555 , 560 , 565 , 570 , 575 , 580 , 585 , 590 , 595 , 
// 600 , 605 , 610 , 615 , 620 , 625 , 630 , 635 , 640 , 645 , 650 , 655 , 660 , 665 , 670 , 675 , 680 , 685 , 690 , 695 , 
// 700 , 705 , 710 , 715 , 720 , 725 , 730 , 735 , 740 , 745 , 750 , 755 , 760 , 765 , 770 , 775 , 780 , 785 , 790 , 795 , 
// 800 , 805 , 810 , 815 , 820 , 825 , 830 , 835 , 840 , 845 , 850 , 855 , 860 , 865 , 870 , 875 , 880 , 885 , 890 , 895 , 
// 900 , 905 , 910 , 915 , 920 , 925 , 930 , 935 , 940 , 945 , 950 , 955 , 960 , 965 , 970 , 975 , 980 , 985 , 990 , 995 , 1000 
//};

//static const double mass[] = {
// 50 , 100 , 150 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 , 650 , 700 , 750 , 800 , 850 , 900 , 950 , 1000 , 
// 1050 , 1100 , 1150 , 1200 , 1250 , 1300 , 1350 , 1400 , 1450 , 1500 , 1550 , 1600 , 1650 , 1700 , 1750 , 1800 , 1850 , 1900 , 1950 , 
// 2000 , 2050 , 2100 , 2150 , 2200 , 2250 , 2300 , 2350 , 2400 , 2450 , 2500 , 2550 , 2600 , 2650 , 2700 , 2750 , 2800 , 2850 , 2900 , 2950 , 
// 3000 , 3050 , 3100 , 3150 , 3200 , 3250 , 3300 , 3350 , 3400 , 3450 , 3500 , 3550 , 3600 , 3650 , 3700 , 3750 , 3800 , 3850 , 3900 , 3950 , 
// 4000 , 4050 , 4100 , 4150 , 4200 , 4250 , 4300 , 4350 , 4400 , 4450 , 4500 , 4550 , 4600 , 4650 , 4700 , 4750 , 4800 , 4850 , 4900 , 4950 , 5000 , 5050 
//   0 , 15 , 30 , 45 , 60 , 75 , 90 , 105 , 120 , 135 , 150 , 165 , 180 , 195 , 210 , 225 , 240 , 255 , 270 , 285 ,
// 300 , 315 , 330 , 345 , 360 , 375 , 390 , 405 , 420 , 435 , 450 , 465 , 480 , 495 , 510 , 525 , 540 , 555 , 570 , 585 , 
// 600 , 615 , 630 , 645 , 660 , 675 , 690 , 705 , 720 , 735 , 750 , 765 , 780 , 795 , 810 , 825 , 840 , 855 , 870 , 885 , 
// 900 , 915 , 930 , 945 , 960 , 975 , 990 , 1005 , 1020 , 1035 , 1050 , 1065 , 1080 , 1095 , 1110 , 1125 , 1140 , 1155 , 1170 , 1185 , 
// 1200 , 1215 , 1230 , 1245 , 1260 , 1275 , 1290 , 1305 , 1320 , 1335 , 1350 , 1365 , 1380 , 1395 , 1410 , 1425 , 1440 , 1455 , 1470 , 1485 , 
// 1500 , 1515 , 1530 , 1545 , 1560 , 1575 , 1590 , 1605 , 1620 , 1635 , 1650 , 1665 , 1680 , 1695 , 1710 , 1725 , 1740 , 1755 , 1770 , 1785 , 
// 1800 , 1815 , 1830 , 1845 , 1860 , 1875 , 1890 , 1905 , 1920 , 1935 , 1950 , 1965 , 1980 , 1995 , 2010 , 2025 , 2040 , 2055 , 2070 , 2085 , 
// 2100 , 2115 , 2130 , 2145 , 2160 , 2175 , 2190 , 2205 , 2220 , 2235 , 2250 , 2265 , 2280 , 2295 , 2310 , 2325 , 2340 , 2355 , 2370 , 2385 , 
// 2400 , 2415 , 2430 , 2445 , 2460 , 2475 , 2490 , 2505 , 2520 , 2535 , 2550 , 2565 , 2580 , 2595 , 2610 , 2625 , 2640 , 2655 , 2670 , 2685 , 
// 2700 , 2715 , 2730 , 2745 , 2760 , 2775 , 2790 , 2805 , 2820 , 2835 , 2850 , 2865 , 2880 , 2895 , 2910 , 2925 , 2940 , 2955 , 2970 , 2985 , 3000
// };
//  { 
//    0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 
//    75 , 80 , 85 , 90 , 95 , 100 , 105 , 110 , 115 , 120 , 125 , 130 , 135 , 140 , 145 , 150
//  };
// static const double pt[] = { 
//   0 , 2 , 4 , 6 , 8 , 10 , 12 , 14 , 16 , 18 , 20 , 22 , 24 , 26 , 28 , 
//   30 , 32 , 34 , 36 , 38 , 40 , 42 , 44 , 46 , 48 , 50 , 52 , 54 , 56 , 58 , 
//   60 , 62 , 64 , 66 , 68 , 70 , 72 , 74 , 76 , 78 , 80 , 82 , 84 , 86 , 88 , 
//   90 , 92 , 94 , 96 , 98 , 100
// };

long unsigned int runs  =  0;
bool isBooked           =  false;
string glabel           =  "";
string outputDirectory  =  "./output/";

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
//  const double *obsBins[_Ngrids] = { eta, pt, eta, pt, eta, pt, mass };
  const double *obsBins[_Ngrids] = { eta, pt, eta, pt };
  
  
  string pdf_function;
  char binsInfo[111];
//  string labels[ _Ngrids ] = { "_eta3", "_pt3", "_eta4", "_pt4", "_eta34", "_pt34", "_mass34" };
  string labels[ _Ngrids ] = { "_eta3", "_pt3", "_eta4", "_pt4" };
  glabel = "grid-";
  
  cout << "Process : " << nproc_.nproc;
  if      ( nproc_.nproc == 1 )  
    {
      cout << " W+ production"; 
      pdf_function = "mcfm-wp";
      sprintf(binsInfo, "%d", nXbins);
      glabel+= binsInfo;
      glabel+="-Wplus";
    }  
  else if ( nproc_.nproc == 6 )  
    {
      cout << " W- production"; 
      pdf_function = "mcfm-wm"; 
      sprintf(binsInfo, "%d", nXbins);
      glabel+= binsInfo;
      glabel+="-Wminus";
    }  
  else if ( nproc_.nproc == 31 ) 
    {
      cout << " Z production"; 
      pdf_function = "mcfm-z"; 
      sprintf(binsInfo, "%d", nXbins);
      glabel+= binsInfo;
      glabel+="-Z0";
      q2Low = 8280.99, q2Up = 8281.01;
    }  
  else if ( (nproc_.nproc == 157) || (nproc_.nproc == 158) || (nproc_.nproc==159) )
    {

      nXbins  = 45;
      xorder=5;
      nQ2bins = 15;
      qorder=3;
      lowest_order = 2;

      sprintf(binsInfo, "%d", nXbins);
      glabel+= binsInfo;

      if (nproc_.nproc == 157)
	{
	  cout << " TTbar production"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar";
          q2Low = std::pow(1.,2), q2Up = std::pow(3500.,2);
	}
      else if(nproc_.nproc == 158)
	{
	  cout << " BBbar production"; 
	  pdf_function = "mcfm-BB"; 
	  glabel+="-BBbar";
          q2Low = std::pow(1.,2), q2Up = std::pow(3500.,2);
	}
      else if (nproc_.nproc == 159)
	{
	  cout << " CCbar production"; 
	  pdf_function = "mcfm-CC"; 
	  glabel+="-CCbar";
          q2Low = std::pow(1.,2), q2Up = std::pow(3500.,2);
	}

    }  
  else                           
    { 
      cerr << "don't know which process" << endl; 
      exit(-1); 
    } 
  cout << endl;
  
  for(int igrid=0; igrid < Ngrids; igrid++) 
    {
    
      bool create_new = false;

      // if the file does not exist, create a new grid...
      if ( !file_exists( outputDirectory + glabel+gridFiles[igrid]) )  create_new = true;

      // or if it does exists but root file is a zombie...
      if ( !create_new ) {  
	TFile testFile( ( outputDirectory  + glabel+gridFiles[igrid]).c_str() );
	if ( testFile.IsZombie() ) create_new = true;
	testFile.Close();
      }

      if ( create_new ) 
	{ 
	  cout<<"Creating NEW grid... "<<endl;
	  
	  mygrid[igrid] = new appl::mcfm_grid( nObsBins[igrid], obsBins[igrid],      // obs bins
					       nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
					       nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
					       pdf_function, lowest_order, nloops ); 
	  
	  mygrid[igrid]->setCMSScale( energy_.sqrts );
	  grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  
	//  cout << "reference histo name = " 
	//       << mygrid[igrid]->getReference()->GetName() << endl;
	  
	//  std::cout<<*mygrid[igrid]<<std::endl;  
	}
      else 
	{
	  std::cout << "Using existing grid file " << (outputDirectory + glabel+gridFiles[igrid]) << std::endl;
	  
	  mygrid[igrid] = new appl::mcfm_grid(outputDirectory + glabel+gridFiles[igrid]); //optimise grid x,Q2 bins
	  grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  mygrid[igrid]->getReference()->Reset();
	  mygrid[igrid]->optimise(nQ2bins, nXbins);
	  
	 // std::cout<<*(mygrid[igrid])<<std::endl;  
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
      
      mygrid[igrid]->Write(outputDirectory + glabel+gridFiles[igrid]);
      
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

  double mass34 = 0;
  mass34 = std::sqrt( std::pow(p3[0] + p4[0],2) - std::pow(p3[1] + p4[1],2) - std::pow(p3[2] + p4[2],2) - std::pow(p3[3] + p4[3],2) );
  
  Observable[ 0 ] = rapidity3;
  Observable[ 1 ] = pt3;
  Observable[ 2 ] = rapidity4;
  Observable[ 3 ] = pt4;
//  Observable[ 4 ] = rapidity34;
//  Observable[ 5 ] = pt34;
//  Observable[ 6 ] = mass34;
  
}

int cuts(int igrid)
{
  int fill;
//  if (Observable[1]<25. || Observable[3]<25. ) 
//    {
//	fill= 0 ;
//        return fill;
//    }	

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
    case(4):
      fill = 1;
      break;
    case(5):
      fill = 1;
      break;
    case(6):
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

