



#include <iostream>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::endl;


//#include <gtttimer.h>
#include "appl_grid/appl_timer.h"
#include "lcg.h"

// io
#include "TFile.h"

// the matrix classes for testing 
#include "appl_grid/SparseMatrix3d.h"
#include "appl_grid/Sparse3d.h"
#include "TH3D.h"
#include "TMatrixDSparse.h"


#include <cmath>
using std::abs;
using std::sqrt;



class SP3d { 

public:

  SP3d(int nx, double lx, double ux, 
       int ny, double ly, double uy, 
       int nz, double lz, double uz) :  m_Nx(nx), mv(NULL) { 

    mv = new TMatrixDSparse*[m_Nx];
    m_lx = lx;
    m_ux = ux;

    for ( int i=0 ; i<m_Nx ; i++ ) {
      mv[i] = new TMatrixDSparse(ny, nz);
    }
  } 

  ~SP3d() {
    for ( int i=0 ; i<m_Nx ; i++ ) delete mv[i];
    delete[] mv;
  }


  double& operator()(int i, int j, int k) { 
    return (*mv[i])(j,k); 
  }


private:
  
  int m_Nx;
  int m_lx;
  int m_ux;
  
  TMatrixDSparse** mv;

};


double e2(double t0, double t1) { 
  return sqrt((t0*0.04)*(t0*0.04)+((t0+t1)*0.04)*((t0+t1)*0.04));
}







struct _time { 

  _time() : N(0), t(0), t2(0) { } 

  double t;
  double t2;
  
  double N;

  double m;
  double s;

  void add(double time) { 
    N++;
    t  += time;
    t2 += time*time;
    m = t/N;
    s = sqrt((t2/N-m*m)/N);
  }

  double mean() { 
    return m;
  }

  double sigma() {
    return s;
  }

  _time& operator=(double d) { 
    t  = d;
    t2 = d*d;
    N  = 1; 
    m  = d;
    s  = 0; 
  }

  _time& operator+=(double d) { 
    
  }

  //  _time& operator-(time_& time) { 
  //    m = (m - time.m);
  //  }

};





// void times(int N, double& t0, double& t1, double& t2, double& t3);
void times(int N, _time& t0, _time& t1, _time& t2, _time& t3);


#include "TPad.h"
#include "TStyle.h"

double quad(double e1, double e2) { 
  return sqrt(e1*e1+e2*e2);
}


int main() { 

  cout << "this runes a series of grid filling excercises" << endl;
  cout << "with different types of grids, and times them"  << endl;
  cout << " - there is the standard TH3D type of grid" << endl;
  cout << " - a grid based on a vector of TSparseMatrix2D" << endl;
  cout << " - and the custom SparseMatrix3d grid class" << endl; 

  TFile f("times.root", "recreate");

  int N=6;

  TH1D* hh0(new TH1D("c0","c0", 30, 0.0018, 0.0024));
  TH1D* h0(new TH1D("t0","t0",  N, 0, N*10)); h0->SetTitle("no filling");
  TH1D* h1(new TH1D("t1","t1",  N, 0, N*10)); h1->SetTitle("TH3D fill");
  TH1D* h2(new TH1D("t2","t2",  N, 0, N*10)); h2->SetTitle("SparseMatrix3D fill");
  TH1D* h3(new TH1D("t3","t3",  N, 0, N*10)); h3->SetTitle("Sparse3D (TSparse) fill");
	
  gStyle->SetErrorX(0.0001);
  
  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(24);
  h3->SetMarkerStyle(25);


  h0->SetMinimum(0);
  h1->SetMinimum(0);
  h2->SetMinimum(0);
  h3->SetMinimum(0);
  
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h1->SetMarkerColor(kRed);
  h2->SetMarkerColor(kBlue);

  //  double t0, t1, t2, t3;
  _time t0, t1, t2, t3;
  for ( int i=0 ; i<N ; i++ ) { 
    int size = h0->GetBinCenter(i+1);
    times(size, t0, t1, t2, t3);
    hh0->Fill(t0.mean());
    h0->SetBinContent(i+1,t0.mean());
#if 0
    h1->SetBinContent(i+1,t1.mean()-t0.mean());
    h2->SetBinContent(i+1,t2.mean()-t0.mean());
    h3->SetBinContent(i+1,t3.mean()-t0.mean());
#else
    h1->SetBinContent(i+1,t1.mean());
    h2->SetBinContent(i+1,t2.mean());
    h3->SetBinContent(i+1,t3.mean());
#endif
    //    h0->SetBinError(i+1,0.04*t0);
    //    h1->SetBinError(i+1,e2(t0,t1));
    //    h2->SetBinError(i+1,e2(t0,t2));
    //    h3->SetBinError(i+1,e2(t0,t3));
    h0->SetBinError(i+1,t0.sigma());
#if 0
    h1->SetBinError(i+1,quad(t1.sigma(),t0.sigma()));
    h2->SetBinError(i+1,quad(t2.sigma(),t0.sigma()));
    h3->SetBinError(i+1,quad(t3.sigma(),t0.sigma()));
#else
    h1->SetBinError(i+1,t1.sigma());
    h2->SetBinError(i+1,t2.sigma());
    h3->SetBinError(i+1,t3.sigma());
#endif

    h2->SetBinError(i+1,t2.sigma());
    h3->SetBinError(i+1,t3.sigma());

    h3->SetMinimum(0);
    h3->SetMaximum(0.01);
 
    h3->DrawCopy("e1");
    h3->DrawCopy("samelhist");
    h0->DrawCopy("samee1");
    h0->DrawCopy("samelhist");
    h1->DrawCopy("samee1");
    h1->DrawCopy("samelhist");
    h2->DrawCopy("samee1");
    h2->DrawCopy("samelhist");

    gPad->Print("plot.eps");
 }

  f.Write();
  f.Close();

  return 0;
}


// void times(int N, double& t0, double& t1, double& t2, double& t3) { 
void times(int N, _time& t0, _time& t1, _time& t2, _time& t3) { 

  TH3D*          h1(new TH3D("h3d", "h3d", N, 0, N, N, 0, N, N, 0, N));
  SparseMatrix3d h2(N, 0, N, N, 0, N, N, 0, N);
  //  SP3d           h3(N, 0, N, N, 0, N, N, 0, N);
  Sparse3d           h3(N, 0, N, N, 0, N, N, 0, N);
  // TMatrixDSparse     h3(20, 20);
  
  t0=t1=t2=t3=0;

  double TH3time = 0;
  double SP3time = 0;
  double SM3time = 0;

  double fac=1;

  int Nev = 10000;
  int Nloop = 100;


  double t0sum = 0;
  for ( int j=0 ; j<Nloop ; j++ ) { 
    struct timeval duff = appl::appl_timer_start();
    for ( int i=0 ; i<Nev ; i++ ) { 
      
      int x = 0.2*N*gauss()+0.5*N;
      int y = 0.2*N*gauss()+0.5*N;
      int z = 0.2*N*gauss()+0.5*N;
      
      if ( x<0 || x>=N ) continue;
      if ( y<0 || y>=N ) continue;
      if ( z<0 || z>=N ) continue;
      
      double w=1;
      
    }
    double dufftime = appl::appl_timer_stop(duff);
    dufftime /= Nev;
    t0.add(dufftime);
  }
  //  t0=dufftime;
  //  t0=t0sum/Nloop;
  
  //  cout << "th3d" << endl;
  double t1sum = 0;
  for ( int j=0 ; j<Nloop ; j++ ) { 
    struct timeval TH3timer = appl::appl_timer_start();
    for ( int i=0 ; i<Nev ; i++ ) { 
      
      int x = 0.2*N*gauss()+0.5*N;
      int y = 0.2*N*gauss()+0.5*N;
      int z = 0.2*N*gauss()+0.5*N;
      
      if ( x<0 || x>=N ) continue;
      if ( y<0 || y>=N ) continue;
      if ( z<0 || z>=N ) continue;
      
      double w=1;
      h1->Fill(x+1,y+1,z+1,w);
    }
    TH3time = appl::appl_timer_stop(TH3timer);
    // t1sum += TH3time/Nev;
    t1.add(TH3time/Nev);

  }
  t1sum /= Nloop;
  delete h1;

  // t1 = t1sum - t0;

  //  cout << "sm3d" << endl;
   double t2sum = 0;
  for ( int j=0 ; j<Nloop ; j++ ) { 
    struct timeval SM3timer = appl::appl_timer_start();
    for ( int i=0 ; i<Nev ; i++ ) { 
      
      int x = 0.02*N*gauss()+0.5*N;
      int y = 0.02*N*gauss()+0.5*N;
      int z = 0.02*N*gauss()+0.5*N;
      
      if ( x<0 || x>=N ) continue;
      if ( y<0 || y>=N ) continue;
      if ( z<0 || z>=N ) continue;
      
      double w=1;
      h2.fill_fast(x,y,z) += w;
      // h2(x,y,z) += w;
    }    
    SM3time  = appl::appl_timer_stop(SM3timer);
    // t2sum += SM3time/Nev;
    t2.add(SM3time/Nev);
  }
  //  t2sum /= Nloop;
  //  t2 =t2sum - t0;
  
  
  h2.trim();
  
  //  h2.print();
  
  
  if ( N>60 ) fac = 4; 
  if ( N>80 ) fac = 5; 
  //  if ( N>90 ) fac = 8; 
  
  Nev /= fac;
  //  if ( N<90 )
  // cout << "sp3d" << endl;
  double t3sum = 0;
  for ( int j=0 ; j<Nloop ; j++ ) { 
    struct timeval SP3timer = appl::appl_timer_start();
    for ( int i=0 ; i<Nev ; i++ ) { 
      
      int x = 0.2*N*gauss()+0.5*N;
      int y = 0.2*N*gauss()+0.5*N;
      int z = 0.2*N*gauss()+0.5*N;
      
      if ( x<0 || x>=N ) continue;
      if ( y<0 || y>=N ) continue;
      if ( z<0 || z>=N ) continue;
      
      double w=1;
      h3(x,y,z) += w;
    }
    SP3time = appl::appl_timer_stop(SP3timer);
    // t3sum += SP3time/Nev;
    t3.add(SP3time/Nev);
  }

  //  t3sum /=Nloop;
  //  t3 = t3sum - t0;
  
  
  //  cout << "T0       = " << N << "  " << dufftime;
  cout << "T0       = " << N << "  " << t0.mean() << " +- " << t0.sigma();
  //  cout << "\tTH3 time = " << TH3time-dufftime;
  // cout << "\tSM3 time = " << SM3time-dufftime;
  // cout << "\tSP3 time = " << SP3time-dufftime << endl;
  cout << "\tTH3 time = " << t1.mean() << " +- " << t1.sigma();
  cout << "\tSM3 time = " << t2.mean() << " +- " << t2.sigma();
  cout << "\tSP3 time = " << t3.mean() << " +- " << t3.sigma() << endl;
  
  return;
}

