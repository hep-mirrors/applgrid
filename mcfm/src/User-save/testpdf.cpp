#include <iostream>
#include <string>
#include <math.h>
#include "TFile.h"
#include "APPLgrid.h"
#include "APPLpdf.h"
//
//
//
extern "C" void dglapeval_(const double *, const double *, double *);
//extern "C" void dglapevalsplit_(const double *, const double *,  const int *, const int *, double *);
extern "C" void initmypdf_(const char *, const int *);
//extern "C" double alphaspdf_(const double *);

void make_pdf_combinations(const double* fA, const double* fB, double* H);

int main(int argc, char **argv) 
{
  //
  // LHAPDF part
  //
  APPL::grid* mygrid = new APPL::grid("weightgrid.root");
  APPL::pdf *mypdf = new APPL::pdf();


  double* H;                    // pdf for subprocess

  //
  int stop;
  while(stop != 0)
    {
      double x1,x2,scale;
      cout<<"Enter x1"<<endl;
      cin>>x1;
      cout<<"Enter x2"<<endl;
      cin>>x2;
      cout<<"Enter scale^2"<<endl;
      cin>>scale;

      H = mypdf->calculate(x1,x2,std::sqrt(scale));
      
      cout<<" H = (";
      for(int isub = 0; isub < mygrid->Nsub(); isub++) 
	{
	  cout<<H[isub]<<" , ";
	}
      cout<<" ) "<<endl;
      cout<<"Stop program? 0 = Yes 1 = No "<<endl;
      cin>>stop;
    }
}

