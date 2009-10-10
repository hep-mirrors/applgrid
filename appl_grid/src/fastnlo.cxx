//
//   @file    fastbnlo.cxx         
//   
//      Implementaion of the fastnlo - applgrid interface
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: fastnlo.cxx, v0.0   Sun 27 Sep 2009 22:39:32 BST sutt $


#include "appl_grid/fastnlo.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>


std::ofstream sout( "/dev/null" );
// std::ostream& sout = std::cout;

template<class A>
void out( A& a ) { sout << a << std::endl; }

template<class A>
void out( const std::string& s, A& a ) { sout << s << ":\t " << a << std::endl; }






fastnlo::fastnlo( const std::string& filename ) {
  
  std::cout << "fastnlo::fastnlo() reading fastnlo grid file " << filename << std::endl;

  // can we open the fastnlo grid file
  struct stat stfileinfo;
  if ( stat(filename.c_str(),&stfileinfo) )   {
    std::cerr << "fastnlo::fastnlo() cannot open fastnlo grid file " << filename << std::endl;
    return;
  }

  // read the fastnlo grid
  std::ifstream faststream( filename.c_str() );

  std::string dummy;

  std::string seperator = "1234567890";


  faststream >> dummy;  

  long long tabversion=0;   faststream >> tabversion; out( "tabversion", tabversion);

  faststream >> dummy; 

  long long ireaction=0;   faststream >> ireaction;   out( "reaction", ireaction );
  double    Ecms=0;        faststream >> Ecms;        out( "Ecms", Ecms );
  
  double ixsecunits = 0;  faststream >> ixsecunits;   out( "ixsecunits", ixsecunits );
  std::string label;      faststream >> label;        out( "label", label );

  faststream >> dummy;  
  faststream >> dummy;
  faststream >> dummy;
  faststream >> dummy; 

  long long iproc;       faststream >> iproc;    out( "iproc", iproc );
  long long ialgo;       faststream >> ialgo;    out( "ialgo", ialgo );
  double    jetres1;     faststream >> jetres1;  out( "jetres1", jetres1 );
  double    jetres2;     faststream >> jetres2;  out( "jetres2", jetres2 );
  long long Norder = 0;  faststream >> Norder;   out( "Norder", Norder );

  std::vector<int> Npow(Norder); 
  for ( int i=0 ; i<Norder ; i++ ) { faststream >> Npow[i]; out( "\t", Npow[i] ); }

  std::vector<std::string> powlabel(Norder);
  std::vector<std::string> codelabel(Norder);
  for ( int i=0 ; i<Norder ; i++ ) { 
    faststream >> powlabel[i];   out( "\tpowlabel:",  powlabel[i] ); 
    faststream >> codelabel[i];  out( "\tcodelabel:", codelabel[i] ); 
  }

  faststream >> dummy;  

  std::vector<long long> Nevt(Norder);   
 
  for ( int i=0 ; i<Norder ; i++ ) { faststream >> Nevt[i]; out( "\t:", Nevt[i] ); }

  int Nxtot    = 0;   faststream >> Nxtot;     out( "Nxtot", Nxtot );
  int ixscheme = 0;   faststream >> ixscheme;  out( "scheme", ixscheme );
  int ipdfwgt  = 0;   faststream >> ipdfwgt;   out( "ipdfwgt", ipdfwgt );
  int iref     = 0;   faststream >> iref;      out( "iref", iref );

  faststream >> dummy; 

  int Nbintot = 0;     faststream >> Nbintot;     out( "Nbintot", Nbintot );
  int Ndimension = 0;  faststream >> Ndimension;  out( "Ndimension", Ndimension );

  std::vector<std::string> dimlabel(Ndimension);
  for ( int i=0 ; i<Ndimension ; i++ ) { faststream >> dimlabel[i]; out( "dimlabel: ", dimlabel[i] ); }

  int Nrapidity = 0;  faststream >> Nrapidity; out( "Nrapidity", Nrapidity );

  // rapidity bins

  std::vector<double> rapbin(Nrapidity+1);
  for ( int i=0 ; i<Nrapidity+1 ; i++ ) { faststream >> rapbin[i]; out( "\trapbin:", rapbin[i] ); }

  std::vector<int>  Npt(Nrapidity); out( "\nNpt:", Nrapidity );
  for ( int i=0 ; i<Nrapidity ; i++ ) {  faststream >> Npt[i]; out( "\tNpt:", Npt[i] );  }


  // pt bins   
  std::vector<std::vector<double> > ptbin(Nrapidity);
  for ( int i=0 ; i<Nrapidity ; i++ ) {
    ptbin[i] = std::vector<double>(Npt[i]+1);
    for ( int j=0 ; j<Npt[i]+1 ; j++ )  { 
      faststream >> ptbin[i][j]; 
      //      std::cout << "\t\tptbin " << ptbin[i][j] << std::endl;
    }
  }

  faststream >> dummy; 
  
  // xlimits 
  std::vector<std::vector<double> > xlimit(Nrapidity);
  for ( int i=0 ; i<Nrapidity ; i++ ) {
    xlimit[i]=std::vector<double>(Npt[i]);
    for ( int j=0 ; j<Npt[i] ; j++ ) faststream >> xlimit[i][j]; 
  }

  
  faststream >> dummy; 

  std::string scalelabel;  faststream >> scalelabel; out( "scalelabel: ", scalelabel );  
  int Nscalebin = 0;       faststream >> Nscalebin;  out( "Nscalebin: ", Nscalebin );  

  std::vector<std::vector<std::vector<double> > > murval(Nrapidity);
  //  murval.reserve(Nrapidity);
  for ( int i=0 ; i<Nrapidity ; i++ ) {
    murval[i] = std::vector<std::vector<double> >(Npt[i]);
    for ( int j=0 ; j<Npt[i] ; j++ ) { 
      murval[i][j]=std::vector<double>(Nscalebin);
      //      std::cout << "rapidity bin " << i << "\n\tpt " << j << "\tNscalebin " << Nscalebin << "\n\t";
      for ( int k=0 ; k<Nscalebin ; k++ ) { 
	faststream >> murval[i][j][k];
	//	std::cout << "\t" << murval[i][j][k];
      }
      //      std::cout << std::endl;
    }
  }
  
  faststream >> dummy; 

  std::vector<std::vector<std::vector<double> > > mufval(Nrapidity);
  for ( int i=0 ; i<Nrapidity ; i++ ) {
    mufval[i] = std::vector<std::vector<double> >(Npt[i]);
    for ( int j=0 ; j<Npt[i] ; j++ ) { 
      mufval[i][j] = std::vector<double>(Nscalebin);
      for ( int k=0 ; k<Nscalebin ; k++ ) faststream >> mufval[i][j][k];
    }
  }
  
  faststream >> dummy; 

  int Nscalevar = 0;  faststream >> Nscalevar;  out( "Nscalevar: ", Nscalevar );

  std::vector<double> murscale(Nscalevar);
  for ( int i=0 ; i<Nscalevar ; i++ ) { faststream >> murscale[i]; out( "mur: ", murscale[i] ); }

  std::vector<double> mufscale(Nscalevar);
  for ( int i=0 ; i<Nscalevar ; i++ ) { faststream >> mufscale[i]; out ( "muf: ", mufscale[i] ); }

  faststream >> dummy;  
  

  int Nxsum = 0;
  int Nsubproc = 0;

 
  if ( (ireaction==2) || (ireaction==3) ) { //  pp or ppbar
    Nxsum = (Nxtot*Nxtot+Nxtot)/2;
    Nsubproc = 7;
  }
  else { 
    if (ireaction==1) { //  ! DIS
      Nxsum = Nxtot;
      Nsubproc = 3;
    }
  }

  out("Nxsum ", Nxsum );
  out("Nsubproc ", Nsubproc );


  out( "Nrapidity: ", Nrapidity );
  

  int Nbins = 0;
  for ( int i=0 ; i<Nrapidity ; i++ ) Nbins += Npt[i];
  
  std::vector<  std::vector<  std::vector<  std::vector<  std::vector<double> > > > > array(Nbins);

  int ibin = 0;
  for ( int irap=0 ; irap<Nrapidity ; irap++ ) { 
    // out( "Npt: ", Npt[irap] );
    for ( int ipt=0 ; ipt<Npt[irap] ; ipt++ ) {
      // out( "Nxsum: ", Nxsum );
      
      array[ibin] = std::vector<  std::vector<  std::vector<  std::vector<double> > > >(Nxsum);

      int ix = 0;
      for ( int ix1=0 ; ix1<Nxtot ; ix1++ ) { 
	for ( int ix2=0 ; ix2<=ix1 ; ix2++ ) { 
	  
	  // out( "Nsubproc: ", Nsubproc );

	  array[ibin][ix] = std::vector<  std::vector<  std::vector<double> > >(Nsubproc);
      
	  for ( int isub=0 ; isub<Nsubproc ; isub++ ) {
	    long long ns=1+Nscalevar*(Norder-1);
	    //	  out( "Nscales: ", ns );
	    
	    array[ibin][ix][isub] = std::vector<  std::vector<double> >(ns);
	    
	    for ( int iscale=0 ; iscale<ns ; iscale++ ) { 
	      // out ( "Nscalebin: ", Nscalebin );
	      
	      array[ibin][ix][isub][iscale] = std::vector<double>(Nscalebin);
	      
	      for ( int i=0 ; i<Nscalebin ; i++ ) { 
		faststream >> array[ibin][ix][isub][iscale][i];

#if 0
		std::cout << "array " 
			  << ibin << " " << ix << " " << isub << " " << iscale << " " << i << " " 
			  << " : " << array[ibin][ix][isub][iscale][i] << std::endl;
	
#endif	

	      }
	    }
	  }

	  ix++;
	}
      }
      ibin++;
    }
  }
    
  faststream >> dummy; 


  std::cout << "making grids..." << std::endl;
  
  // get the central bin since we use our own scale variations
  // we don't have seperate grids for each different scale
  
  int icentral = 0;
  for ( int imu=0 ; imu<Nscalevar ; imu++ ) { 
    if ( std::fabs(murscale[imu]-1)<1e-5 ) icentral = imu;
  } 
  
  
  int iposition[3] = { 0, 1+icentral, 1+icentral+Nscalevar };
  ibin = 0;
  


  //  appl::grid** m_grid = new appl::grid*[Nrapidity];
  m_grid = std::vector<appl::grid*>(Nrapidity);

  // need to implement other schemes if need be...
  std::string transform;
  if ( ixscheme==2 ) transform = "f3";
  

  std::string pdfname; 
  if ( ireaction==2 ) pdfname = "nlojet";
  if ( ireaction==3 ) pdfname = "nlojetpp";


  for ( int irap=0 ; irap<Nrapidity ; irap++ ) { 

    //    std::cout << "rap " << irap << std::endl;
    
    std::vector<double> ptlims;
    
    for ( int ipt=0 ; ipt<Npt[irap]+1 ; ipt++ ) {
      ptlims.push_back(ptbin[irap][ipt]);
    }

    
    m_grid[irap] = new appl::grid( ptlims, pdfname, Npow[0], Npow.back()-Npow[0], transform );
    m_grid[irap]->symmetrise(true);
    m_grid[irap]->setCMSScale(Ecms);
    
    for ( int ipt=0 ; ipt<Npt[irap] ; ipt++ ) {
      
      // need this width, since the grid is already divided by it, 
      // so when the applgrid normalises to this, it will be wrong.
      double width = ptbin[irap][ipt+1]-ptbin[irap][ipt];

      //      std::cout << "width " << ipt << "\t" << width << std::endl;

      // need to generalise this??
      double hxlim = -std::sqrt( -std::log10( xlimit[irap][ipt] ) );
      
      double hxl = hxlim;
      double hxu = hxlim/Nxtot;
      double xl = std::pow(10, -hxl*hxl);
      double xu = std::pow(10, -hxu*hxu);
      
      //      std::cout << "xl " << xl << "\txu " << xu << std::endl; 
      
      // this has made the x values - only need to make the limits, 
      // using transform f3
      
      for ( int iord=0 ; iord<Norder ; iord++ ) { 

	//	std::cout << "Rap " << irap << "  pt " << ipt << "  order " << iord << std::endl;

	double Q2min = murval[irap][ipt][0];            Q2min *= Q2min;
	double Q2max = murval[irap][ipt][Nscalebin-1];  Q2max *= Q2max;

	//	std::cout << "Q2 " << Q2min << "\t" << Q2max << "\tNscale " << Nscalebin << std::endl;  


	appl::igrid* g = new appl::igrid( Nscalebin, Q2min, Q2max, 0, 
					  Nxtot, xl, xu, 0, 
					  transform, Nsubproc );
	
	g->symmetrise(true);
	if ( ipdfwgt ) g->reweight(true);

	m_grid[irap]->add_igrid(ipt, iord, g);

	//	std::cout << "grid ";
	for ( int isub=0 ; isub<Nsubproc ; isub++ ) { 
	  
	  
	  SparseMatrix3d& sgrid = *g->weightgrid()[isub];
	  
	  
 	  // now need to fill the grid;
	  
	  int ix=0;
	  
	  int NQ  = g->weightgrid()[isub]->Nx();
	  int Nx1 = g->weightgrid()[isub]->Ny();
	  int Nx2 = g->weightgrid()[isub]->Nz();
	  
	  //	  std::cout << "ord " << iord << "\tipos " << iposition[iord] << std::endl; 
	  
	  for ( int ix1=0 ; ix1<Nxtot ; ix1++ ) { 
	    for ( int ix2=0 ; ix2<=ix1 ; ix2++ ) { 
	      for ( int iQbin=0 ; iQbin<Nscalebin ; iQbin++ ) {       
		// leading order
		//		double d = array[ibin][ix][isub][iposition[iord]][iQbin];
		double d = array[ibin][ix][isub][iposition[iord]][iQbin];
#if 0
		std::cout << "\taray[" << ibin << "][" << ix << "][" << isub << "][" << iposition[iord] << "]" << std::endl;
		std::cout << "  Q " << iQbin << "  x1 " <<  ix1 << "  x1 " << ix2 << "  ix " << ix 
			  << "  g " << d << std::endl;
		
		std::cout << "array " 
			  << ibin << " " << ix << " " << isub << " " << iposition[iord] << " " << iQbin << " " 
			  << " : " << array[ibin][ix][isub][iposition[iord]][iQbin] << std::endl;
#endif
		
		sgrid(iQbin, Nx1-1-ix1, Nx2-1-ix2) = d*width;
	      }
	      ix++;
	    }
	  }
	  // loops over x1, x2 and Q2 bins
	  
	} // Nsubproc
      }

      ibin++;	
      
    }
  }
  
}






