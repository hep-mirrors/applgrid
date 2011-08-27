//
//   @file    mcfm_grid.cxx         
//             
//
//   @author M.Sutton
// 
//   Copyright (C) 2011 M.Sutton (sutt@cern.ch)    
//
//   $Id: mcfm_grid.cxx, v0.0   Fri 12 Aug 2011 08:03:37 CEST sutt $



#include "appl_grid/mcfm_grid.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


/// fortran structures and common blocks...

static const int __nf__   = 5;
static const int __nf2__  = 11;
static const int __maxd__ = 41;


typedef struct {
  double weightfactor;
  double weightb [ __nf2__ ][ __nf2__ ];
  double weightv [ __nf2__ ][ __nf2__ ];
  double weightv1[ __nf2__ ][ __nf2__ ];
  double weightv2[ __nf2__ ][ __nf2__ ];
  double weightr [ __nf2__ ][ __nf2__ ][ __maxd__ ];
} __gridweight__;


typedef struct {
  double vsq[ __nf2__ ][ __nf2__ ], vsum[ __nf2__ ];
} __ckm__;


typedef struct  {
  double ag_xx1,ag_xx2,ag_x1z,ag_x2z,ag_scale,refwt,refwt2;
  int    contrib, dipole;
} __gridevent__;


typedef struct  { 
  int nproc;
} __nproc__;


typedef struct {
  int nflav;
} __nflav__;


extern "C" __ckm__ ckm_;
extern "C" __gridevent__  gridevent_;
extern "C" __gridweight__ gridweight_;
extern "C" __nproc__ nproc_;
extern "C" __nflav__ nflav_;


#if 0
/// these are defined as extern, why are these needed???
__ckm__ ckm_;
__gridweight__ gridweight_;
__gridevent__ gridevent_;
__nproc__ nproc_;
__nflav__ nflav_;
#endif


/// special mcfm specific methods

void appl::mcfm_grid::fillMCFM(double obs)
{
  //  std::cout << "   --------------------------------------------- " << std::endl;
  //  double* weight = new double[ subProcesses() ];

  std::vector<double> _weight( subProcesses() );

  double* weight = &_weight[0];

  double  scale2 =  gridevent_.ag_scale * gridevent_.ag_scale;
  
  //   std::cout <<" x1 = "<< gridevent_.ag_xx1
  // 	    <<" x2 = "<< gridevent_.ag_xx2
  // 	    <<" Q  = "<< gridevent_.ag_scale
  // 	    <<" CON = "<< gridevent_.contrib
  // 	    <<std::endl;
  
  int flag = 0;

  if ( gridevent_.contrib == 100 )  // BORN
    {
      flag = 0;
      collectWeight( gridevent_.contrib , flag, weight );
      
      if (isOptimised())  {
	fill( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 0 );
	getReference()->Fill( obs, gridevent_.refwt );
      }
      else { 
	fill_phasespace( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 0 );
      }
    }
  else if ( gridevent_.contrib == 200 ) { 
    //REAL
    
    collectWeight( gridevent_.contrib, gridevent_.dipole, weight );
      
    if (isOptimised()) {
      fill( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 1 );
      getReference()->Fill( obs, gridevent_.refwt );
    }
    else { 
      fill_phasespace( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 1 );
    }
  }
  else if ( gridevent_.contrib == 300 ) { 
    //VIRTUAL
   
    // BORN
    int flag = 0;
    collectWeight( gridevent_.contrib, flag, weight );
    if (isOptimised())
      fill( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 0 );
    else
      fill_phasespace( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 0 );
    
    //VIRT X1X2
    flag = -1;
    collectWeight( gridevent_.contrib, flag, weight );
    if (isOptimised())
      fill( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 1 );
    else
      fill_phasespace( gridevent_.ag_xx1, gridevent_.ag_xx2, scale2, obs, weight, 1 );
    
    //VIRT X1onZ
    flag = -2;
    collectWeight( gridevent_.contrib, flag, weight );
    if (isOptimised())
      fill( gridevent_.ag_x1z, gridevent_.ag_xx2, scale2, obs, weight, 1 );
    else
      fill_phasespace( gridevent_.ag_x1z, gridevent_.ag_xx2, scale2, obs, weight, 1 );
    
    //VIRT X2onZ
    flag = -3;
    collectWeight( gridevent_.contrib, flag, weight );
    if (isOptimised())
      fill( gridevent_.ag_xx1, gridevent_.ag_x2z, scale2, obs, weight, 1 );
    else
      fill_phasespace( gridevent_.ag_xx1, gridevent_.ag_x2z, scale2, obs, weight, 1 );
    
    if (isOptimised()) getReference()->Fill( obs, gridevent_.refwt );
    
  }
  else  { 
    // UNKNOWN
    std::cout <<__PRETTY_FUNCTION__<<" Unknown contribution!!!!!!!!!"<< std::endl; 
    std::cout <<__PRETTY_FUNCTION__<<" Unknown contribution!!!!!!!!!"<< std::endl; 
    std::cout <<__PRETTY_FUNCTION__<<" Unknown contribution!!!!!!!!!"<< std::endl; 
  }
  
  //  std::cout << "   --------------------------------------------- " << std::endl;
  
  //  delete [] weight;
  return;
}


void appl::mcfm_grid::collectWeight(const int& order, const int& id, double* wt)
{
  double factor = 1.0;
  int iproc = -1;
  for ( int i = 0 ;i <subProcesses(order); i++ ) wt[i] =  0.0;

  //   std::cout<<" \n( ";
  //   for (int jj =  0 ; jj < subProcesses(order);jj++) std::cout << wt[jj] <<" , ";
  //   std::cout<<")\n";
  
  for ( int iflav = -nflav_.nflav; iflav <= nflav_.nflav; iflav++)
    for ( int jflav = -nflav_.nflav; jflav <= nflav_.nflav; jflav++)
      {

	decideSubProcess( iflav, jflav, iproc, factor);
	if ( iproc < 0 ) continue;
	
	if ( order == 100 ) 
	  wt[ iproc ] +=  factor * gridweight_.weightb[ __nf__ + jflav ][ __nf__ + iflav ];
	if ( order == 200 ) 
	  wt[ iproc ] +=  factor * gridweight_.weightr[ __nf__ + jflav ][ __nf__ + iflav ][ id ];

	
	// 	if (true)
	// 	  {
	// 	    if (0 != gridweight_.weightb[ __nf__ + jflav ][ __nf__ + iflav ]) 
	// 	      std::cout <<" ( i= "<<iflav<<" j= "<<jflav<<" : "
	// 			<< gridweight_.weightb[ __nf__ + jflav ][ __nf__ + iflav ]<<" ) proc = "<<iproc;
	
	// 	    std::cout<<" \n( ";
	// 	    for (int jj =  0 ; jj < subProcesses(order);jj++) std::cout << wt[jj] <<" , ";
	// 	    std::cout<<")\n";
	
	// 	  }

	if ( order == 300 )
	  {
	    if      ( id ==  0 )
	      wt[ iproc ] +=  factor * gridweight_.weightb [ __nf__ + jflav ][ __nf__ + iflav ];
	    else if ( id == -1 )
	      wt[ iproc ] +=  factor * gridweight_.weightv [ __nf__ + jflav ][ __nf__ + iflav ];
	    else if ( id == -2 )
	      wt[ iproc ] +=  factor * gridweight_.weightv1[ __nf__ + jflav ][ __nf__ + iflav ];
	    else if ( id == -3 )
	      wt[ iproc ] +=  factor * gridweight_.weightv2[ __nf__ + jflav ][ __nf__ + iflav ];
	  }
	

      }

  for (int jj=0 ; jj<subProcesses(order) ; jj++ )  wt[jj] *= gridweight_.weightfactor ;

  //   std::cout<<" factor = "<< gridweight_.weightfactor<<std::endl;
  
  //   std::cout << "PROC = " << nproc_.nproc <<" W = ( dip = "<<id<<" , ";
  //   for (int jj =  0 ; jj < subProcesses(order);jj++) std::cout << wt[jj] <<" , ";
  //   std::cout <<" ) if = "<<factor<<" psw = "<< gridweight_.weightfactor <<" refwt = "<<gridevent_.refwt  << std::endl;
  //   std::cout <<" me(2,-1,0)"<< gridweight_.weightr[ __nf__ - 1][ __nf__ + 2][0] <<std::endl;
  //   std::cout <<" me(2,-1,1)"<< gridweight_.weightr[ __nf__ - 1][ __nf__ + 2][1] <<std::endl;
  
  return ;
}


void appl::mcfm_grid::decideSubProcess(const int& iflav1, const int& iflav2, int & iProcess , double &factor )
{
  iProcess = -1;
  factor   = 0.;

  if ( nproc_.nproc == 1 )
    {
      if( (iflav1 == 0) && (iflav2 ==  2) ) 
	{
	  factor = 1.0/ckm_.vsum[ __nf__ + iflav2];
	  iProcess = 5;
	}        
      else if( (iflav1 == 0) && (iflav2 == -1) )
	{
	  factor = 1.0/ckm_.vsum[ __nf__ + iflav2];
	  iProcess = 4;
	}
      else if( (iflav2 == 0) && (iflav1 == 2) )
	{
	  factor = 1.0/ckm_.vsum[ __nf__ + iflav1];
	  iProcess = 3;
	}
      else if( (iflav2 == 0) && (iflav1 == -1) )
	{
	  factor = 1.0/ckm_.vsum[ __nf__ + iflav1];
	  iProcess = 2;
	}
      else if( (iflav1 == 2) && (iflav2 == -1) )
	{
	  factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + iflav1];
	  iProcess = 1;
	}
      else if( (iflav1 == -1) && (iflav2 == 2) )
	{
	  factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + iflav1];
	  iProcess = 0;
	}
    }
  else if ( nproc_.nproc == 6 )
    {

      if((iflav1 == 0) && (iflav2 == -2))
	{
	  factor = 1.0/ckm_.vsum[__nf__ + iflav2];
	  iProcess = 5;
	}
      else if((iflav1 == 0) && (iflav2 == 1))
	{
	  factor = 1.0/ckm_.vsum[__nf__ + iflav2];
	  iProcess = 4;
	}
      else if((iflav2 == 0) && (iflav1 == -2))
	{
	  factor = 1.0/ckm_.vsum[__nf__ + iflav1];
	  iProcess = 3;
	}
      else if((iflav2 == 0) && (iflav1 == 1))
	{
	  factor = 1.0/ckm_.vsum[__nf__ + iflav1];
	  iProcess = 2;
	}
      else if ((iflav1 == -2) && (iflav2 == 1))
	{
	  factor = 1.0/ckm_.vsq[__nf__ + iflav2][__nf__ + iflav1];
	  iProcess = 1;
	}
      else if ((iflav1 == 1) && (iflav2 == -2))
	{
	  factor = 1.0/ckm_.vsq[__nf__ + iflav1][ __nf__ + iflav2];
	  iProcess = 0;
	}
    }
  else if ( nproc_.nproc == 31 )
    {
      if      ( iflav2 == 0 )
	{
	  if      (iflav1 == -1) {iProcess = 11;}
	  else if (iflav1 ==  1) {iProcess = 10;}
	  else if (iflav1 == -2) {iProcess = 9;}
	  else if (iflav1 ==  2) {iProcess = 8;}
	}
      else if ( iflav1 == 0 )
	{
	  if      (iflav2 == -1) iProcess = 7;
	  else if (iflav2 ==  1) iProcess = 6;
	  else if (iflav2 == -2) iProcess = 5;
	  else if (iflav2 ==  2) iProcess = 4;
	}
      else if ( (iflav1 != 0 ) && ( iflav2 != 0 ) )
	{
	  if      (iflav1 == -1) iProcess = 3;
	  else if (iflav1 == -2) iProcess = 2;
	  else if (iflav1 ==  1) iProcess = 1;
	  else if (iflav1 ==  2) iProcess = 0;
	}
      factor = 1.0;
    }
  //   else
  //     {
  //       iProcess = -1;
  //       factor   = 0.;
  //     }
  //  std::cout << "\t\t *** \t"<<iflav1<<" <> "<<iflav2<<" iproc = "<<iProcess<< std::endl; 
  
  return;
}

