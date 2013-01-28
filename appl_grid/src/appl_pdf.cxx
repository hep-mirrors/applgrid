// emacs: this is -*- c++ -*-
//
//   appl_pdf.cxx        
//
//   pdf transform functions                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: appl_pdf.cxx, v1.0   Mon Dec 10 01:36:04 GMT 2007 sutt $


#include "appl_grid/appl_pdf.h" 
using namespace appl;

#include "appl_grid/mcfmz_pdf.h"
#include "appl_grid/mcfmw_pdf.h"
#include "appl_grid/mcfmQQ_pdf.h"
#include "appl_grid/nlojet_pdf.h"
#include "appl_grid/jetrad_pdf.h"
#include "appl_grid/nlojetpp_pdf.h"
#include "appl_grid/dis_pdf.h"
#include "appl_grid/vrapz_pdf.h"
#include "appl_grid/generic_pdf.h"


namespace appl { 

// initialise the map with some default instances
// although the user could create these themselves
// if they wanted
pdfmap appl_pdf::m_pdfmap; 

mcfmz_pdf     mcfmzpdf;
mcfmwp_pdf    mcfmwppdf;
mcfmwm_pdf    mcfmwmpdf;

mcfmCC_pdf    mcfmccpdf;
mcfmBB_pdf    mcfmbbpdf;
mcfmTT_pdf    mcfmttpdf;

nlojet_pdf    nlojetpdf;
nlojetpp_pdf  nlojetpppdf;
jetrad_pdf    jetradpdf;
dis_pdf       dispdf;
vrapzLO_pdf   vrapzLOpdf;
vrapzNLO_pdf  vrapzNLOpdf;
vrapzNNLO_pdf vrapzNNLOpdf;

generic_pdf   genericpdf;


};

