//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include <cstring>

//   nlo includes
#include "bits/proc-hhc.h"
#include "bits/nlo-color.h"
#include "bits/nlo-flkern.h"

#include "ampg4.h"
#include "ampg5.h"
#include "ampg6.h"

#include "ampq2g2.h"
#include "ampq2g3.h"
#include "ampq2g4.h"

#include "ampq4.h"
#include "ampq4g1.h"
#include "ampq4g2.h"

#include "ampq6.h"




namespace nlo {
  

  //  0. gg  -> gg    (1/2 * 1/Na * 1/Na)
  //  1. gg  -> qqb   (1/Na * 1/Na * nf)
  //  2. qg  -> qg    (1/Nc * 1/Na)
  //  3. gq  -> qg    (1/Na * 1/Nc)
  //  4. qqb -> gg    (1/2 * 1/Nc * 1/Nc)
  //  5. qr  -> qr    (1/Nc * 1/Nc)
  //  6. qq  -> qq    (1/2 * 1/Nc * 1/Nc)
  //  7. qqb -> rrb   (1/Nc * 1/Nc * (nf-1))
  //  8. qqb -> qqb   (1/Nc * 1/Nc)
  //  9. qrb -> qrb   (1/Nc * 1/Nc)
  void _hhc_jet_base::
  amp_tree(ampg4 *amp1, ampq2g2 *amp2, ampq4 *amp3, double *subp, double *out)
  {
    static double __sub[10];
    double *sub = (subp ? subp : __sub);
    double tmp[2], nf = Nf, nf1 = Nf-1;
    std::memset(out, 0, 7*sizeof(double));

    //----- FOUR GLUONS PROCESS -----
    if(amp1) {
      //  gg -> gg    (1/2 * 1/Na * 1/Na)
      out[0] += (sub[0] = (amp1 -> su3_tree(-1,0,1,2))/(2.0*Na2));
    }
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    if(amp2) {
      //  gg -> qqb   (1/Na * 1/Na * nf)
      out[0] += (sub[1] = nf*(amp2 -> su3_tree(1,2,-1,0))/Na2);
      
      //  qg -> qg    (1/Nc * 1/Na)
      out[1] += (sub[2] = -(amp2 -> su3_tree(1,-1,2,0))/(Nc*Na));
      
      //  gq -> qg    (1/Na * 1/Nc)
      out[2] += (sub[3] = -(amp2 -> su3_tree(1,0,2,-1))/(Na*Nc));
      
      //  qqb -> gg   (1/2 * 1/Nc * 1/Nc)
      out[5] += (sub[4] = (amp2 -> su3_tree(0,-1,1,2))/(2.0*Nc2));
    }
    
    //----- FOUR QUARKS PROCESSES -----
    if(amp3) {     
      //  qr -> qr    (1/Nc * 1/Nc)
      //  qq -> qq    (1/2 * 1/Nc * 1/Nc)
      amp3 -> su3_tree(1,-1,2,0, sub+5);
      out[3] += (sub[5] /= Nc2);
      out[4] += (sub[6] /= 2.0*Nc2);
      
      //  qqb -> rrb  (1/Nc * 1/Nc * (nf-1))
      //  qqb -> qqb  (1/Nc * 1/Nc)
      amp3 -> su3_tree(0,-1,1,2, sub+7); 
      out[5] += (sub[7] *= nf1/Nc2) + (sub[8] /= Nc2);
      
      //  qrb -> qrb  (1/Nc * 1/Nc)
      amp3 -> su3_tree(1,-1,0,2, tmp);
      out[6] += (sub[9] = tmp[0]/Nc2);
    }
  }

  void _hhc_jet_base::
  amp_1loop(ampg4 *amp1, ampq2g2 *amp2, ampq4 *amp3, double *out)
  {
    double tmp[2], nf = Nf, nf1 = Nf-1;

    //----- FOUR GLUONS PROCESS -----
    //  gg -> gg    (1/2 * 1/Na * 1/Na)
    out[0] = (amp1 -> su3_1loop(Nf, -1,0,1,2))/(2.0*Na2);
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    //  gg -> qqb   (1/Na * 1/Na * nf)
    out[0] += nf*(amp2 -> su3_1loop(1,2,-1,0))/Na2;
    
    //  qg -> qg    (1/Nc * 1/Na)
    out[1] = (amp2 -> su3_1loop(1,-1,2,0))/(Nc*Na);
    
    //  gq -> qg    (1/Na * 1/Nc)
    out[2] = (amp2 -> su3_1loop(1,0,2,-1))/(Na*Nc);
    
    //  qqb -> gg   (1/2 * 1/Nc * 1/Nc)
    out[5] = (amp2 -> su3_1loop(0,-1,1,2))/(2.0*Nc2);

    //----- FOUR QUARKS PROCESSES -----
    //  qr -> qr    (1/Nc * 1/Nc)
    //  qq -> qq    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_1loop(Nf, 1,-1,2,0, tmp);
    out[3] = tmp[0]/Nc2;
    out[4] = tmp[1]/(2.0*Nc2);
    
    //  qqb -> rrb  (1/Nc * 1/Nc * (nf-1))
    //  qqb -> qqb  (1/Nc * 1/Nc)
    amp3 -> su3_1loop(Nf, 0,-1,1,2, tmp); 
    out[5] += tmp[0]*nf1/Nc2 + tmp[1]/Nc2;
    
    //  qrb -> qrb  (1/Nc * 1/Nc)
    amp3 -> su3_1loop(Nf, 1,-1,0,2, tmp);
    out[6] = tmp[0]/Nc2;
  }

  void _hhc_jet_base::
  amp_1loop_mch(ampg4 *amp1, ampq2g2 *amp2, ampq4 *amp3, double *out)
  {
    double tmp[2], nf = Nf, nf1 = Nf-1;

    //----- FOUR GLUONS PROCESS -----
    //  gg -> gg    (1/2 * 1/Na * 1/Na)
    out[0] = (amp1 -> su3_1loop_mch(Nf, -1,0,1,2))/(2.0*Na2);
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    //  gg -> qqb   (1/Na * 1/Na * nf)
    out[0] += nf*(amp2 -> su3_1loop_mch(1,2,-1,0))/Na2;
    
    //  qg -> qg    (1/Nc * 1/Na)
    out[1] = (amp2 -> su3_1loop_mch(1,-1,2,0))/(Nc*Na);
    
    //  gq -> qg    (1/Na * 1/Nc)
    out[2] = (amp2 -> su3_1loop_mch(1,0,2,-1))/(Na*Nc);
    
    //  qqb -> gg   (1/2 * 1/Nc * 1/Nc)
    out[5] = (amp2 -> su3_1loop_mch(0,-1,1,2))/(2.0*Nc2);

    //----- FOUR QUARKS PROCESSES -----
    //  qr -> qr    (1/Nc * 1/Nc)
    //  qq -> qq    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_1loop_mch(Nf, 1,-1,2,0, tmp);
    out[3] = tmp[0]/Nc2;
    out[4] = tmp[1]/(2.0*Nc2);
    
    //  qqb -> rrb  (1/Nc * 1/Nc * (nf-1))
    //  qqb -> qqb  (1/Nc * 1/Nc)
    amp3 -> su3_1loop_mch(Nf, 0,-1,1,2, tmp); 
    out[5] += tmp[0]*nf1/Nc2 + tmp[1]/Nc2;
    
    //  qrb -> qrb  (1/Nc * 1/Nc)
    amp3 -> su3_1loop_mch(Nf, 1,-1,0,2, tmp);
    out[6] = tmp[0]/Nc2;
  }

  void _hhc_jet_base::
  amp_ccgg(ampg4 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[0] = (amp -> su3_cc(i,j, -1,0,p1,p2))/Na2;
  }

  void _hhc_jet_base::
  amp_ccgg(ampq2g2 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[0] = (amp -> su3_cc(i,j, p1,p2,-1, 0))/Na2;
  }

  void _hhc_jet_base::
  amp_ccag(ampq2g2 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[1] = -(amp -> su3_cc(i,j, -1,p1,p2,0))/(Na*Nc);
  }
  
  void _hhc_jet_base::
  amp_ccga(ampq2g2 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[2] = -(amp -> su3_cc(i,j, 0,p1,p2,-1))/(Na*Nc);
  }
  
  void _hhc_jet_base::
  amp_ccqg(ampq2g2 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[1] = -(amp -> su3_cc(i,j, p1,-1,p2,0))/(Na*Nc);
  }
  
  void _hhc_jet_base::
  amp_ccgq(ampq2g2 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[2] = -(amp -> su3_cc(i,j, p1,0,p2,-1))/(Na*Nc);
  }

  void _hhc_jet_base::
  amp_ccqa(ampq2g2 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[5] = (amp -> su3_cc(i,j, -1,0,p1,p2))/Nc2;
  }

  void _hhc_jet_base::
  amp_ccaq(ampq2g2 *amp, int i, int j, int p1, int p2, double *cc) {
    cc[5] = (amp -> su3_cc(i,j, 0,-1,p1,p2))/Nc2;
  }

  void _hhc_jet_base::
  amp_cc(ampq4 *amp, int i, int j, int p1, int p2, int p3, int p4, double *cc)
  {
    amp -> su3_cc(i,j, p1,p2,p3,p4, cc);
    cc[0] /= Nc2;
    cc[1] /= Nc2;
  }
  
  void _hhc_jet_base::
  amp_kp(double al, ampg4 *amp1, ampq2g2 *amp2, ampq4 *amp3, su3_kp_i2 *out)
  {
    static su3_kp_i2 tmp[2];
    double nf = Nf, nf1 = Nf-1;
    
    //----- FOUR GLUONS PROCESS -----
    //  gg -> gg    (1/2 * 1/Na * 1/Na)
    amp1 -> su3_kp(Nf, -1,0, -1,0,1,2, tmp, al); 
    out[0] = tmp[0]/(2.0*Na2);
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    //  gg -> qqb   (1/Na * 1/Na * nf)
    amp2 -> su3_kp(Nf, -1,0, 1,2,-1,0, tmp, al);
    out[0] += tmp[0]*(nf/Na2);
    
    //  qg -> qg    (1/Nc * 1/Na)
    amp2 -> su3_kp(Nf, -1,0, 1,-1,2,0, tmp, al);
    out[1] = tmp[0]/(-Nc*Na);
    
    //  gq -> qg    (1/Na * 1/Nc)
    amp2 -> su3_kp(Nf, -1,0, 1,0,2,-1, tmp, al);
    out[2] = tmp[0]/(-Na*Nc);
    
    //  qqb -> gg   (1/2 * 1/Nc * 1/Nc)
    amp2 -> su3_kp(Nf, -1,0, 0,-1,1,2, tmp, al);
    out[5] = tmp[0]/(2.0*Nc2);
    
    //----- FOUR QUARKS PROCESSES -----
    //  qr -> qr    (1/Nc * 1/Nc)
    //  qq -> qq    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_kp(-1,0, 1,-1,2,0, tmp, al);
    out[3] = tmp[0]/Nc2;
    out[4] = tmp[1]/(2.0*Nc2);
    
    //  qqb -> rrb  (1/Nc * 1/Nc * (nf-1))
    //  qqb -> qqb  (1/Nc * 1/Nc)
    amp3 -> su3_kp(-1,0, 0,-1,1,2, tmp, al); 
    out[5] += tmp[0]*(nf1/Nc2) + tmp[1]/Nc2;
    
    //  qrb -> qrb  (1/Nc * 1/Nc)
    amp3 -> su3_kp(-1,0, 1,-1,0,2, tmp, al);
    out[6] = tmp[0]/Nc2;
  }
 

  //  0.  gg  -> ggg    (1/6 * 1/Na * 1/Na)
  //  1.  gg  -> qqbg   (1/Na * 1/Na * nf)
  //  2.  qg  -> qgg    (1/2 * 1/Nc * 1/Na)
  //  3.  gq  -> qgg    (1/2 * 1/Na * 1/Nc)
  //  4.  qqb -> ggg    (1/6 * 1/Nc * 1/Nc)
  //  5.  qg  -> rrbq   (1/Nc * 1/Na * (nf-1))
  //  6.  qg  -> qqbq   (1/2 * 1/Nc * 1/Na)
  //  7.  gq  -> rrbq   (1/Na * 1/Nc * (nf-1))
  //  8.  gq  -> qqbq   (1/2 * 1/Na * 1/Nc)
  //  9.  qr  -> qrg    (1/Nc * 1/Nc)
  //  10. qq  -> qqg    (1/2 * 1/Nc * 1/Nc)
  //  11. qqb -> rrbg   (1/Nc * 1/Nc * (nf-1))
  //  12. qqb -> qqbg   (1/Nc * 1/Nc)
  //  13. qrb -> qrbg   (1/Nc * 1/Nc)
  void _hhc_jet_base::
  amp_tree(ampg5 *amp1, ampq2g3 *amp2, ampq4g1 *amp3, double *subp, double *out)
  {
    static double __sub[14];
    double *sub = (subp ? subp : __sub);
    double tmp[2], nf = Nf, nf1 = Nf-1;
    std::memset(out, 0,  7*sizeof(double));

    //----- FIVE GLUONS PROCESS -----
    if(amp1) {
      //  gg -> ggg    (1/6 * 1/Na * 1/Na)
      out[0] += (sub[0] = (amp1 -> su3_tree(-1,0,1,2,3))/(6.0*Na2));
    }
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    if(amp2) {
      //  gg -> qqbg   (1/Na * 1/Na * nf)
      out[0] += (sub[1] = nf*(amp2 -> su3_tree(1,2,3,-1,0))/Na2);
      
      //  qg -> qgg    (1/2 * 1/Nc * 1/Na)
      out[1] += (sub[2] = -(amp2 -> su3_tree(1,-1,2,3,0))/(2.0*Nc*Na));
      
      //  gq -> qgg    (1/2 * 1/Na * 1/Nc)
      out[2] += (sub[3] = -(amp2 -> su3_tree(1,0,2,3,-1))/(2.0*Na*Nc));
      
      //  qqb -> ggg   (1/6 * 1/Nc * 1/Nc)
      out[5] += (sub[4] = (amp2 -> su3_tree(0,-1,1,2,3))/(6.0*Nc2));
    }
    
    //----- FOUR QUARKS ONE GLUON PROCESSES -----
    if(amp3) {     
      //  qg -> rrbq   (1/Nc * 1/Na * (nf-1))
      //  qg -> qqbq   (1/2 * 1/Nc * 1/Na)
      amp3 -> su3_tree(3,-1,1,2,0, sub+5);
      out[1] += (sub[5] *= -nf1/(Nc*Na)) + (sub[6] /= -2.0*Nc*Na);
      
      //  gq -> rrbq   (1/Na * 1/Nc * (nf-1))
      //  gq -> qqbq   (1/2 * 1/Na * 1/Nc)
      amp3 -> su3_tree(3,0,1,2,-1, sub+7);
      out[2] += (sub[7] *= -nf1/(Nc*Na)) + (sub[8] /= -2.0*Nc*Na);
      
      //  qr -> qrg    (1/Nc * 1/Nc)
      //  qq -> qqg    (1/2 * 1/Nc * 1/Nc)
      amp3 -> su3_tree(1,-1,2,0,3, sub+9);
      out[3] += (sub[9] /= Nc2);
      out[4] += (sub[10] /= 2.0*Nc2);
      
      //  qqb -> rrbg  (1/Nc * 1/Nc * (nf-1))
      //  qqb -> qqbg  (1/Nc * 1/Nc)
      amp3 -> su3_tree(0,-1,1,2,3, sub+11); 
      out[5] += (sub[11] *= nf1/Nc2) + (sub[12] /= Nc2);
      
      //  qrb -> qrbg  (1/Nc * 1/Nc)
      amp3 -> su3_tree(1,-1,0,2,3, tmp);
      out[6] += (sub[13] = tmp[0]/Nc2);
    }
  }
  

  void _hhc_jet_base::
  amp_1loop(ampg5 *amp1, ampq2g3 *amp2, ampq4g1 *amp3, double *out)
  {
    double sub[2], nf = Nf, nf1 = Nf-1;
    
    //----- FIVE GLUONS PROCESS -----
    //  gg -> ggg    (1/6 * 1/Na * 1/Na)
    out[0] = (amp1 -> su3_1loop(Nf, -1,0,1,2,3))/(6.0*Na2);
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    //  gg -> qqbg   (1/Na * 1/Na * nf)
    out[0] += nf*(amp2 -> su3_1loop(Nf, 1,2,3,-1,0))/Na2;
    
    //  qg -> qgg    (1/2 * 1/Nc * 1/Na)
    out[1] = (amp2 -> su3_1loop(Nf, 1,-1,2,3,0))/(2.0*Nc*Na);
    
    //  gq -> qgg    (1/2 * 1/Na * 1/Nc)
    out[2] = (amp2 -> su3_1loop(Nf, 1,0,2,3,-1))/(2.0*Na*Nc);
      
    //  qqb -> ggg   (1/6 * 1/Nc * 1/Nc)
    out[5] = (amp2 -> su3_1loop(Nf, 0,-1,1,2,3))/(6.0*Nc2);
    
    //----- FOUR QUARKS ONE GLUON PROCESSES -----
    //  qg -> rrbq   (1/Nc * 1/Na * (nf-1))
    //  qg -> qqbq   (1/2 * 1/Nc * 1/Na)
    amp3 -> su3_1loop(Nf, 3,-1,1,2,0, sub);
    out[1] += nf1*sub[0]/(Nc*Na) + sub[1]/(2.0*Nc*Na);
    
    //  gq -> rrbq   (1/Na * 1/Nc * (nf-1))
    //  gq -> qqbq   (1/2 * 1/Na * 1/Nc)
    amp3 -> su3_1loop(Nf, 3,0,1,2,-1, sub);
    out[2] += nf1*sub[0]/(Nc*Na) + sub[1]/(2.0*Nc*Na);
      
    //  qr -> qrg    (1/Nc * 1/Nc)
    //  qq -> qqg    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_1loop(Nf, 1,-1,2,0,3, sub);
    out[3] = sub[0]/Nc2;
    out[4] = sub[1]/(2.0*Nc2);
    
    //  qqb -> rrbg  (1/Nc * 1/Nc * (nf-1))
    //  qqb -> qqbg  (1/Nc * 1/Nc)
    amp3 -> su3_1loop(Nf, 0,-1,1,2,3, sub); 
    out[5] += nf1*sub[0]/Nc2 + sub[1]/Nc2;
    
    //  qrb -> qrbg  (1/Nc * 1/Nc)
    amp3 -> su3_1loop(Nf, 1,-1,0,2,3, sub);
    out[6] = sub[0]/Nc2;
  }

  void _hhc_jet_base::
  amp_1loop_mch(ampg5 *amp1, ampq2g3 *amp2, ampq4g1 *amp3, double *out)
  {
    double sub[2], nf = Nf, nf1 = Nf-1;
    
    //----- FIVE GLUONS PROCESS -----
    //  gg -> ggg    (1/6 * 1/Na * 1/Na)
    out[0] = (amp1 -> su3_1loop_mch(Nf, -1,0,1,2,3))/(6.0*Na2);

    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    //  gg -> qqbg   (1/Na * 1/Na * nf)
    out[0] += nf*(amp2 -> su3_1loop_mch(Nf, 1,2,3,-1,0))/Na2;
    
    //  qg -> qgg    (1/2 * 1/Nc * 1/Na)
    out[1] = (amp2 -> su3_1loop_mch(Nf, 1,-1,2,3,0))/(2.0*Nc*Na);
    
    //  gq -> qgg    (1/2 * 1/Na * 1/Nc)
    out[2] = (amp2 -> su3_1loop_mch(Nf, 1,0,2,3,-1))/(2.0*Na*Nc);
      
    //  qqb -> ggg   (1/6 * 1/Nc * 1/Nc)
    out[5] = (amp2 -> su3_1loop_mch(Nf, 0,-1,1,2,3))/(6.0*Nc2);
        
    //----- FOUR QUARKS ONE GLUON PROCESSES -----
    //  qg -> rrbq   (1/Nc * 1/Na * (nf-1))
    //  qg -> qqbq   (1/2 * 1/Nc * 1/Na)
    amp3 -> su3_1loop_mch(Nf, 3,-1,1,2,0, sub);
    out[1] += nf1*sub[0]/(Nc*Na) + sub[1]/(2.0*Nc*Na);
    
    //  gq -> rrbq   (1/Na * 1/Nc * (nf-1))
    //  gq -> qqbq   (1/2 * 1/Na * 1/Nc)
    amp3 -> su3_1loop_mch(Nf, 3,0,1,2,-1, sub);
    out[2] += nf1*sub[0]/(Nc*Na) + sub[1]/(2.0*Nc*Na);
      
    //  qr -> qrg    (1/Nc * 1/Nc)
    //  qq -> qqg    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_1loop_mch(Nf, 1,-1,2,0,3, sub);
    out[3] = sub[0]/Nc2;
    out[4] = sub[1]/(2.0*Nc2);
    
    //  qqb -> rrbg  (1/Nc * 1/Nc * (nf-1))
    //  qqb -> qqbg  (1/Nc * 1/Nc)
    amp3 -> su3_1loop_mch(Nf, 0,-1,1,2,3, sub); 
    out[5] += nf1*sub[0]/Nc2 + sub[1]/Nc2;
    
    //  qrb -> qrbg  (1/Nc * 1/Nc)
    amp3 -> su3_1loop_mch(Nf, 1,-1,0,2,3, sub);
    out[6] = sub[0]/Nc2;
  }
  

  void _hhc_jet_base::
  amp_ccgg(ampg5 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[0] = amp -> su3_cc(i,j, -1,0,p1,p2,p3);
    cc[0].first /= Na2;
    cc[0].second /= Na2;
  }

  void _hhc_jet_base::
  amp_ccgg(ampq2g3 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[0] = amp -> su3_cc(i,j, p1,p2,-1, 0, p3);
    cc[0].first /= Na2;
    cc[0].second /= Na2;
  }
 
  void _hhc_jet_base::
  amp_ccag(ampq2g3 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[1] = amp -> su3_cc(i,j, -1, p1, 0, p2, p3);
    cc[1].first /= -Na*Nc;
    cc[1].second /= -Na*Nc;
  }

  void _hhc_jet_base::
  amp_ccqg(ampq2g3 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[1] = amp -> su3_cc(i,j, p1, -1, 0, p2, p3);
    cc[1].first /= -Na*Nc;
    cc[1].second /= -Na*Nc;
  }

  void _hhc_jet_base::
  amp_ccga(ampq2g3 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[2] = amp -> su3_cc(i,j, 0, p1, -1, p2, p3);
    cc[2].first /= -Na*Nc;
    cc[2].second /= -Na*Nc;
  }

  void _hhc_jet_base::
  amp_ccgq(ampq2g3 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[2] = amp -> su3_cc(i,j, p1, 0, -1, p2, p3);
    cc[2].first /= -Na*Nc;
    cc[2].second /= -Na*Nc;
  }

  void _hhc_jet_base::
  amp_ccaq(ampq2g3 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[5] = amp -> su3_cc(i,j, -1, 0, p1, p2, p3);
    cc[5].first /= Nc2;
    cc[5].second /= Nc2;
  }

  void _hhc_jet_base::
  amp_ccqa(ampq2g3 *amp, int i, int j, int p1, int p2, int p3, _Pair *cc) 
  {
    cc[5] = amp -> su3_cc(i,j, 0, -1, p1, p2, p3);
    cc[5].first /= Nc2;
    cc[5].second /= Nc2;
  }

  void _hhc_jet_base::
  amp_cc(ampq4g1 *amp, int i, int j, int p1, int p2, int p3, int p4, int p5, _Pair *cc)
  {
    amp -> su3_cc(i,j, p1,p2,p3,p4,p5, cc);
    if(p5 == 0 || p5 == -1) {
      cc[0].first /= -Nc*Na; cc[0].second /= -Nc*Na;
      cc[1].first /= -Nc*Na; cc[1].second /= -Nc*Na;
    } else { 
      cc[0].first /= Nc2; cc[0].second /= Nc2;
      cc[1].first /= Nc2; cc[1].second /= Nc2;
    }
  }
  
  void _hhc_jet_base::
  amp_kp(double al, ampg5 *amp1, ampq2g3 *amp2, ampq4g1 *amp3, su3_kp_i2 *out)
  {
    static su3_kp_i2 tmp[2];
    double nf = Nf, nf1 = Nf-1;
    
    //----- FIVE GLUONS PROCESS -----
    //  gg -> ggg    (1/6 * 1/Na * 1/Na)
    amp1 -> su3_kp(Nf, -1,0, -1,0,1,2,3, tmp, al);
    out[0] = tmp[0]/(6.0*Na2);
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    //  gg -> qqbg   (1/Na * 1/Na * nf)
    amp2 -> su3_kp(Nf, -1,0, 1,2,3,-1,0, tmp, al);
    out[0] += tmp[0]*(nf/Na2);
      
    //  qg -> qgg    (1/2 * 1/Nc * 1/Na)
    amp2 -> su3_kp(Nf, -1,0, 1,-1,2,3,0, tmp, al);
    out[1] = tmp[0]/(-2.0*Nc*Na);
    
    //  gq -> qgg    (1/2 * 1/Na * 1/Nc)
    amp2 -> su3_kp(Nf, -1,0, 1,0,2,3,-1, tmp, al);
    out[2] = tmp[0]/(-2.0*Na*Nc);
    
    //  qqb -> ggg   (1/6 * 1/Nc * 1/Nc)
    amp2 -> su3_kp(Nf, -1,0, 0,-1,1,2,3, tmp, al);
    out[5] = tmp[0]/(6.0*Nc2);
      
    //----- FOUR QUARKS ONE GLUON PROCESSES -----
    //  qg -> rrbq   (1/Nc * 1/Na * (nf-1))
    //  qg -> qqbq   (1/2 * 1/Nc * 1/Na)
    amp3 -> su3_kp(Nf, -1,0, 3,-1,1,2,0, tmp, al);
    out[1] -= tmp[0]*nf1/(Nc*Na) + tmp[1]/(2.0*Nc*Na);
    
    //  gq -> rrbq   (1/Na * 1/Nc * (nf-1))
    //  gq -> qqbq   (1/2 * 1/Na * 1/Nc)
    amp3 -> su3_kp(Nf, -1,0, 3,0,1,2,-1, tmp, al);
    out[2] -= tmp[0]*nf1/(Nc*Na) + tmp[1]/(2.0*Nc*Na);
    
    //  qr -> qrg    (1/Nc * 1/Nc)
    //  qq -> qqg    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_kp(Nf, -1,0, 1,-1,2,0,3, tmp, al);
    out[3] = tmp[0]/Nc2;
    out[4] = tmp[1]/(2.0*Nc2);
    
    //  qqb -> rrbg  (1/Nc * 1/Nc * (nf-1))
    //  qqb -> qqbg  (1/Nc * 1/Nc)
    amp3 -> su3_kp(Nf, -1,0, 0,-1,1,2,3, tmp, al); 
    out[5] += tmp[0]*(nf1/Nc2) + tmp[1]/Nc2;
    
    //  qrb -> qrbg  (1/Nc * 1/Nc)
    amp3 -> su3_kp(Nf, -1,0, 1,-1,0,2,3, tmp, al);
    out[6] = tmp[0]/Nc2;
  }

  // 0.  gg  -> gggg    (1/24 * 1/Na * 1/Na)
  // 1.  gg  -> qqbgg   (1/2 * 1/Na * 1/Na * nf)
  // 2.  qg  -> qggg    (1/6 * 1/Nc * 1/Na)
  // 3.  gq  -> qggg    (1/6 * 1/Na * 1/Nc)
  // 4.  qqb -> gggg    (1/24 * 1/Nc * 1/Nc)
  // 5.  gg  -> qqbrrb  (1/Na * 1/Na * nf * (nf-1) * 1/2)
  // 6.  gg  -> qqbqqb  (1/2 * 1/2 * 1/Na * 1/Na * nf)
  // 7.  qg  -> rrbqg   (1/Nc * 1/Na * (nf-1))
  // 8.  qg  -> qqbqg   (1/2 * 1/Nc * 1/Na)
  // 9.  gq  -> rrbqg   (1/Na * 1/Nc * (nf-1))
  // 10. gq  -> qqbqg   (1/2 * 1/Na * 1/Nc)
  // 11. qr  -> qrgg    (1/2 * 1/Nc * 1/Nc)
  // 12. qq  -> qqgg    (1/2 *1/2 * 1/Nc * 1/Nc)
  // 13. qqb -> rrbgg   (1/2 * 1/Nc * 1/Nc * (nf-1))
  // 14. qqb -> qqbgg   (1/2 * 1/Nc * 1/Nc)
  // 15. qrb -> qrbgg   (1/2 * 1/Nc * 1/Nc)
  // 16. qqb -> rrbssb  (1/Nc * 1/Nc * (nf-1) * (nf-2) * 1/2)
  // 17. qqb -> qqbrrb  (1/Nc * 1/Nc * (nf-1))
  // 18. qqb -> rrbrrb  (1/2 * 1/2 * 1/Nc * 1/Nc * (nf-1))
  // 19. qqb -> qqbqqb  (1/2 * 1/2 * 1/Nc * 1/Nc)    
  // 20. qr  -> qrssb   (1/Nc * 1/Nc * (nf-2))
  // 21. qr  -> qrqqb   (1/2 * 1/Nc * 1/Nc)
  // 22. qq  -> qqrrb   (1/2 * 1/Nc * 1/Nc * (nf-1))
  // 23. qq  -> qqqqb   (1/6 * 1/Nc * 1/Nc)
  // 24. qrb -> qrbssb  (1/Nc * 1/Nc * (nf-2))
  // 25. qrb -> qrbqqb  (1/2 * 1/Nc * 1/Nc)
  // 26. qrb -> qrbrrb  (1/2 * 1/Nc * 1/Nc)
  void _hhc_jet_base::
  amp_tree(ampg6 *amp1, ampq2g4 *amp2, ampq4g2 *amp3, 
	   ampq6 *amp4, double *subp, double *out)
  {
    static double __sub[27];
    double *sub = (subp ? subp : __sub);
    double tmp[5], nf = Nf, nf1 = Nf-1, nf2 = Nf-2;
    std::memset(out, 0, 7*sizeof(double));

    //----- SIX GLUONS PROCESS -----
    if(amp1) {
      //  gg -> gggg    (1/24 * 1/Na * 1/Na)
      out[0] += (sub[0] = (amp1 -> su3_tree(-1,0,1,2,3,4))/(24.0*Na2));
    }
    
    //----- TWO QUARKS FOUR GLUONS PROCESSES -----
    if(amp2) {
      //  gg -> qqbgg   (1/2 * 1/Na * 1/Na * nf)
      out[0] += (sub[1] = nf*(amp2 -> su3_tree(1,2,3,4,-1,0))/(2.0*Na2));
      
      //  qg -> qggg    (1/6 * 1/Nc * 1/Na)
      out[1] += (sub[2] = (amp2 -> su3_tree(1,-1,2,3,4,0))/(6.0*Nc*Na));
      
      //  gq -> qggg    (1/6 * 1/Na * 1/Nc)
      out[2] += (sub[3] = (amp2 -> su3_tree(1,0,2,3,4,-1))/(6.0*Na*Nc));
      
      //  qqb -> gggg   (1/24 * 1/Nc * 1/Nc)
      out[5] += (sub[4] = (amp2 -> su3_tree(0,-1,1,2,3,4))/(24.0*Nc2));
    }
    
    //----- FOUR QUARKS TWO GLUON PROCESSES -----
    if(amp3) {
      //  gg -> qqbrrb  (1/Na * 1/Na * nf * (nf-1) * 1/2)
      //  gg -> qqbqqb  (1/2 * 1/2 * 1/Na * 1/Na * nf)
      amp3 -> su3_tree(1,2,3,4,0,-1, sub+5);
      out[0] += (sub[5] *= nf*nf1/(2.0*Na2)) + (sub[6] *= nf/(4.0*Na2)); 
      
      //  qg -> rrbqg   (1/Nc * 1/Na * (nf-1))
      //  qg -> qqbqg   (1/2 * 1/Nc * 1/Na)
      amp3 -> su3_tree(3,-1,1,2,0,4, sub+7);
      out[1] += (sub[7] *= nf1/(Nc*Na)) +  (sub[8] /= 2.0*Nc*Na);
      
      //  gq -> rrbqg   (1/Na * 1/Nc * (nf-1))
      //  gq -> qqbqg   (1/2 * 1/Na * 1/Nc)
      amp3 -> su3_tree(3,0,1,2,-1,4, sub+9);
      out[2] += (sub[9] *= nf1/(Nc*Na)) + (sub[10] /= 2.0*Nc*Na);
      
      //  qr -> qrgg    (1/2 * 1/Nc * 1/Nc)
      //  qq -> qqgg    (1/2 *1/2 * 1/Nc * 1/Nc)
      amp3 -> su3_tree(1,-1,2,0,3,4, sub+11);
      out[3] += (sub[11] /= 2.0*Nc2);
      out[4] += (sub[12] /= 4.0*Nc2);
      
      //  qqb -> rrbgg  (1/2 * 1/Nc * 1/Nc * (nf-1))
      //  qqb -> qqbgg  (1/2 * 1/Nc * 1/Nc)
      amp3 -> su3_tree(0,-1,1,2,3,4, sub+13); 
      out[5] += (sub[13] *= nf1/(2.0*Nc2)) + (sub[14] /= 2.0*Nc2);
      
      //  qrb -> qrbgg  (1/2 * 1/Nc * 1/Nc)
      amp3 -> su3_tree(1,-1,0,2,3,4, sub+15);
      out[6] += (sub[15] /= 2.0*Nc2);
    }
    
    //----- SIX QUARKS PROCESS -----
    if(amp4) {
      //  qqb -> rrbssb  (1/Nc * 1/Nc * (nf-1) * (nf-2) * 1/2)
      //  qqb -> qqbrrb  (1/Nc * 1/Nc * (nf-1))
      //  qqb -> rrbrrb  (1/2 * 1/2 * 1/Nc * 1/Nc * (nf-1))
      //  qqb -> qqbqqb  (1/2 * 1/2 * 1/Nc * 1/Nc)    
      amp4 -> su3_tree(0,-1,1,2,3,4, "11011", tmp); 
      out[5] += (sub[16] = nf1*nf2*tmp[0]/(2.0*Nc2));
      out[5] += (sub[17] = nf1*tmp[1]/Nc2);
      out[5] += (sub[18] = nf1*tmp[3]/(4.0*Nc2)); 
      out[5] += (sub[19] = tmp[4]/(4.0*Nc2));
      
      //  qr -> qrssb    (1/Nc * 1/Nc * (nf-2))
      //  qr -> qrqqb    (1/2 * 1/Nc * 1/Nc)
      //  qq -> qqrrb    (1/2 * 1/Nc * 1/Nc * (nf-1))
      //  qq -> qqqqb    (1/6 * 1/Nc * 1/Nc)
      amp4 -> su3_tree(1,-1,2,0,3,4, "11101", tmp);
      out[3] += (sub[20] = nf2*tmp[0]/Nc2) + (sub[21] = tmp[2]/(2.0*Nc2));
      out[4] += (sub[22] = nf1*tmp[1]/(2.0*Nc2)) + (sub[23] = tmp[4]/(6.0*Nc2));
      
      //  qrb -> qrbssb  (1/Nc * 1/Nc * (nf-2))
      //  qrb -> qrbqqb  (1/2 * 1/Nc * 1/Nc)
      //  qrb -> qrbrrb  (1/2 * 1/Nc * 1/Nc)
      amp4 -> su3_tree(1,-1,0,2,3,4, "10110", tmp);
      out[6] += (sub[24] = nf2*tmp[0]/Nc2); 
      out[6] += (sub[25] = tmp[2]/(2.0*Nc2));
      out[6] += (sub[26] = tmp[3]/(2.0*Nc2));
    }
  }



#define A(i) kp[i].tree
#define C(i) kp[i].cca
#define P(i) kp[i].pa
#define G(i) kp[i].ga

  void _hhc_jet_base::
  __conv_x1(double eta, double x, double xjac, double al, const su3_kp_i2 *kp, weight_hhc *S) 
  {
    unsigned int nf1 = Nf-1;
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);
        
    //----- K term -----
    double k[4][2];
    Kgg(x, xjac, Nf, al, k[0]); Kgq(x, xjac, al, k[1]); 
    Kqg(x, xjac, al, k[2]); Kqq(x, xjac, al, k[3]); 
    k[0][1] += Ca*lie; k[3][1] += Cf*lie;
    
    S[0][0] = k[0][0]*A(0) + k[1][0]*2.0*Nf*A(1);
    S[0][1] = k[3][0]*A(1) + k[2][0]*A(0);
    S[0][2] = k[0][0]*A(2) + k[1][0]*(A(4)+A(5)+nf1*(A(3)+A(6)));
    S[0][3] = k[3][0]*A(3) + k[2][0]*A(2);
    S[0][4] = k[3][0]*A(4) + k[2][0]*A(2);
    S[0][5] = k[3][0]*A(5) + k[2][0]*A(2);
    S[0][6] = k[3][0]*A(6) + k[2][0]*A(2);
    
    //----- Ktilde term -----
    double t[4][2];
    Tgg(x, xjac, al, t[0]); Tgq(x, xjac, al, t[1]); 
    Tqg(x, xjac, al, t[2]); Tqq(x, xjac, al, t[3]); 
    t[0][1] += Ca*lne*lne; t[3][1] += Cf*lne*lne; 
    
    S[0][0] += t[0][0]*C(0) + t[1][0]*2.0*Nf*C(1);
    S[0][1] += t[3][0]*C(1) + t[2][0]*C(0);
    S[0][2] += t[0][0]*C(2) + t[1][0]*(C(4)+C(5)+nf1*(C(3)+C(6)));
    S[0][3] += t[3][0]*C(3) + t[2][0]*C(2);
    S[0][4] += t[3][0]*C(4) + t[2][0]*C(2);
    S[0][5] += t[3][0]*C(5) + t[2][0]*C(2);
    S[0][6] += t[3][0]*C(6) + t[2][0]*C(2);
    
    //----- P term -----
    double p[4][2];
    Pgg(x, xjac, Nf, p[0]); Pgq(x, xjac, p[1]); 
    Pqg(x, xjac, p[2]); Pqq(x, xjac, p[3]); 
    p[0][1] += 2.0*Ca*lne; p[3][1] += 2.0*Cf*lne; 

    S[0][0] += p[0][0]*P(0) + p[1][0]*2.0*Nf*P(1);
    S[0][1] += p[3][0]*P(1) + p[2][0]*P(0);
    S[0][2] += p[0][0]*P(2) + p[1][0]*(P(4)+P(5)+nf1*(P(3)+P(6)));
    S[0][3] += p[3][0]*P(3) + p[2][0]*P(2);
    S[0][4] += p[3][0]*P(4) + p[2][0]*P(2);
    S[0][5] += p[3][0]*P(5) + p[2][0]*P(2);
    S[0][6] += p[3][0]*P(6) + p[2][0]*P(2);
    
    //----- G term -----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[0][0] += g[0]*G(0); 
    S[0][1] += g[0]*G(1); 
    S[0][2] += g[0]*G(2); 
    S[0][3] += g[0]*G(3); 
    S[0][4] += g[0]*G(4); 
    S[0][5] += g[0]*G(5); 
    S[0][6] += g[0]*G(6);

    //----- dirac delta and + description terms -----
    S[2][0] = g[1]*G(0) + k[0][1]*A(0) + t[0][1]*C(0) + p[0][1]*P(0);
    S[2][1] = g[1]*G(1) + k[3][1]*A(1) + t[3][1]*C(1) + p[3][1]*P(1);
    S[2][2] = g[1]*G(2) + k[0][1]*A(2) + t[0][1]*C(2) + p[0][1]*P(2);
    S[2][3] = g[1]*G(3) + k[3][1]*A(3) + t[3][1]*C(3) + p[3][1]*P(3);
    S[2][4] = g[1]*G(4) + k[3][1]*A(4) + t[3][1]*C(4) + p[3][1]*P(4);
    S[2][5] = g[1]*G(5) + k[3][1]*A(5) + t[3][1]*C(5) + p[3][1]*P(5);
    S[2][6] = g[1]*G(6) + k[3][1]*A(6) + t[3][1]*C(6) + p[3][1]*P(6);

    //---- factorization scale dependent term ----
    S[3][0] = -p[0][0]*A(0) - p[1][0]*2.0*Nf*A(1);
    S[3][1] = -p[3][0]*A(1) - p[2][0]*A(0);
    S[3][2] = -p[0][0]*A(2) - p[1][0]*(A(4)+A(5)+nf1*(A(3)+A(6)));
    S[3][3] = -p[3][0]*A(3) - p[2][0]*A(2);
    S[3][4] = -p[3][0]*A(4) - p[2][0]*A(2);
    S[3][5] = -p[3][0]*A(5) - p[2][0]*A(2);
    S[3][6] = -p[3][0]*A(6) - p[2][0]*A(2);

    S[5][0] = -p[0][1]*A(0);
    S[5][1] = -p[3][1]*A(1);
    S[5][2] = -p[0][1]*A(2);
    S[5][3] = -p[3][1]*A(3);
    S[5][4] = -p[3][1]*A(4);
    S[5][5] = -p[3][1]*A(5);
    S[5][6] = -p[3][1]*A(6);
  }

#undef C
#undef P
#undef G

#define C(i) kp[i].ccb
#define P(i) kp[i].pb
#define G(i) kp[i].gb


  void _hhc_jet_base::
  __conv_x2(double eta, double x, double xjac, double al, const su3_kp_i2 *kp, weight_hhc *S) 
  {
    unsigned int nf1 = Nf-1;
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);
        
    //----- K term -----
    double k[4][2];
    Kgg(x, xjac, Nf, al, k[0]); Kgq(x, xjac, al, k[1]); 
    Kqg(x, xjac, al, k[2]); Kqq(x, xjac, al, k[3]); 
    k[0][1] += Ca*lie; k[3][1] += Cf*lie;
    
    S[1][0] = k[0][0]*A(0) + k[1][0]*2.0*Nf*A(2);
    S[1][1] = k[0][0]*A(1) + k[1][0]*(A(4)+A(5)+nf1*(A(3)+A(6)));
    S[1][2] = k[3][0]*A(2) + k[2][0]*A(0);
    S[1][3] = k[3][0]*A(3) + k[2][0]*A(1);
    S[1][4] = k[3][0]*A(4) + k[2][0]*A(1);
    S[1][5] = k[3][0]*A(5) + k[2][0]*A(1);
    S[1][6] = k[3][0]*A(6) + k[2][0]*A(1);
    
    //----- Ktilde term -----
    double t[4][2];
    Tgg(x, xjac, al, t[0]); Tgq(x, xjac, al, t[1]); 
    Tqg(x, xjac, al, t[2]); Tqq(x, xjac, al, t[3]); 
    t[0][1] += Ca*lne*lne; t[3][1] += Cf*lne*lne; 
    
    S[1][0] += t[0][0]*C(0) + t[1][0]*2.0*Nf*C(2);
    S[1][1] += t[0][0]*C(1) + t[1][0]*(C(4)+C(5)+nf1*(C(3)+C(6)));
    S[1][2] += t[3][0]*C(2) + t[2][0]*C(0);
    S[1][3] += t[3][0]*C(3) + t[2][0]*C(1);
    S[1][4] += t[3][0]*C(4) + t[2][0]*C(1);
    S[1][5] += t[3][0]*C(5) + t[2][0]*C(1);
    S[1][6] += t[3][0]*C(6) + t[2][0]*C(1);
    
    //----- P term -----
    double p[4][2];
    Pgg(x, xjac, Nf, p[0]); Pgq(x, xjac, p[1]); 
    Pqg(x, xjac, p[2]); Pqq(x, xjac, p[3]); 
    p[0][1] += 2.0*Ca*lne; p[3][1] += 2.0*Cf*lne; 

    S[1][0] += p[0][0]*P(0) + p[1][0]*2.0*Nf*P(2);
    S[1][1] += p[0][0]*P(1) + p[1][0]*(P(4)+P(5)+nf1*(P(3)+P(6)));
    S[1][2] += p[3][0]*P(2) + p[2][0]*P(0);
    S[1][3] += p[3][0]*P(3) + p[2][0]*P(1);
    S[1][4] += p[3][0]*P(4) + p[2][0]*P(1);
    S[1][5] += p[3][0]*P(5) + p[2][0]*P(1);
    S[1][6] += p[3][0]*P(6) + p[2][0]*P(1);
    
    //----- G term -----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[1][0] += g[0]*G(0); 
    S[1][1] += g[0]*G(1); 
    S[1][2] += g[0]*G(2); 
    S[1][3] += g[0]*G(3); 
    S[1][4] += g[0]*G(4); 
    S[1][5] += g[0]*G(5); 
    S[1][6] += g[0]*G(6);

    //----- dirac delta and + description terms -----
    S[2][0] += g[1]*G(0) + k[0][1]*A(0) + t[0][1]*C(0) + p[0][1]*P(0);
    S[2][1] += g[1]*G(1) + k[0][1]*A(1) + t[0][1]*C(1) + p[0][1]*P(1);
    S[2][2] += g[1]*G(2) + k[3][1]*A(2) + t[3][1]*C(2) + p[3][1]*P(2);
    S[2][3] += g[1]*G(3) + k[3][1]*A(3) + t[3][1]*C(3) + p[3][1]*P(3);
    S[2][4] += g[1]*G(4) + k[3][1]*A(4) + t[3][1]*C(4) + p[3][1]*P(4);
    S[2][5] += g[1]*G(5) + k[3][1]*A(5) + t[3][1]*C(5) + p[3][1]*P(5);
    S[2][6] += g[1]*G(6) + k[3][1]*A(6) + t[3][1]*C(6) + p[3][1]*P(6);

    //---- factorization scale dependent term ----
    S[4][0] = -p[0][0]*A(0) - p[1][0]*2.0*Nf*A(2);
    S[4][1] = -p[0][0]*A(1) - p[1][0]*(A(4)+A(5)+nf1*(A(3)+A(6)));
    S[4][2] = -p[3][0]*A(2) - p[2][0]*A(0);
    S[4][3] = -p[3][0]*A(3) - p[2][0]*A(1);
    S[4][4] = -p[3][0]*A(4) - p[2][0]*A(1);
    S[4][5] = -p[3][0]*A(5) - p[2][0]*A(1);
    S[4][6] = -p[3][0]*A(6) - p[2][0]*A(1);

    S[5][0] += -p[0][1]*A(0);
    S[5][1] += -p[0][1]*A(1);
    S[5][2] += -p[3][1]*A(2);
    S[5][3] += -p[3][1]*A(3);
    S[5][4] += -p[3][1]*A(4);
    S[5][5] += -p[3][1]*A(5);
    S[5][6] += -p[3][1]*A(6);
  }
}
