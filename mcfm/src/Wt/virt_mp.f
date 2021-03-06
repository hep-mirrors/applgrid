      double complex function virt_mp(mQ,ig,is,ie,in,ic,p)
      implicit none
************************************************************************
*     Author: F. Tramontano                                            *
*     October, 2004                                                    *
************************************************************************
c---- One-loop -+ helicity amplitude for W+c production
c---- hQ=-1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
c---- hg=+1  with s(mu) as gauge vector
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + Qbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ Q(p5) 
c----
      include 'constants.f'
      include 'epinv.f'
      include 'b0.f'
      include 'scheme.f'
      include 'scale.f'
      include 'zprods_decl.f'
      integer is,ig,ic,ie,in,nu,jpart
      double precision p(mxpart,4),q(mxpart,4)
      double precision dot,taucg,taugs,taucs,msq,mQ
      double precision qsq,tcg,tcs,qsqhat,qcg,qcs,cgs1
      double precision ddilog,epin,epin2,xlog
      double complex smp
      double complex amp1,amp2,amp3,amp4
      double complex ffcg1,ffcg2,ffcg3,ffcs1,ffcs2,fun1,fun2,fun3,fun4
      double complex fL6m1,fL6m2,fL6m3,L6m1,L6m2,L6m3
      double complex lnrat,C0fa2m,C0fb2m,I3me
c      double precision mass2,width2,mass3,width3
c      common/breit/n2,n3,mass2,width2,mass3,width3
      scheme='dred'
      msq=mQ**2
      xlog=log(musq/msq)
      epin=epinv+xlog
      epin2=epinv**2+epinv*xlog+half*xlog**2
      taugs=+2d0*dot(p,is,ig)
      taucs=+2d0*dot(p,ic,is)
      taucg=+2d0*dot(p,ic,ig)
      qsqhat=taucs+taucg+taugs
      qsq=qsqhat+msq
      tcg=taucg+msq
      tcs=taucs+msq
      qcg=qsqhat-taucg
      qcs=qsqhat-taucs
      cgs1=taucs-msq*taugs/taucg
      ffcg2=fun3(tcg,qsq,msq)
      ffcs1=fun4(tcs,qsq,msq)
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
      ffcs2=fun2(tcs,qsq,msq)
      fL6m1=L6m1(taucg,taucs,taugs,msq,qsq)
      fL6m2=L6m2(taucg,taucs,taugs,msq,qsq)
      fL6m3=L6m3(taucg,taucs,taugs,msq,qsq)
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart.eq.ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      amp1=za(ig,ie)*zb(ig,in)/za(ig,is)**2/zb(is,ic)
      amp2=za(is,ie)*zb(is,in)/za(ig,is)**2/zb(is,ic)
      amp3=za(ig,ie)*za(is,ic)*zb(is,in)/za(ig,is)**2
     . /za(ig,ic)/zb(is,ic)
      amp4=za(ie,ic)*zb(ig,in)/za(ig,is) + za(is,ic)*
     . za(ie,ic)*zb(is,in)/za(ig,is)/za(ig,ic)
      smp=  + taucs*taucg**(-1)*taugs*xn**(-1)*amp1 - 1.D0/2.D0*taucs*
     &    taucg*taugs**(-1)*xn**(-1)*amp3*fL6m2 + 1.D0/2.D0*taucs*
     &    taugs**(-1)*xn**(-1)*cgs1*amp1*fL6m2 - 1.D0/2.D0*taucs*
     &    xn**(-1)*amp3*fL6m1 + 1.D0/2.D0*taucs*xn*amp3*fL6m3 - 
     &    taucg**(-2)*taugs**2*tcg*cf*amp2 - 2.D0*taucg**(-2)*taugs**2*
     &    mQ**2*cf*amp2*ffcg3 - 2.D0*taucg**(-1)*taugs*tcg*cf*amp1 + 1.D
     &    0/2.D0*taucg**(-1)*taugs*tcg*xn**(-1)*amp1*fL6m1 - 1.D0/2.D0*
     &    taucg**(-1)*taugs*tcg*xn*amp1*fL6m3 + taucg**(-1)*taugs*
     &    qsqhat*cf*amp2 + 2.D0*taucg**(-1)*taugs*mQ**2*cf*amp2 - 2.D0*
     &    taucg**(-1)*taugs*mQ**2*cf*amp3*ffcg3 - 1.D0/2.D0*taucg**(-1)
     &    *taugs*mQ**2*xn**(-1)*amp2*fL6m2 + taucg**(-1)*taugs*mQ**2*xn
     &    *amp2*ffcg3 + taucg**(-1)*taugs*mQ**2*xn*amp2*fL6m3 - 
     &    taucg**(-1)*taugs*mQ**2*xn*amp3*fL6m3 - taucg**(-1)*taugs**2*
     &    tcs*xn**(-1)*amp1*ffcs1 - taucg**(-1)*taugs**2*qsq**(-1)*
     &    mQ**2*cf*amp2*ffcg1 - taucg**(-1)*taugs**2*cf*amp2*ffcg2 - 2.D
     &    0*taucg**(-1)*taugs**2*xn**(-1)*amp1*ffcs2
      smp = smp + taucg**(-1)*qsqhat*cf*cgs1*amp1 - 1.D0/2.D0*taucg*
     &    xn**(-1)*amp3*fL6m1 + 1.D0/2.D0*taucg*xn*amp3*fL6m3 - 3.D0*
     &    taugs*tcg*qsq**(-1)*cf*amp3*ffcg3 - taugs*tcg*qsq**(-1)*
     &    xn**(-1)*amp2*ffcg3 + 2.D0*taugs*tcg*qsq**(-1)*xn*amp2*ffcg3
     &     - 2.D0*taugs*cf*amp3 + taugs*xn**(-1)*amp3*ffcs2 - taugs*
     &    xn**(-1)*amp3 - taugs**2*qsq**(-1)*cf*amp2*ffcg3 + 1.D0/2.D0*
     &    tcs*xn**(-1)*amp3*fL6m2 - 2.D0*qsqhat*cf*amp2 + 1.D0/2.D0*
     &    qsqhat*xn**(-1)*amp2*fL6m1 - 1.D0/2.D0*qsqhat*xn*amp2*fL6m3
     &     + cf*qcs*amp3 + cf*cgs1*amp3 - 3.D0/2.D0*cf*amp4*epin - 5.D0/
     &    2.D0*cf*amp4 - b0*amp4*epin - xn**(-1)*cgs1*amp3*fL6m2 + 1.D0/
     &    2.D0*xn**(-1)*pisqo6*amp4 + 1.D0/2.D0*xn**(-1)*amp4*epin + 1.D
     &    0/2.D0*xn**(-1)*amp4*epin2 - 1.D0/2.D0*xn*pisqo6*amp4 - 1.D0/
     &    2.D0*xn*amp4*epin - 3.D0/2.D0*xn*amp4*epin2 + ddilog(
     &    msq**(-1)*tcs)*xn**(-1)*amp4 - ddilog(msq**(-1)*tcg)*xn*amp4
     &     - I3me(msq,taugs,qsq)*taucs**(-1)*taucg**(-1)*taugs**2*tcg*
     &    mQ**2*xn**(-1)*amp1
      smp = smp + I3me(msq,taugs,qsq)*taucs**(-1)*taucg*taugs*mQ**2*
     & xn**(-1)*amp3 - I3me(msq,taugs,qsq)*taucs**(-1)*taugs*mQ**2*
     &    xn**(-1)*qcs*amp2 + I3me(msq,taugs,qsq)*taucg**(-2)*taugs**2*
     &    tcg*mQ**2*xn*amp1 + I3me(msq,taugs,qsq)*taucg**(-1)*taugs*
     &    mQ**2*xn*qcs*amp2 + 2.D0*I3me(msq,taugs,qsq)*taucg**(-1)*
     &    taugs*mQ**2*xn*cgs1*amp2 - 2.D0*I3me(msq,taugs,qsq)*
     &    taucg**(-1)*taugs*mQ**2*xn*cgs1*amp3 - I3me(msq,taugs,qsq)*
     &    taugs*mQ**2*xn*amp3 + lnrat( - taucs,msq)*taucs*taucg**(-1)*
     &    taugs*xn**(-1)*amp1 + lnrat( - taucs,msq)*taucs*taucg**(-1)*
     &    xn**(-1)*qcg*amp1 - lnrat( - taucs,msq)*taucs*xn**(-1)*amp3
     &     + lnrat( - taucs,msq)*taucg**(-1)*taugs*tcs**(-1)*qsq*mQ**2*
     &    xn**(-1)*amp2 + lnrat( - taucs,msq)*taucg**(-1)*taugs**2*
     &    tcs**(-1)*mQ**2*xn**(-1)*amp1 - lnrat( - taucs,msq)*taugs*
     &    tcs**(-1)*mQ**2*xn**(-1)*amp3 - lnrat( - taucs,msq)*xn**(-1)*
     &    amp4*epin + lnrat( - taucs,msq)**2*xn**(-1)*amp4 + 2.D0*
     &    lnrat(
     &  - taucg,msq)*taucg**(-1)*taugs*mQ**2*cf*amp1 + 2.D0*lnrat( - 
     &    taucg,msq)*taucg**(-1)*taugs*cf*qcs*amp2 + lnrat( - taucg,msq
     &    )*taucg*taugs*tcg**(-1)*xn*amp2 - lnrat( - taucg,msq)*taucg*
     &    xn**(-1)*amp3 + 2.D0*lnrat( - taucg,msq)*taugs*tcg**(-1)*
     &    mQ**2*cf*amp3 + 3.D0*lnrat( - taucg,msq)*taugs*cf*amp1 - 2.D0
     &    *lnrat( - taucg,msq)*taugs**2*tcg**(-1)*cf*amp2 + lnrat( - 
     &    taucg,msq)*xn**(-1)*qcg*amp1 + lnrat( - taucg,msq)*xn*amp4*
     &    epin - lnrat( - taucg,msq)**2*xn*amp4 + lnrat( - taugs,msq)*
     &    xn*amp4*epin - 1.D0/2.D0*lnrat( - taugs,msq)**2*xn*amp4 - 
     &    lnrat( - qsqhat,msq)*taucg**(-2)*tcg*qsq**(-1)*qsqhat*cf*
     &    qcs**2*amp2 + lnrat( - qsqhat,msq)*taucg**(-1)*taugs**2*
     &    qsq**(-1)*qsqhat*cf*amp2 + lnrat( - qsqhat,msq)*taucg**(-1)*
     &    taugs**2*qsq**(-1)*qsqhat*xn**(-1)*amp1 - lnrat( - qsqhat,msq
     &    )*taucg**(-1)*qsqhat*xn**(-1)*qcg*amp1 + lnrat( - qsqhat,msq)
     &    *taucg*qsq**(-1)*qsqhat*cf*amp2 - 3.D0*lnrat( - qsqhat,msq)*
     &    taucg*qsq**(-1)*qsqhat*cf*amp3
      smp = smp + lnrat( - qsqhat,msq)*taucg*qsq**(-1)*qsqhat*xn**(-1)*
     & amp3 - lnrat( - qsqhat,msq)*taugs*qsq**(-1)*qsqhat*xn*amp2 - 2.D0
     &    *lnrat( - qsqhat,msq)*tcs*qsq**(-1)*qsqhat*cf*amp3 + lnrat(
     &     - qsqhat,msq)*tcs*qsq**(-1)*qsqhat*xn**(-1)*amp3 - lnrat( - 
     &    qsqhat,msq)*qsq**(-1)*qsqhat*cf*cgs1*amp3 + 3.D0*lnrat( - 
     &    qsqhat,msq)*qsqhat*cf*amp2 + xlog*b0*amp4 + 
     &    C0fb2m(tcg,msq)*taucs**(-1)*taucg*mQ**2*xn**(-1)*cgs1*amp3 - 
     &    C0fb2m(tcg,msq)*taucs**(-1)*taugs*mQ**2*xn**(-1)*cgs1*amp2 + 
     &    2.D0*C0fb2m(tcg,msq)*taucs**(-1)*taugs*mQ**2*xn**(-1)*cgs1*
     &    amp3 - C0fb2m(tcg,msq)*mQ**2*xn**(-1)*cgs1*amp1 - C0fa2m(tcs,
     &    qsq,msq)*taucs**(-1)*taucg**(-2)*taugs**2*mQ**4*xn**(-1)*qcs*
     &    amp2 - 2.D0*C0fa2m(tcs,qsq,msq)*taucs**(-1)*taucg**(-1)*taugs
     &    *mQ**2*xn**(-1)*qcs*cgs1*amp3 - C0fa2m(tcs,qsq,msq)*
     &    taucs**(-1)*mQ**2*xn**(-1)*qcs*cgs1*amp3 - C0fa2m(tcs,qsq,msq
     &    )*taucg**(-1)*taugs**2*mQ**2*xn**(-1)*amp1 + C0fa2m(tcs,qsq,
     &    msq)*taucg**(-1)*mQ**2*xn**(-1)*qcs*cgs1*amp1
      smp = smp + C0fa2m(tcs,qsq,msq)*taugs*mQ**2*xn**(-1)*amp3

c      lomp=za(ie,ic)
c    . /za(ig,is)/za(ig,ic)*(za(is,ic)*zb(is,in)+za(ig,ic)*zb(ig,in)) 
      virt_mp=smp
      return
      end
