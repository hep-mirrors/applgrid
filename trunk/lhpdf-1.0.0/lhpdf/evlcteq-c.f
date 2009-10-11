      function GetNx()
      implicit double precision (a-h, o-z)
      parameter (mxx = 105)
      common / xxaray / xcr, xmin, xv(0:mxx), lstx, nx
      integer getnx
      getnx = nx
      return
      end

      subroutine GetXv(xxv, xvpow)
      implicit double precision (a-h, o-z)
      parameter (mxx = 105, mxq = 25, mxf = 6)
      parameter (mxpn = mxf * 2 + 2)
      parameter (mxqx= mxq * mxx,   mxpqx = mxqx * mxpn)
      common / xxaray / xcr, xmin, xv(0:mxx), lstx, nx
      common / qaray1 / qini,qmax, qv(0:mxq),tv(0:mxq), nt,jt,ng
      common / evlpac / al, iknl, ipd0, ihdn, nfmx
      common / pevldt / upd(mxpqx), kf, nelmt
      common / comqms / valqms(9)
      double precision xxv(0:*), xvpow(0:*)
      data xpow / 0.3d0 /

      do i = 0,nx
         xxv(i) = xv(i)
         xvpow(i) = xv(i)**xpow
      end do
      return
      end
      
      function GetNt()
      implicit double precision (a-h, o-z)
      parameter (mxq = 25)
      common / qaray1 / qini,qmax, qv(0:mxq),tv(0:mxq), nt,jt,ng
      integer getnt
      getnt = nt
      return
      end

      subroutine GetTv(ttv)
      implicit double precision (a-h, o-z)
      parameter (mxq = 25)
      common / qaray1 / qini,qmax, qv(0:mxq),tv(0:mxq), nt,jt,ng
      double precision ttv(0:*) 
      do i = 0,nt
         ttv(i) = tv(i)
      end do
      return
      end

      function GetAl()
      implicit double precision (a-h, o-z)
      common / evlpac / al, iknl, ipd0, ihdn, nfmx
      double precision getal
      getal = al
      return
      end

      subroutine GetValQms(maxnf, qms)
      implicit double precision (a-h, o-z)
      common / evlpac / al, iknl, ipd0, ihdn, nfmx
      common / comqms / valqms(9)
      double precision qms(*)
      maxnf = nfmx
      do i = 1,6
         qms(i) = valqms(i)
      end do
      return
      end
  
      subroutine GetUpd(uupd)
      implicit double precision (a-h, o-z)
      parameter (mxx = 105, mxq = 25, mxf = 6)
      parameter (mxpn = mxf * 2 + 2)
      parameter (mxqx= mxq * mxx,   mxpqx = mxqx * mxpn)
      common / xxaray / xcr, xmin, xv(0:mxx), lstx, nx
      common / qaray1 / qini,qmax, qv(0:mxq),tv(0:mxq), nt,jt,ng
      common / evlpac / al, iknl, ipd0, ihdn, nfmx
      common / pevldt / upd(mxpqx), kf, nelmt
      double precision uupd(*)
      
      nn = (nx+1)*(nt+1)*(nfmx+3)
      do i = 1,nn
         uupd(i) = upd(i)
      end do
      return
      end
      
