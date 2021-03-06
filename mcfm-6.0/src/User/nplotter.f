      subroutine nplotter(p,wt,wt2,nd)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c---     nd:  an integer specifying the dipole number of this contribution
c---          (if applicable), otherwise equal to zero
      implicit none
      include 'constants.f'
      include 'process.f'
c   P.S. use of grids
      include 'ptilde.f'
      include 'APPLinclude.f'
c   P.S. end

      double precision p(mxpart,4),wt,wt2
      integer switch,nproc,nd
      common/nproc/nproc

c--- This routine simply picks out a process-specific plotting routine
c---  (if available) and falls back to the generic routine otherwise
c---  so far available: W_only, Z_only, Wbbbar, Wbbmas, WpWp2j, WpWp3j

c--- switch:  an integer equal to either 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      if (nd.gt.0) then 
         switch=1
      else
         switch=0
      endif

      if     (case .eq. 'W_only') then
        call nplotter_W_only(p,wt,wt2,switch)
      elseif (case .eq. 'Z_only') then
        call nplotter_Z_only(p,wt,wt2,switch)
      elseif (case .eq. 'Wbbbar') then
        call nplotter_Wbbmas(p,wt,wt2,switch)
      elseif (case .eq. 'Wbbmas') then
        call nplotter_Wbbmas(p,wt,wt2,switch)
      elseif ((case .eq. 'WpWp2j') .or. (case .eq. 'WpWp3j'))then
        call nplotter_WpWp(p,wt,wt2,switch)
      elseif ((case.eq.'W_1jet') .or. (case.eq.'W_2jet') 
     &   .or. (case.eq.'W_3jet')) then
        call nplotter_Wjets(p,wt,wt2,switch)
c--- photon processes also need to know the dipole number
      elseif ((case.eq.'Wgamma') .or. (case.eq.'Zgamma')
     &   .or. (case.eq.'Wgajet') .or. (case.eq.'Zgajet')) then
        call nplotter_Vgamma(p,wt,wt2,switch,nd)
      elseif ((case.eq.'gamgam') .or. (case .eq. 'gmgmjt'))  then
        call nplotter_gamgam(p,wt,wt2,switch,nd)
      else
        call nplotter_generic(p,wt,wt2,switch)
      endif

c P.S. filling applgrid
      if (creategrid) call fill_grid(p)
c P.S. end

      return
      end
      
