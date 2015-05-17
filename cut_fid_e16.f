      logical function fidu_e_sub(pel,thetael,phiel,in_fid_reg)
      REAL phiel_s,del_phie
      REAL THETACUT
      REAL EXPON,PI
      REAL pel,thetael,phiel
      logical in_fid_reg 
      integer mod6
      REAL pshift
      DATA pshift/0.14/
      REAL c1,c2,c3,c4
      real factor
      DATA c1,c2,c3,c4/12.0,18.5,0.25,25.0/
      DATA PI/3.141592/
      data factor/0.416667/


      phiel_s = phiel+180.
      mod6=(phiel_s+30.)/60
      phiel_s = phiel_s - mod6*60.

      thetacut=c1+c2/((pel+pshift))
      expon=c3*(pel**factor)

      del_phie = c4*sin((thetael-thetacut)*pi/180.)**expon 
      if (abs(phiel_s).le.del_phie.and.thetael.ge.thetacut) then
      in_fid_reg=.TRUE.
      else
      in_fid_reg=.FALSE.
      endif

      !write(*,*) 'in_fid_reg=',in_fid_reg
      fidu_e_sub=in_fid_reg
      END


       logical function hadronfid(result,theta,phi,s)
       implicit none
       real theta,phi,phic
       real a0p(6),a0m(6),a1p(6),a1m(6),a2p(6),a2m(6),a3p(6),a3m(6)
       logical result
       integer s
       real phi_min,phi_max
       data a0p/24.,24.,23.,23.5,24.5,24.5/
       data a0m/25.,26.,26.,25.5,27.,26./
       data a1p/0.22,0.23,0.2,0.2,0.22,0.22/
       data a1m/0.22,0.22,0.22,0.22,0.16,0.16/
       data a2p/8.,8.,8.,8.,8.,8./
       data a2m/8.,8.,8.,8.,8.,8./
       data a3p/1.,1.,1.,1.,1.,1./
       data a3m/1.,1.,1.,1.,1.,1./

       phic=phi + 180.
       phic = phic - INT((phic+30.)/60)*60
       phi_min=-a0m(s)*(1-exp(-a1m(s)*(theta-a2m(s))))+a3m(s)
       phi_max= a0p(s)*(1-exp(-a1p(s)*(theta-a2p(s))))+a3p(s)

       if (phic.ge.phi_min.and.phic.le.phi_max) then
       result=.TRUE.
       else 
       result=.FALSE.
       endif

       !write(*,*) 'result=',result
       hadronfid=result 
       end
