      real function elas(eb,theta)
      
c     Returns the e-p differential cross section
c     as a function of the electron scattering angle theta.

c     eb: beam energy (GeV)
c     theta: scattering angle (deg)

      real*4 e_prime              !Energy of the scattered electron
      real*4 gep                  !Electric form factor of proton
      real*4 gmp                  !Magnetic form factor of proton
      real*4 q2,q1	        !Momentum transfer squared
      real*4 alpha                !Mott scattering constant, = e^2/4pi
      real*4 tau                  !q^2/4M^2
      real*4 eb			!beam energy
      real*4 theta                !electron scattering angle
      real*4 rm_p			!proton mass
      real*4 s2,c2,theta2

      data rm_p/.938/

      theta2 = theta*3.1415/180./2.
      s2 = sin(theta2)
      c2 = cos(theta2)
      
c     Energy of the scattered electron from recoil of the target: 
      e_prime=eb/(1.+(2.*eb/rm_p)*s2**2)

c     Q squared, neglecting electron mass:
      q2=4.*eb*e_prime*s2**2
      q1=sqrt(q2)
      
c     Gep and Gmn dipole:
c      gep=1./(1.+q2/0.71)**2
c      gmp=2.793*gep
      
c     Bosted parameterization of form factors
c     Phys. Rev. C 51, 409

      corr1	= 1+0.62*q1+0.68*q2+2.8*q1**3+0.83*q1**4
      corr2	= 1+0.35*q1+2.44*q2+0.5*q1**3+1.04*q1**4+0.34*q1**5
      gep	= 1./corr1
      gmp	= 2.7928/corr2
      
c     Mott
      rmott=389*(1./137.)**2/(4.*eb**2*s2**4)	! microbarns/sr
      
c     Recoil correction 
      rmott=rmott*(e_prime/eb)

c     Form factor
      tau=q2/(4.*rm_p**2)
      elas=rmott*
     $     ( ((gep**2+tau*gmp**2)/(1.+tau))*c2**2 +
     $     2.*tau*gmp**2*s2**2)
     
      return
      end
c------------------------------------------------      
      real function spence(x)
      
      real*8 x
      
      if (abs(x).lt.0.1) then
        spence = x+x**2/4.
      elseif (x.gt.0.99.and.x.lt.1.01) then
        spence = pi*pi/6.
      elseif (x.gt.-1.01.and.x.lt.-0.99) then
        spence = -pi*pi/12.
      elseif (x.gt.0) then
        spence = 0.1025+sintp(x)
      else
        spence = -0.0975+sintn(x)
      endif
      
      end
c---------------------------------------------------      
      real function sintp(x)
      implicit none
      real*8 x,arg
      real*4 xstep,sum,y
      integer i

      xstep	= (x-.1)/100.
      sum	= 0.
      y		= 0.1-xstep/2.
      do i=1,100
        y   	= y+xstep
        arg 	= abs(1.-y)
        sum 	= sum-dlog(arg)/y
      enddo
      sintp	= sum*xstep
      
      end
c---------------------------------------------------
      real function sintn(x)
      implicit none
      real*8 x,y
      real*4 xa,ystep,sum
      integer i

      xa	= abs(x)
      ystep	= (xa-0.1)/100.
      sum	= 0.
      y		= 0.1-ystep/2.
      do i=1,100
        y	= y+ystep
        sum	= sum-dlog(1.+y)/y
      enddo
      sintn	= sum*ystep
      
      end
c--------------------------------------------------      
      real function bfunc(z)
      implicit none
      real*4 z,xi

      xi = log(1440.)-2.*alog(z)/3.
      xi = xi/(alog(183.)-alog(z)/3.)
      bfunc = (4./3.)*(1.+((z+1.)/(z+xi))/(alog(183.)-alog(z)/3.)/9.)
      
      end
c-----------------------------------------------     
      real function elasrad( es, theta_d, t, wcut )  
      
c     Returns radiated cross section using equation II.6 of
c     Mo and Tsai, Rev. Mod. Phys. 41, 205-235 (1969).

c     es: incident electron energy (GeV)
c     theta_d: scattered electron angle (deg)
c     t: target thickness (radiation lengths)
c     wcut: upper limit of W (should be > 0.938+W resolution - e.g. 1.0) (GeV)

c     This program calls function elas defined above.

      real*8 arg
      real*4 me,mp,pi,alpha,z,t,es,theta_d,theta,delta,wcut,spence
      real*4 b,cst1,eta,bt,sigunp,sigel,eel,qs
      real*4 mpi
      real*4 znuc,deltac(27),del_mo,delta_t,idel
      real*4 arg11,arg15,arg19,arg23
      real*4 epr,e1,e3,e4,gamma4,beta4
      real*4 radcor
      real*4 snth
      real*4 bfunc

      data pi/3.1415926/
      data mp/.9382/
      data mpi/.13957/
      data me/0.000511/
      data alpha/7.2993e-3/

c      write(6,*) 'Hello' 
c      write(6,*) 'Received ', es, theta_d, t, wcut
      
      theta	= theta_d*pi/180.

      znuc	= 1.
      z		= 1.
      
      b		= bfunc(z)
      snth	= sin(theta)
      cst1	= 1.-cos(theta)
      eel	= es/(1.+es/mp*cst1)
      epr	= es+mp-eel
      
      epcut	= (es+(mp*mp-wcut*wcut)/2./mp)*eel/es
      delta	= eel-epcut
      
      gamma4	= epr/mp
      beta4	= sqrt(1.-1/gamma4**2)
      e1	= es
      e3	= eel
      e4	= epr
      eta	= es/eel
      qs	= 2.*es*eel*cst1
      
      deltac(1)=28./9.-13./6.*alog(qs/me**2)
      deltac(2)= (alog(qs/me**2) - 1. + 2.*alog(eta) )
     *     * ( 2.*alog(e1/delta) -      3.*alog(eta) )
      arg=(e3-e1)/e3
      deltac(3)=-spence(arg)
      deltac(4)=znuc**2*alog(e4/mp)
      deltac(5)=znuc**2*alog(mp/eta/delta)*
     * (alog((1.+beta4)/(1.-beta4))/beta4-2.)
      deltac(6)=znuc**2*0.5/beta4*(alog((1+beta4)/(1.-beta4))*
     * alog((e4+mp)/2/mp))
      arg=sqrt((e4-mp)*(1.+beta4)/(e4+mp)/(1.-beta4))
      deltac(7)=-znuc**2*spence(-arg)/beta4
      arg=(e3-mp)/e1
      deltac(8)=znuc*spence(arg)
      arg=(mp-e3)*mp/(2.*e3*e4-mp*e1)
      deltac(9)=-znuc*spence(arg)
      arg=2.*e3*(mp-e3)/(2.*e3*e4-mp*e1)
      deltac(10)=znuc*spence(arg)
      arg11=(2.*e3*e4-mp*e1)/e1/(mp-2.*e3)
      arg11=abs(arg11)
      deltac(11)=znuc*alog(arg11)*alog(mp/2./e3)
      arg=(e3-e4)/e3
      deltac(12)=-znuc*spence(arg)
      arg=(e4-e3)*mp/(2.*e1*e4-mp*e3)
      deltac(13)=znuc*spence(arg)
      arg=2.*e1*(e4-e3)/(2*e1*e4-mp*e3)
      deltac(14)=-znuc*spence(arg)
      arg15=(2.*e1*e4-mp*e3)/e3/(mp-2.*e1)
      arg15=abs(arg15)
      deltac(15)=-znuc*alog(arg11)*alog(mp/2./e1)
      arg=(e1-mp)/e1
      deltac(16)=-znuc*spence(arg)
      arg=(mp-e1)/e1
      deltac(17)=znuc*spence(arg)
      arg=2.*(mp-e1)/mp
      deltac(18)=-znuc*spence(arg)
      arg19=abs(mp/(2.*e1-mp))
      deltac(19)=-znuc*alog(arg19)*alog(mp/2./e1)
      arg=(e3-mp)/e3
      deltac(20)=znuc*spence(arg)
      arg=(mp-e3)/e3
      deltac(21)=-znuc*spence(arg)
      arg=2*(mp-e3)/mp
      deltac(22)=znuc*spence(arg)
      arg23=abs(mp/(2.*e3-mp))
      deltac(23)=znuc*alog(arg23)*alog(mp/2./e3)
      arg=(e1-e3)/e1
      deltac(24)=-spence(arg)
      arg=(e4-mp)*(1-beta4)/(e4+mp)/(1+beta4)
      arg=sqrt(arg)
      deltac(25)=znuc**2*spence(arg)/beta4
      arg=(e4-mp)/(e4+mp)
      arg=sqrt(arg)
      deltac(26)=-znuc**2*spence(arg)/beta4
      deltac(27)=znuc**2*spence(-arg)/beta4

      del_mo  = 0.
      do idel=1,2
         del_mo = del_mo + deltac(idel)
      enddo
      del_mo  = -alpha*del_mo/pi

      sigunp  = elas(es,theta_d)

c     Straggling correction delta_t
      
      bt      = b*t
      delta_t = -0.5*bt*(alog(es/eta**2/delta)+alog(eel/delta))
c      print *, qs,del_mo,delta_t
c      radcor  = 1.+del_mo+delta_t
      radcor = exp(del_mo+delta_t)
      elasrad = radcor*sigunp
      
      end

