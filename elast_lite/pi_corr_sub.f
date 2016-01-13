      Subroutine pi_corr_sub
     &   (thetapid,phipid,mom_pip,torcur,secth,nthetapid,newpip)
****************************************************************
*           secte : detector sector                            *
*           torcur : torus current [A]
*           p_pip  : uncorrected (+) pion momentum             *
*  INPUT    phipid, thetapid : dimension : [degree]            *
*                           : uncorrected pion(+)  angles      *
*--------------------------------------------------------------*
*  OUTPUT   thetapi_corr : corrected pion(+) angle [degree]    *
*           ppi_corr     : corrected pion(+) momenutm          *
****************************************************************

      REAL nthetapid,mom_pip
      REAL phipi,thetapi,phipid,phipid2,thetapid
      REAL p_pip(10),e_pip(10),px_pip(10),py_pip(10),pz_pip(10)
      REAL torcur,tormax,Bratio
      DATA tormax/3375.0/

      INTEGER secth

C     ########### pion(+) angles correction constants ##########

      REAL piconst,pilinear,piquadratic,pithirdorder,pifourthorder,pinew

      REAL as1,as2,as3,as4,as5,as6
      REAL bs1,bs2,bs3,bs4,bs5,bs6
      REAL cs1,cs2,cs3,cs4,cs5,cs6

      REAL las1,las2,las3,las4,las5,las6
      REAL lbs1,lbs2,lbs3,lbs4,lbs5,lbs6
      REAL lcs1,lcs2,lcs3,lcs4,lcs5,lcs6

      REAL qas1,qas2,qas3,qas4,qas5,qas6
      REAL qbs1,qbs2,qbs3,qbs4,qbs5,qbs6
      REAL qcs1,qcs2,qcs3,qcs4,qcs5,qcs6

      REAL has1,has2,has3,has4,has5,has6
      REAL hbs1,hbs2,hbs3,hbs4,hbs5,hbs6
      REAL hcs1,hcs2,hcs3,hcs4,hcs5,hcs6

      REAL tas1,tas2,tas3,tas4,tas5,tas6
      REAL tbs1,tbs2,tbs3,tbs4,tbs5,tbs6
      REAL tcs1,tcs2,tcs3,tcs4,tcs5,tcs6

C23456789012345678901234567890123456789012345678901234567890123456789

      DATA tas1,tas2,tas3,tas4,tas5,tas6/-5.1299142E-03,-2.0352014E-03,
     &-2.826383E-04,-1.58515365E-03,-9.6490117E-04,-8.60422808E-04/
      DATA has1,has2,has3,has4,has5,has6/6.37836137E-04,-3.910457E-04,
     &1.17432835E-03,7.9181789E-04,2.44439224E-03,3.71918054E-05/
      DATA qas1,qas2,qas3,qas4,qas5,qas6/2.72359806E-03,-1.8383547E-03,
     &1.71147716E-03,8.92225567E-05,-4.205488E-03,3.39126187E-03/
      DATA las1,las2,las3,las4,las5,las6/-1.13440559E-03,4.10083333E-03,
     &1.1502628E-02,8.71534207E-03,6.3294693E-04,1.38427978E-02/
      DATA as1,as2,as3,as4,as5,as6/-5.16479755E-03,8.16171478E-04,
     &-3.45627949E-03,-2.6545698E-03,1.33297453E-03,-3.70724347E-03/

      DATA tbs1,tbs2,tbs3,tbs4,tbs5,tbs6/2.24826958E-01,9.0687908E-02,
     &1.7845019E-02,8.09892278E-02,4.9708738E-02,4.63797833E-02/
      DATA hbs1,hbs2,hbs3,hbs4,hbs5,hbs6/-3.31684806E-02,1.2504244E-02,
     &-5.12792186E-02,-3.6355055E-02,-1.08187405E-01,-7.9916894E-03/
      DATA qbs1,qbs2,qbs3,qbs4,qbs5,qbs6/-1.221424E-01,5.6957746E-02,
     &-6.00117927E-02,-1.39340483E-02,1.72421567E-01,-1.28512835E-01/
      DATA lbs1,lbs2,lbs3,lbs4,lbs5,lbs6/4.20764036E-02,-1.61183679E-01,
     &-4.27954855E-01,-3.10226155E-01,-1.8280028E-02,-5.06583454E-01/
      DATA bs1,bs2,bs3,bs4,bs5,bs6/2.07143724E-01,-4.2716178E-02,
     &1.09249069E-01,1.0517602E-01,-2.003366E-02,1.59486435E-01/


      DATA tcs1,tcs2,tcs3,tcs4,tcs5,tcs6/-2.44489241,-1.0028315,
     &-2.7030428E-01,-9.95714015E-01,-6.3760666E-01,-6.0376152E-01/
      DATA hcs1,hcs2,hcs3,hcs4,hcs5,hcs6/4.07094475E-01,-6.89544888E-02,
     &5.47101163E-01,4.23272698E-01,1.20123351,1.41910734E-01/
      DATA qcs1,qcs2,qcs3,qcs4,qcs5,qcs6/1.31174283,-5.1219438E-01,
     &3.5334639E-01,1.62542337E-01,-1.77085757,9.3981116E-01/
      DATA lcs1,lcs2,lcs3,lcs4,lcs5,lcs6/-2.8053815E-01,1.66281611,
     &3.4657869,2.3246248,4.5728914E-02,4.00965538/
      DATA cs1,cs2,cs3,cs4,cs5,cs6/-2.09355677,5.64114112E-01,
     &-7.70204323E-01,-1.14963156,-5.01741685E-01,-1.83039564/


C     ########### pi^+ momentum correction constants ##########

      REAL hconst,hlinear,hquadratic,hthirdorder,hfourthorder
      REAL hnew,newpip

      REAL ppas1,ppas2,ppas3,ppas4,ppas5,ppas6
      REAL ppbs1,ppbs2,ppbs3,ppbs4,ppbs5,ppbs6
      REAL ppcs1,ppcs2,ppcs3,ppcs4,ppcs5,ppcs6

      REAL pplas1,pplas2,pplas3,pplas4,pplas5,pplas6
      REAL pplbs1,pplbs2,pplbs3,pplbs4,pplbs5,pplbs6
      REAL pplcs1,pplcs2,pplcs3,pplcs4,pplcs5,pplcs6

      REAL ppqas1,ppqas2,ppqas3,ppqas4,ppqas5,ppqas6
      REAL ppqbs1,ppqbs2,ppqbs3,ppqbs4,ppqbs5,ppqbs6
      REAL ppqcs1,ppqcs2,ppqcs3,ppqcs4,ppqcs5,ppqcs6

      REAL pptas1,pptas2,pptas3,pptas4,pptas5,pptas6
      REAL pptbs1,pptbs2,pptbs3,pptbs4,pptbs5,pptbs6
      REAL pptcs1,pptcs2,pptcs3,pptcs4,pptcs5,pptcs6

      REAL pphas1,pphas2,pphas3,pphas4,pphas5,pphas6
      REAL pphbs1,pphbs2,pphbs3,pphbs4,pphbs5,pphbs6
      REAL pphcs1,pphcs2,pphcs3,pphcs4,pphcs5,pphcs6

C23456789012345678901234567890123456789012345678901234567890123456789


      DATA pphas1,pphas2,pphas3,pphas4,pphas5,pphas6/
     &-7.1844219079881D-05,-4.8440373290278D-05,-7.5425590169320D-05,
     &-4.4554868249908D-05,-1.2137659841733D-04,-7.0998407879210D-05/
      DATA pptas1,pptas2,pptas3,pptas4,pptas5,pptas6/
     &6.3269024406430D-06,-2.4069056941098D-05,-7.0215912241867D-05,
     &-2.0278991993274D-05,1.3774882073097D-04,1.7314651554110D-05/
      DATA ppqas1,ppqas2,ppqas3,ppqas4,ppqas5,ppqas6/
     &6.9390123917417D-05,-2.7535249014628D-06,1.1662582466434D-04,
     &4.3454299924587D-05,8.7070678161577D-05,5.6084084931323D-05/
      DATA pplas1,pplas2,pplas3,pplas4,pplas5,pplas6/
     &-1.0869417734953D-04,-5.8143222676711D-05,1.9792243033956D-04,
     &1.1593141811510D-04,-3.7830937738190D-05,4.9455233402118D-05/
      DATA ppas1,ppas2,ppas3,ppas4,ppas5,ppas6/
     &-4.3424651403068D-07,1.0630214597195D-06,-3.5053116445780D-07,
     &-5.8560815086425D-07,-2.5695606257633D-07,-4.0057393727035D-07/


      DATA pphbs1,pphbs2,pphbs3,pphbs4,pphbs5,pphbs6/
     &4.6949962334791D-03,3.2695665809276D-03,4.9441619886524D-03,
     &2.8710131678481D-03,7.8106249434667D-03,4.4864346138715D-03/
      DATA pptbs1,pptbs2,pptbs3,pptbs4,pptbs5,pptbs6/
     &-1.3143424370788D-04,1.7806239924986D-03,4.6940828938751D-03,
     &1.4721139409307D-03,-8.6812638324776D-03,-8.1883768503485D-04/
      DATA ppqbs1,ppqbs2,ppqbs3,ppqbs4,ppqbs5,ppqbs6/
     &-4.9511941857830D-03,-2.6993010700726D-04,-7.6944128537099D-03,
     &-2.9053080429949D-03,-6.3117358280142D-03,-4.0186323969971D-03/
      DATA pplbs1,pplbs2,pplbs3,pplbs4,pplbs5,pplbs6/
     &8.8234617050698D-03,4.5095253789595D-03,-1.3832456668735D-02,
     &-7.9456250292464D-03,3.1575061611213D-03,-4.7771875186452D-03/
      DATA ppbs1,ppbs2,ppbs3,ppbs4,ppbs5,ppbs6/
     &2.5117588263152D-05,-8.8638293953487D-05,7.3721662420257D-06,
     &4.4341862119425D-05,1.2818991016199D-05,3.5424024673774D-05/

      DATA pphcs1,pphcs2,pphcs3,pphcs4,pphcs5,pphcs6/
     &-7.9400009868257D-02,-5.9553204368528D-02,-8.4112643988478D-02,
     &-4.8536135903168D-02,-0.12494311516901,-7.4126030849555D-02/
      DATA pptcs1,pptcs2,pptcs3,pptcs4,pptcs5,pptcs6/
     &-8.0304529031240D-03,-3.8842746925294D-02,-8.1306337783908D-02,
     &-3.0499202793116D-02,0.13249260356287,2.1991659224153D-03/
      DATA ppqcs1,ppqcs2,ppqcs3,ppqcs4,ppqcs5,ppqcs6/
     &0.10637200070218,3.4246805206286D-02,0.14194708949095,
     &5.8621184194944D-02,0.11934450352563,9.7327218289753D-02/
      DATA pplcs1,pplcs2,pplcs3,pplcs4,pplcs5,pplcs6/
     &-0.14511389911125,-5.2945117280107D-02,0.27043480204968,
     &0.15979137314589,-4.7681356529085D-02,0.15144397990750/
      DATA ppcs1,ppcs2,ppcs3,ppcs4,ppcs5,ppcs6/
     &9.9380861161087D-02,1.0165916607859D-01,1.0024929124968D-01,
     &9.8929456850176D-02,9.9827501650252D-02,9.9153515589292D-02/
 
c
C#########################################################################


C     ##################  PHI&THETA PION ANGLE CORRECTION  #############BEGIN
      if(secth.eq.1)then 
         phipid2=phipid
      else if(secth.eq.2)then 
         phipid2=phipid-60.
      else if(secth.eq.3)then 
         phipid2=phipid-120.
      else if(secth.eq.4)then 
         IF(phipid.ge.150.and.phipid.le.180)THEN
            phipid2=phipid-180
         ELSE IF(phipid.ge.-180.and.phipid.le.-150)THEN
            phipid2=phipid+180
         ENDIF
      else if(secth.eq.5)then 
         phipid2=phipid+120.
      else if(secth.eq.6)then 
         phipid2=phipid+60.
      endif
                  
      if(thetapid.gt.12.and.thetapid.lt.25.)then

       if(secth.eq.1)then 
        if(phipid2.gt.-25.0.and.phipid2.lt.25.0)then ! constraint phi angle
         const=(cs1+bs1*thetapid+as1*thetapid**2)/10.
         linear=(lcs1+lbs1*thetapid+las1*thetapid**2)/100.
         quadratic=(qcs1+qbs1*thetapid+qas1*thetapid**2)/1000.
         thirdorder=(tcs1+tbs1*thetapid+tas1*thetapid**2)/10000.
         fourthorder=(hcs1+hbs1*thetapid+has1*thetapid**2)/100000.

         pinew=(const+linear*phipid2+quadratic*phipid2**2
     &        +thirdorder*phipid2**3+fourthorder*phipid2**4)
         nthetapid=thetapid+pinew
        endif  ! constraint phi angle

       else if(secth.eq.2)then 
        if(phipid2.gt.-25.0.and.phipid2.lt.25.0)then ! constraint phi angle
         const=(cs2+bs2*thetapid+as2*thetapid**2)/10.
         linear=(lcs2+lbs2*thetapid+las2*thetapid**2)/100.
         quadratic=(qcs2+qbs2*thetapid+qas2*thetapid**2)/1000.
         thirdorder=(tcs2+tbs2*thetapid+tas2*thetapid**2)/10000.
         fourthorder=(hcs2+hbs2*thetapid+has2*thetapid**2)/100000.

         pinew=(const+linear*phipid2+quadratic*phipid2**2
     &        +thirdorder*phipid2**3+fourthorder*phipid2**4)
         nthetapid=thetapid+pinew
        endif  ! constraint phi angle

	 
       else if(secth.eq.3)then 
        if(phipid2.gt.-25.0.and.phipid2.lt.25.0)then ! constraint phi angle
         const=(cs3+bs3*thetapid+as3*thetapid**2)/10.
         linear=(lcs3+lbs3*thetapid+las3*thetapid**2)/100.
         quadratic=(qcs3+qbs3*thetapid+qas3*thetapid**2)/1000.
         thirdorder=(tcs3+tbs3*thetapid+tas3*thetapid**2)/10000.
         fourthorder=(hcs3+hbs3*thetapid+has3*thetapid**2)/100000.

         pinew=(const+linear*phipid2+quadratic*phipid2**2
     &        +thirdorder*phipid2**3+fourthorder*phipid2**4)
         nthetapid=thetapid+pinew
        endif  ! constraint phi angle


       else if(secth.eq.4)then 
        if(phipid2.gt.-25.0.and.phipid2.lt.25.0)then ! constraint phi angle
         const=(cs4+bs4*thetapid+as4*thetapid**2)/10.
         linear=(lcs4+lbs4*thetapid+las4*thetapid**2)/100.
         quadratic=(qcs4+qbs4*thetapid+qas4*thetapid**2)/1000.
         thirdorder=(tcs4+tbs4*thetapid+tas4*thetapid**2)/10000.
         fourthorder=(hcs4+hbs4*thetapid+has4*thetapid**2)/100000.

         pinew=(const+linear*phipid2+quadratic*phipid2**2
     &        +thirdorder*phipid2**3+fourthorder*phipid2**4)
         nthetapid=thetapid+pinew
        endif  ! constraint phi angle


       else if(secth.eq.5)then 
        if(phipid2.gt.-25.0.and.phipid2.lt.25.0)then ! constraint phi angle
         const=(cs5+bs5*thetapid+as5*thetapid**2)/10.
         linear=(lcs5+lbs5*thetapid+las5*thetapid**2)/100.
         quadratic=(qcs5+qbs5*thetapid+qas5*thetapid**2)/1000.
         thirdorder=(tcs5+tbs5*thetapid+tas5*thetapid**2)/10000.
         fourthorder=(hcs5+hbs5*thetapid+has5*thetapid**2)/100000.

         pinew=(const+linear*phipid2+quadratic*phipid2**2
     &        +thirdorder*phipid2**3+fourthorder*phipid2**4)
         nthetapid=thetapid+pinew
        endif  ! constraint phi angle


       else if(secth.eq.6)then 
        if(phipid2.gt.-25.0.and.phipid2.lt.25.0)then ! constraint phi angle
         const=(cs6+bs6*thetapid+as6*thetapid**2)/10.
         linear=(lcs6+lbs6*thetapid+las6*thetapid**2)/100.
         quadratic=(qcs6+qbs6*thetapid+qas6*thetapid**2)/1000.
         thirdorder=(tcs6+tbs6*thetapid+tas6*thetapid**2)/10000.
         fourthorder=(hcs6+hbs6*thetapid+has6*thetapid**2)/100000.

         pinew=(const+linear*phipid2+quadratic*phipid2**2
     &        +thirdorder*phipid2**3+fourthorder*phipid2**4)
         nthetapid=thetapid+pinew
        endif  ! constraint phi angle

       endif

      else if(thetapid.le.12..or.thetapid.ge.25.)then

         if(secth.eq.1)then 
            nthetapid = thetapid
         else if(secth.eq.2)then 
            nthetapid = thetapid
         else if(secth.eq.3)then 
            nthetapid = thetapid
         else if(secth.eq.4)then 
            nthetapid = thetapid
         else if(secth.eq.5)then 
            nthetapid = thetapid
         else if(secth.eq.6)then 
            nthetapid = thetapid
         endif
         
      endif                     ! electronn theta angle limitation

      nthetapid = nthetapid

C     ##################  PHI&THETA PION ANGLE CORRECTION  #############END

C     ************** magnetic field consideration **************

      Bratio = torcur/tormax


C     ##########  PHI&THETA PIONS MOMENTUM CORRECTION  #############BEGIN
      newpip=mom_pip

      if(nthetapid.ge.12..and.nthetapid.le.30.)then
                  
        if(secth.eq.1)then 
           if(phipid2.gt.-25..and.phipid2.lt.25.)then
         hconst=(ppcs1+ppbs1*nthetapid+ppas1*nthetapid**2)*10.
         hlinear=(pplcs1+pplbs1*nthetapid+pplas1*nthetapid**2)/100.
         hquadratic=(ppqcs1+ppqbs1*nthetapid+ppqas1*nthetapid**2)/1000.
         hthirdorder=(pptcs1+pptbs1*nthetapid+pptas1*nthetapid**2)
     &        /10000.
        hfourthorder=(pphcs1+pphbs1*nthetapid+pphas1*nthetapid**2)
     &        /100000.

         hnew=(hconst+hlinear*phipid2+hquadratic*phipid2**2
     &        +hthirdorder*phipid2**3+hfourthorder*phipid2**4)
         newpip=mom_pip*hnew*Bratio
            else
              newpip=mom_pip 
            endif

        else if(secth.eq.2)then 
           if(phipid2.gt.-25..and.phipid2.lt.25.)then
         hconst=(ppcs2+ppbs2*nthetapid+ppas2*nthetapid**2)*10.
         hlinear=(pplcs2+pplbs2*nthetapid+pplas2*nthetapid**2)/100.
         hquadratic=(ppqcs2+ppqbs2*nthetapid+ppqas2*nthetapid**2)/1000.
         hthirdorder=(pptcs2+pptbs2*nthetapid+pptas2*nthetapid**2)
     &        /10000.
        hfourthorder=(pphcs2+pphbs2*nthetapid+pphas2*nthetapid**2)
     &        /100000.

         hnew=(hconst+hlinear*phipid2+hquadratic*phipid2**2
     &        +hthirdorder*phipid2**3+hfourthorder*phipid2**4)
         newpip=mom_pip*hnew*Bratio
            else
              newpip=mom_pip 
            endif
 
        else if(secth.eq.3)then 
           if(phipid2.gt.-25..and.phipid2.lt.25.)then
         hconst=(ppcs3+ppbs3*nthetapid+ppas3*nthetapid**2)*10.
         hlinear=(pplcs3+pplbs3*nthetapid+pplas3*nthetapid**2)/100.
         hquadratic=(ppqcs3+ppqbs3*nthetapid+ppqas3*nthetapid**2)/1000.
         hthirdorder=(pptcs3+pptbs3*nthetapid+pptas3*nthetapid**2)
     &        /10000.
        hfourthorder=(pphcs3+pphbs3*nthetapid+pphas3*nthetapid**2)
     &        /100000.

         hnew=(hconst+hlinear*phipid2+hquadratic*phipid2**2
     &        +hthirdorder*phipid2**3+hfourthorder*phipid2**4)
         newpip=mom_pip*hnew*Bratio
            else
              newpip=mom_pip 
            endif

        else if(secth.eq.4)then 
           if(phipid2.gt.-25..and.phipid2.lt.25.)then
         hconst=(ppcs4+ppbs4*nthetapid+ppas4*nthetapid**2)*10.
         hlinear=(pplcs4+pplbs4*nthetapid+pplas4*nthetapid**2)/100.
         hquadratic=(ppqcs4+ppqbs4*nthetapid+ppqas4*nthetapid**2)/1000.
         hthirdorder=(pptcs4+pptbs4*nthetapid+pptas4*nthetapid**2)
     &        /10000.
        hfourthorder=(pphcs4+pphbs4*nthetapid+pphas4*nthetapid**2)
     &        /100000.

         hnew=(hconst+hlinear*phipid2+hquadratic*phipid2**2
     &        +hthirdorder*phipid2**3+hfourthorder*phipid2**4)
         newpip=mom_pip*hnew*Bratio
            else
              newpip=mom_pip 
            endif

        else if(secth.eq.5)then 
           if(phipid2.gt.-25..and.phipid2.lt.25.)then
         hconst=(ppcs5+ppbs5*nthetapid+ppas5*nthetapid**2)*10.
         hlinear=(pplcs5+pplbs5*nthetapid+pplas5*nthetapid**2)/100.
         hquadratic=(ppqcs5+ppqbs5*nthetapid+ppqas5*nthetapid**2)/1000.
         hthirdorder=(pptcs5+pptbs5*nthetapid+pptas5*nthetapid**2)
     &        /10000.
        hfourthorder=(pphcs5+pphbs5*nthetapid+pphas5*nthetapid**2)
     &        /100000.

         hnew=(hconst+hlinear*phipid2+hquadratic*phipid2**2
     &        +hthirdorder*phipid2**3+hfourthorder*phipid2**4)
         newpip=mom_pip*hnew*Bratio
            else
              newpip=mom_pip 
            endif

        else if(secth.eq.6)then 
           if(phipid2.gt.-25..and.phipid2.lt.25.)then
         hconst=(ppcs6+ppbs6*nthetapid+ppas6*nthetapid**2)*10.
         hlinear=(pplcs6+pplbs6*nthetapid+pplas6*nthetapid**2)/100.
         hquadratic=(ppqcs6+ppqbs6*nthetapid+ppqas6*nthetapid**2)/1000.
         hthirdorder=(pptcs6+pptbs6*nthetapid+pptas6*nthetapid**2)
     &        /10000.
         hfourthorder=(pphcs6+pphbs6*nthetapid+pphas6*nthetapid**2)
     &        /100000.

         hnew=(hconst+hlinear*phipid2+hquadratic*phipid2**2
     &        +hthirdorder*phipid2**3+hfourthorder*phipid2**4)
         newpip=mom_pip*hnew*Bratio
            else
              newpip=mom_pip 
            endif

        endif
      else
         
         if(secth.eq.1)then 
            newpip=mom_pip
         else if(secth.eq.2)then 
            newpip=mom_pip
         else if(secth.eq.3)then 
            newpip=mom_pip
         else if(secth.eq.4)then 
            newpip=mom_pip
         else if(secth.eq.5)then 
            newpip=mom_pip
         else if(secth.eq.6)then 
            newpip=mom_pip
         endif
         
      endif

C     ##########  PHI&THETA PIONS MOMENTUM CORRECTION  #############END
      newpip=newpip
      return
      end


