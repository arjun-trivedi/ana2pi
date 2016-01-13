      Subroutine e_corr_sub
     &     (thetaeld,phield,pel,torcur,secte,thetaeldnew,newpel)
****************************************************************
*           secte : detector sector                            *
*           torcur : torus current [A]
*           pel  : uncorrected electron momentum               *
*  INPUT    phield, thetaeld : dimension : [degree]            *
*                           : uncorrected electron angles      *
*--------------------------------------------------------------*
*  OUTPUT   thetae_corr : corrected electron angle [degree]    *
*           pe_corr     : corrected electron momenutm          *
****************************************************************

      REAL pel,thetael,phiel,thetaeld,phield,phield2
      REAL tormax,torcur,Bratio
      DATA tormax/3375./

      INTEGER secte
      
C     ########### theta & phi correction constants ##########
      REAL const,linear,quadratic,thirdorder,fourthorder
      REAL new,thetaeldnew
      
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
      
C     23456789012345678901234567890123456789012345678901234567890123456789
      
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


C     ########### momentum correction constants ##########
      REAL sangsu,eilcha,jeggop,sejeggop,corr_fact,newpel

      REAL ps1,ps2,ps3,ps4,ps5,ps6
      REAL ks1,ks2,ks3,ks4,ks5,ks6
      REAL js1,js2,js3,js4,js5,js6

      REAL lps1,lps2,lps3,lps4,lps5,lps6
      REAL lks1,lks2,lks3,lks4,lks5,lks6
      REAL ljs1,ljs2,ljs3,ljs4,ljs5,ljs6

      REAL qps1,qps2,qps3,qps4,qps5,qps6
      REAL qks1,qks2,qks3,qks4,qks5,qks6
      REAL qjs1,qjs2,qjs3,qjs4,qjs5,qjs6

      REAL dps1,dps2,dps3,dps4,dps5,dps6
      REAL dks1,dks2,dks3,dks4,dks5,dks6
      REAL djs1,djs2,djs3,djs4,djs5,djs6

      DATA ps1,ps2,ps3,ps4,ps5,ps6/9.5112707858653D-08,
     &1.0033722952831D-07,1.0797150380251D-07,1.7100120165989D-07,
     &1.3082672597941D-08,8.7085945937647D-08/
      DATA ks1,ks2,ks3,ks4,ks5,ks6/-7.9599239534114D-06,
     &-5.6599480422299D-06,-5.0800069598611D-06,-1.1728924887767D-05,
     &-2.4757597510374D-06,-5.3018979653485D-06/
      DATA js1,js2,js3,js4,js5,js6/1.0176579118581D-02,
     &1.0087801836821D-02,1.0073933449548D-02,1.0221497899831D-02,
     &1.0110107335248D-02,1.0112945449548D-02/

      DATA lps1,lps2,lps3,lps4,lps5,lps6/-6.9400370675523D-05,
     &-1.5960488468465D-04,-2.1899832926782D-04,-6.3507481408119D-05,
     &-9.3393940238598D-05,-1.0356199140906D-04/
      DATA lks1,lks2,lks3,lks4,lks5,lks6/7.9655917387783D-03,
     &1.1661447536158D-02,1.1211638296420D-02,2.4612021589962D-03,
     &6.1878054901081D-03,4.6198689235790D-03/
      DATA ljs1,ljs2,ljs3,ljs4,ljs5,ljs6/-0.20495101625537,
     &-0.22412829804132,-1.0019870422334D-01,-2.0318433751035D-02,
     &-0.11955549846234,-1.8772394756644D-02/

      DATA qps1,qps2,qps3,qps4,qps5,qps6/-1.8086561701112D-04,
     &-1.8929500388527D-04,-1.1416945384378D-04,-9.1297828333589D-05,
     &-1.1784209351117D-04,5.3194937440755D-05/
      DATA qks1,qks2,qks3,qks4,qks5,qks6/1.0743244205188D-02,
     &1.0942470368942D-02,6.8099881808583D-03,5.2998598339562D-03,
     &6.8477708809656D-03,-3.2324591761590D-03/
      DATA qjs1,qjs2,qjs3,qjs4,qjs5,qjs6/-0.15348562918845,
     &-0.15076163008159,-9.9580311244930D-02,-7.6365686745569D-02,
     &-9.6360858729669D-02,3.1792598257448D-02/

      DATA dps1,dps2,dps3,dps4,dps5,dps6/1.3019443527636D-04,
     &1.2506972055652D-04,-7.3395343512598D-06,-1.4020818540861D-05,
     &4.4662480048095D-05,1.8214642092420D-05/
      DATA dks1,dks2,dks3,dks4,dks5,dks6/-7.9967328974761D-03,
     &-7.6499984847730D-03,2.5650746001025D-04,4.7628893452639D-04,
     &-2.7046265634257D-03,-1.3819042499386D-03/
      DATA djs1,djs2,djs3,djs4,djs5,djs6/0.11761993283213,
     &0.11317471787380,-1.2820491669835D-04,4.7326314122685D-04,
     &3.9383613408158D-02,2.6160738941275D-02/

C     ###########  ELECTRON PHI&THETA ANGLE CORRECTION  #########BEGIN
      if(thetaeld.gt.12.and.thetaeld.lt.25.)then
         
         if(secte.eq.1)then 
            phield2=phield
c          if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
            const=(cs1+bs1*thetaeld+as1*thetaeld**2)/10.
            linear=(lcs1+lbs1*thetaeld+las1*thetaeld**2)/100.
            quadratic=(qcs1+qbs1*thetaeld+qas1*thetaeld**2)/1000.
            thirdorder=(tcs1+tbs1*thetaeld+tas1*thetaeld**2)/10000.
            fourthorder=(hcs1+hbs1*thetaeld+has1*thetaeld**2)/100000.
            
            new=(const+linear*phield2+quadratic*phield2**2
     &           +thirdorder*phield2**3+fourthorder*phield2**4)
            thetaeldnew=thetaeld+new
c           endif        ! constraint phi angle 
            
         else if(secte.eq.2)then
            phield2=phield-60
c          if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
            const=(cs2+bs2*thetaeld+as2*thetaeld**2)/10.
            linear=(lcs2+lbs2*thetaeld+las2*thetaeld**2)/100.
            quadratic=(qcs2+qbs2*thetaeld+qas2*thetaeld**2)/1000.
            thirdorder=(tcs2+tbs2*thetaeld+tas2*thetaeld**2)/10000.
            fourthorder=(hcs2+hbs2*thetaeld+has2*thetaeld**2)/100000.

            new=(const+linear*phield2+quadratic*phield2**2
     &           +thirdorder*phield2**3+fourthorder*phield2**4)
            thetaeldnew=thetaeld+new
c           endif        ! constraint phi angle 

         else if(secte.eq.3)then
            phield2=phield-120
c           if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
           const=(cs3+bs3*thetaeld+as3*thetaeld**2)/10.
            linear=(lcs3+lbs3*thetaeld+las3*thetaeld**2)/100.
            quadratic=(qcs3+qbs3*thetaeld+qas3*thetaeld**2)/1000.
            thirdorder=(tcs3+tbs3*thetaeld+tas3*thetaeld**2)/10000.
            fourthorder=(hcs3+hbs3*thetaeld+has3*thetaeld**2)/100000.
            
            new=(const+linear*phield2+quadratic*phield2**2
     &           +thirdorder*phield2**3+fourthorder*phield2**4)
            thetaeldnew=thetaeld+new
c           endif        ! constraint phi angle 
            
         else if(secte.eq.4)then 
            IF(phield.ge.150.and.phield.le.180)THEN
               phield2=phield-180
            ELSE IF(phield.ge.-180.and.phield.le.-150)THEN
               phield2=phield+180
            ENDIF
c          if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
            const=(cs4+bs4*thetaeld+as4*thetaeld**2)/10.
            linear=(lcs4+lbs4*thetaeld+las4*thetaeld**2)/100.
            quadratic=(qcs4+qbs4*thetaeld+qas4*thetaeld**2)/1000.
            thirdorder=(tcs4+tbs4*thetaeld+tas4*thetaeld**2)/10000.
            fourthorder=(hcs4+hbs4*thetaeld+has4*thetaeld**2)/100000.
            
            new=(const+linear*phield2+quadratic*phield2**2
     &           +thirdorder*phield2**3+fourthorder*phield2**4)
            thetaeldnew=thetaeld+new
c           endif        ! constraint phi angle 
            
         else if(secte.eq.5)then
            phield2=phield+120
c          if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
            const=(cs5+bs5*thetaeld+as5*thetaeld**2)/10.
            linear=(lcs5+lbs5*thetaeld+las5*thetaeld**2)/100.
            quadratic=(qcs5+qbs5*thetaeld+qas5*thetaeld**2)/1000.
            thirdorder=(tcs5+tbs5*thetaeld+tas5*thetaeld**2)/10000.
            fourthorder=(hcs5+hbs5*thetaeld+has5*thetaeld**2)/100000.
            
            new=(const+linear*phield2+quadratic*phield2**2
     &           +thirdorder*phield2**3+fourthorder*phield2**4)
            thetaeldnew=thetaeld+new
c           endif        ! constraint phi angle 

         else if(secte.eq.6)then
            phield2=phield+60
c          if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
            const=(cs6+bs6*thetaeld+as6*thetaeld**2)/10.
            linear=(lcs6+lbs6*thetaeld+las6*thetaeld**2)/100.
            quadratic=(qcs6+qbs6*thetaeld+qas6*thetaeld**2)/1000.
            thirdorder=(tcs6+tbs6*thetaeld+tas6*thetaeld**2)/10000.
            fourthorder=(hcs6+hbs6*thetaeld+has6*thetaeld**2)/100000.
            
            new=(const+linear*phield2+quadratic*phield2**2
     &           +thirdorder*phield2**3+fourthorder*phield2**4)
            thetaeldnew=thetaeld+new
c           endif        ! constraint phi angle 
            
         endif
      else if(thetaeld.le.12..or.thetaeld.ge.25.)then
         if(secte.eq.1)then  
            phield2=phield 
            thetaeldnew = thetaeld
         elseif(secte.eq.2)then  
            phield2=phield-60 
            thetaeldnew = thetaeld
         elseif(secte.eq.3)then 
            phield2=phield-120 
            thetaeldnew = thetaeld
         elseif(secte.eq.4)then 
            IF(phield.ge.150.and.phield.le.180)THEN 
               phield2=phield-180 
            ELSE IF(phield.ge.-180.and.phield.le.-150)THEN 
               phield2=phield+180 
            ENDIF 
            thetaeldnew = thetaeld
         elseif(secte.eq.5)then 
            phield2=phield+120 
            thetaeldnew = thetaeld
         elseif(secte.eq.6)then 
            phield2=phield+60 
            thetaeldnew = thetaeld
         endif 
      endif                     ! proton angle limitation

      thetaeldnew=thetaeldnew

C     ######### ELECTRON PHI&THETA ANGLE CORRECTION  ##########END
                 

C     ************** magnetic field consideration **************

      Bratio = torcur/tormax
    
C     ************ Electron Momentum Correction Parts ************BEGIN*
      newpel=pel

      if(secte.eq.1)then 
       if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
         sangsu=(js1+ks1*thetaeldnew+ps1*thetaeldnew**2)*100.
         eilcha=(ljs1+lks1*thetaeldnew+lps1*thetaeldnew**2)/100.
         jeggop=(qjs1+qks1*thetaeldnew+qps1*thetaeldnew**2)/1000.
         sejeggop=(djs1+dks1*thetaeldnew+dps1*thetaeldnew**2)/10000.

         corr_fact=(sangsu+eilcha*phield2+
     +        jeggop*phield2**2+sejeggop*phield2**3)
         newpel=pel*corr_fact*Bratio
        endif        ! constraint phi angle 

         
      else if(secte.eq.2)then 
       if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
         sangsu=(js2+ks2*thetaeldnew+ps2*thetaeldnew**2)*100.
         eilcha=(ljs2+lks2*thetaeldnew+lps2*thetaeldnew**2)/100.
         jeggop=(qjs2+qks2*thetaeldnew+qps2*thetaeldnew**2)/1000.
         sejeggop=(djs2+dks2*thetaeldnew+dps2*thetaeldnew**2)/10000.

         corr_fact=(sangsu+eilcha*phield2+
     +        jeggop*phield2**2+sejeggop*phield2**3)
         newpel=pel*corr_fact*Bratio
        endif        ! constraint phi angle 

	 
      else if(secte.eq.3)then 
       if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
         sangsu=(js3+ks3*thetaeldnew+ps3*thetaeldnew**2)*100.
         eilcha=(ljs3+lks3*thetaeldnew+lps3*thetaeldnew**2)/100.
         jeggop=(qjs3+qks3*thetaeldnew+qps3*thetaeldnew**2)/1000.
         sejeggop=(djs3+dks3*thetaeldnew+dps3*thetaeldnew**2)/10000.

         corr_fact=(sangsu+eilcha*phield2+
     +        jeggop*phield2**2+sejeggop*phield2**3)
         newpel=pel*corr_fact*Bratio
        endif        ! constraint phi angle 


      else if(secte.eq.4)then 
       if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
         sangsu=(js4+ks4*thetaeldnew+ps4*thetaeldnew**2)*100.
         eilcha=(ljs4+lks4*thetaeldnew+lps4*thetaeldnew**2)/100.
         jeggop=(qjs4+qks4*thetaeldnew+qps4*thetaeldnew**2)/1000.
         sejeggop=(djs4+dks4*thetaeldnew+dps4*thetaeldnew**2)/10000.

         corr_fact=(sangsu+eilcha*phield2+
     +        jeggop*phield2**2+sejeggop*phield2**3)
         newpel=pel*corr_fact*Bratio
        endif        ! constraint phi angle 


      else if(secte.eq.5)then 
       if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
         sangsu=(js5+ks5*thetaeldnew+ps5*thetaeldnew**2)*100.
         eilcha=(ljs5+lks5*thetaeldnew+lps5*thetaeldnew**2)/100.
         jeggop=(qjs5+qks5*thetaeldnew+qps5*thetaeldnew**2)/1000.
         sejeggop=(djs5+dks5*thetaeldnew+dps5*thetaeldnew**2)/10000.

         corr_fact=(sangsu+eilcha*phield2+
     +        jeggop*phield2**2+sejeggop*phield2**3)
         newpel=pel*corr_fact*Bratio
        endif        ! constraint phi angle 


      else if(secte.eq.6)then 
       if(phield2.gt.-20.0.and.phield2.lt.20.0)then ! constraint phi angle
         sangsu=(js6+ks6*thetaeldnew+ps6*thetaeldnew**2)*100.
         eilcha=(ljs6+lks6*thetaeldnew+lps6*thetaeldnew**2)/100.
         jeggop=(qjs6+qks6*thetaeldnew+qps6*thetaeldnew**2)/1000.
         sejeggop=(djs6+dks6*thetaeldnew+dps6*thetaeldnew**2)/10000.

         corr_fact=(sangsu+eilcha*phield2+
     +        jeggop*phield2**2+sejeggop*phield2**3)
         newpel=pel*corr_fact*Bratio
        endif        ! constraint phi angle 
      endif
      newpel=newpel

      return
      end


