
         
      Double precision function th_par(par,theta_pr) 
      implicit none
      integer par
      real theta_pr
      DOUBLE PRECISION thpar(7,5)
      data thpar / 
     *1.05029d-8,1.00,7.79734d-11,-1.14271d-12,5.23112d-15,0.0,0.0 ,
     *0.0553791,-0.0247709,3.27582d-5,-1.27335d-6,7.25717d-9,0.0,0.0,
     *-0.0127642,0.0119281,-4.17922d-5,4.55358d-7,-1.57514d-9,0.,0.,
     * -0.0482976,0.0183194,-0.000121302,2.03026d-6,-8.96678d-9,0.,0.,
     * 0.,0.,0.,0.,0.,0.,0. / 
      SAVE thpar
      
c      write
      
      th_par=((thpar(1,par)+(thpar(2,par)*theta_pr)+
     & (thpar(3,par)*theta_pr*theta_pr)+
     & (thpar(4,par)*theta_pr*theta_pr*theta_pr)+
     &(thpar(5,par)*
     & theta_pr*theta_pr*theta_pr*theta_pr))/theta_pr);
  
       return 
       end



       subroutine pr_ec(p_meas,theta_me,corrf)
       implicit none
       DOUBLE PRECISION th_par
       real p_meas,theta_me,corrf
       corrf =
     & th_par(1,theta_me)+
     &  (th_par(2,theta_me)/(p_meas))+
     &  (th_par(3,theta_me)/(p_meas*p_meas))+
     &  (th_par(4,theta_me)/Sqrt(p_meas))+
     &  (th_par(5,theta_me)/(p_meas*p_meas*p_meas*p_meas));
c      p_cor=corrf*p_meas

       return 
       end
       
       subroutine preloss(p_me,thet,pcorr)
       implicit none
       real p_me,thet,pcorr,corr,cor1
       
       corr=1.
       
      if (p_me.gt.0.25.and.p_me.lt.6.and.thet.gt.6.and.thet.lt.120) then
         if (p_me.lt.2.) then
	 call pr_ec(p_me,thet,corr)
	 pcorr=corr*p_me
	  else ! p_me.gt.2.
	  
      corr = 2*0.905*Sqrt((p_me*p_me)+0.938*0.938)/(p_me*p_me)
      call pr_ec(2.,thet,cor1)
      corr = 1 + (corr*(cor1-1))
      pcorr=corr*p_me	  
        endif	
	  else
	  pcorr=p_me 
      endif
c      pcorr=corr*p_me

       
       return
       end


      

