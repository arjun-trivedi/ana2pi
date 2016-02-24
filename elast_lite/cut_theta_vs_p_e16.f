       subroutine theta_vs_p_el(sect,theta,p,good,at_mod)
       implicit none
       real theta,p,up,down,up2,down2,up3,down3
       real up5_1,down5_1,up5_2,down5_2
       real up6_1,down6_1
       integer good,sect
       integer at_mod
       
       good=1
       
       if (sect.eq.2.and.p.gt.1.5) then
       
       up=  35.5-4.*sqrt(p-1.5)
       down=34.-4.*sqrt(p-1.5)
       
       up2=   27.3-4.8*sqrt(p-1.5)
       down2= 26.7-4.8*sqrt(p-1.5)
       
       up3=   28.5-4.8*sqrt(p-1.5)
       down3= 26.1-4.8*sqrt(p-1.5)

c      at_mod 
       if (at_mod.ne.1) then
         if(theta.lt.up.and.theta.gt.down.and.p.gt.2.0.and.p.lt.3.05) good=0
         if(theta.lt.up2.and.theta.gt.down2.and.p.gt.2.6.and.p.lt.4)  good=0
         if(theta.lt.up3.and.theta.gt.down3.and.p.gt.2.6.and.p.lt.4)  good=0       
       endif

       endif


       if (sect.eq.5.and.p.gt.2.) then      
       up5_1 = 20.5
       down5_1 = 21.3-2.4*sqrt(p-2.)
       
       up5_2 =   26.3
       down5_2 = 24.7
       
       
       if(theta.lt.up5_1.and.theta.gt.down5_1.and.p.gt.2.9.and.p.lt.4.2)  good=0       
       if(theta.lt.up5_2.and.theta.gt.down5_2.and.p.gt.2.4.and.p.lt.3.4)  good=0       
       
       endif
       
       if (sect.eq.6.and.p.gt.1.5) then      
       up6_1 =31.-5.*sqrt(p-1.5)
       down6_1 = 29.7-5.*sqrt(p-1.5)
       
     
c      at_mod 
       if (at_mod.ne.1) then
         if(theta.lt.up6_1.and.theta.gt.down6_1.and.p.gt.2.4.and.p.lt.3.7)  good=0       
       endif

       
       endif
        
       
       
       
       return
       end
       
 	subroutine theta_vs_p_pr(sect,theta,p,good,at_mod)
 	implicit none
        real theta,p,up,down,up2,down2,up3,down3
        real up5_1,down5_1,up5_2,down5_2,up5_3,down5_3
        real up6_1,down6_1
        integer good,sect
        integer at_mod
	
       good=1
       
       if (sect.eq.2) then
       up=  27.2
       down=24.7
       
       if(theta.lt.up.and.theta.gt.down.and.p.gt.0.6.and.p.lt.4) good=0
       
       endif
       

       if (sect.eq.3.and.p.gt.0.9) then
       up2=  37.7+6.*sqrt(p-0.9)
       down2=31.5+9.*sqrt(p-0.9)

c      at_mod 
       if (at_mod.ne.1) then
         if(theta.lt.up2.and.theta.gt.down2.and.p.gt.0.9.and.p.lt.2) good=0
       endif
       if(p.lt.0.9) good=0
       
       endif
       
       
       if (sect.eq.5.and.p.gt.0.9) then
c       up5_1=25.2  
c       down5_2=24.
       
       up5_2=11.2  
       down5_2=10.1

       up5_3=11.+3.5*sqrt(p-0.9) 
       down5_3=13.5-0.38*(p-4)*(p-4)

       if(theta.lt.25.2.and.theta.gt.24..and.p.gt.1..and.p.lt.3.8) good=0       
       if(theta.lt.up5_2.and.theta.gt.down5_2.and.p.gt.1..and.p.lt.4.4) good=0
       if(theta.lt.up5_3.and.theta.gt.down5_3.and.p.gt.0.9.and.p.lt.4.0) good=0 
       
       endif
       
       if (sect.eq.6) then
       up6_1=29.2-1.9*(p-3.2)*(p-3.2)  
       down6_1=26.-1.3*(p-3.5)*(p-3.5)
       
      if(theta.lt.up6_1.and.theta.gt.down6_1.and.p.gt.0.9.and.p.lt.2.6) good=0       
        
       
       endif              
       
 


       
      return
       end
 
 	subroutine theta_vs_p_pip(sect,theta,p,good,at_mod)
 	implicit none
        real theta,p,up,down,up2,down2,up3,down3
        real up5,down5,up5_2,down5_2,up5_3,down5_3,up3_3
        real up6_1,down6_1,up2_1,down2_1,up3_1,down3_1,up3_2,down3_2
	real up4,down4,down6,up6_2,down6_2,up6_4,down6_4,up6_3,down6_3
        integer good,sect
        integer at_mod

       good=1
          
       if (sect.eq.2) then
       
       up= 0.75 + 0.003116*(theta-24.64)*(theta-24.64) 
       down=0.68 + 0.00207*(theta-24.96)*(theta-24.96)
       
       up2=-2.65+0.035*theta
       down2=-1.93 + 0.0235*theta
       
       up2_1=-2.07+0.0235*theta
       down2_1=-2.21 + 0.0235*theta
       
               
      
       if(p.lt.0.72.and.p.gt.0.64.and.theta.gt.6.and.theta.lt.26.) good=0	
       if(p.lt.up2.and.p.gt.down2.and.theta.gt.80.and.theta.lt.98.) good=0	
c      at_mod 
       if (at_mod.ne.1) then
         if(p.lt.up2_1.and.p.gt.down2_1.and.theta.gt.94.and.theta.lt.107.) good=0          	
       endif
       if(p.lt.up.and.p.gt.down.and.theta.gt.26.and.theta.lt.46.) good=0	

       endif
       
        if (sect.eq.3) then
	up3 = 0.7 + 0.003116*(theta-24.64)*(theta-24.64)
	down3 = 0.611 + 0.00207*(theta-24.96)*(theta-24.96)
	
	up3_1 = 0.308 + 0.001466*(theta-54.77)*(theta-54.77)
	down3_1 = 0.193 + 0.00084*(theta-54.04)*(theta-54.04)
	
	up3_2 = 0.19 + 0.00091*(theta-68.8)*(theta-68.8)
	down3_2 = -0.68 + 0.011*theta
	
	up3_3 = -2.02 + 0.0235*theta
	
	
c      at_mod 
       if (at_mod.ne.1) then
         if(p.lt.0.69.and.p.gt.0.6.and.theta.gt.6.and.theta.lt.26.3) good=0	
         if(p.lt.up3.and.p.gt.down3.and.theta.gt.26.and.theta.lt.46.) good=0
         if(p.lt.up3_1.and.p.gt.down3_1.and.theta.gt.52.and.theta.lt.79.) good=0
         if(p.lt.up3_2.and.p.gt.down3_2.and.theta.gt.73.and.theta.lt.98.) good=0       
                
         if(p.lt.up3_3.and.theta.gt.94.) good=0	 
       endif
	 
	 endif	      
       
       
        if (sect.eq.4) then
	up4   = 0.317 + 0.00105*(theta-42.6)*(theta-42.6)
	down4 = 0.266 + 0.000766*(theta-43.1)*(theta-43.1)
c      at_mod 
       if (at_mod.ne.1) then
         if(p.lt.up4.and.p.gt.down4.and.theta.gt.42.and.theta.lt.69.) good=0	
       endif
	
	endif      
       
        if (sect.eq.5) then
	up5   = 0.135 + 0.00064*(theta-64.2)*(theta-64.2)
	down5 = 0.101 + 0.0002*(theta-56.7)*(theta-56.7)
c      at_mod 
       if (at_mod.ne.1) then
         if(p.lt.up5.and.p.gt.down5.and.theta.gt.76.and.theta.lt.95.) good=0	
       endif
	
	endif      
       	
        if (sect.eq.6) then
	down6 = -0.12 + 0.06*theta 
	 up6_1 =   0.62 + 0.004*(theta-22.7)*(theta-22.7)
	 down6_1 = 0.47 + 0.0036*(theta-22.7)*(theta-22.7)
	 
	 up6_2 =   0.46 + 0.0044*(theta-32.8)*(theta-32.8)
	 down6_2 = 0.34 + 0.004*(theta-32.8)*(theta-32.8)
	 
	 up6_3 =   -1.075 + 0.015*theta
	 down6_3 = -1.25 + 0.015*theta
	 
	 up6_4 =   -0.64   + 0.027*theta
	 down6_4 = 0.25 + 0.00054*(theta-48.5)*(theta-48.5)
	 
	 

       if(p.gt.down6.and.theta.gt.8.and.theta.lt.36.) good=0
       if(p.lt.up6_1.and.p.gt.down6_1.and.theta.gt.23.and.theta.lt.37.) good=0       
       if(p.lt.up6_2.and.p.gt.down6_2.and.theta.gt.33.and.theta.lt.45.) good=0 	
c      at_mod 
       if (at_mod.ne.1) then
         if(p.lt.up6_3.and.p.gt.down6_3.and.theta.gt.82.and.theta.lt.106.) good=0
         if(p.lt.up6_4.and.p.gt.down6_4.and.theta.gt.35.5.and.theta.lt.77.) good=0       
       endif
        	       

	endif
	
	      
        return
        end
 
