 	subroutine theta_vs_p_pr(sect,theta,p,good)
 	implicit none
        real theta,p,up,down,up2,down2,up3,down3
        real up5_1,down5_1,up5_2,down5_2,up5_3,down5_3
        real up6_1,down6_1
        integer good,sect
	
       good=1
       
       if (sect.eq.2) then
       up=  27.2
       down=24.7
       
       if(theta.lt.up.and.theta.gt.down.and.p.gt.0.6.and.p.lt.4) good=0
       
       endif
       

       if (sect.eq.3.and.p.gt.0.9) then
       up2=  37.7+6.*sqrt(p-0.9)
       down2=31.5+9.*sqrt(p-0.9)
       
       if(theta.lt.up2.and.theta.gt.down2.and.p.gt.0.9.and.p.lt.2) good=0
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
 
