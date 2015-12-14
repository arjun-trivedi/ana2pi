       subroutine theta_vs_p_el(sect,theta,p,good)
       implicit none
       real theta,p,up,down,up2,down2,up3,down3
       real up5_1,down5_1,up5_2,down5_2
       real up6_1,down6_1
       integer good,sect
    
       
       good=1
       
       if (sect.eq.2.and.p.gt.1.5) then
       
       up=  35.5-4.*sqrt(p-1.5)
       down=34.-4.*sqrt(p-1.5)
       
       up2=   27.3-4.8*sqrt(p-1.5)
       down2= 26.7-4.8*sqrt(p-1.5)
       
       up3=   28.5-4.8*sqrt(p-1.5)
       down3= 26.1-4.8*sqrt(p-1.5)
       
       
       if(theta.lt.up.and.theta.gt.down.and.p.gt.2.0.and.p.lt.3.05) good=0
       if(theta.lt.up2.and.theta.gt.down2.and.p.gt.2.6.and.p.lt.4)  good=0
       if(theta.lt.up3.and.theta.gt.down3.and.p.gt.2.6.and.p.lt.4)  good=0       

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
       
     
       
       if(theta.lt.up6_1.and.theta.gt.down6_1.and.p.gt.2.4.and.p.lt.3.7)  good=0       

       
       endif
        
       
       
       
       return
       end
       
