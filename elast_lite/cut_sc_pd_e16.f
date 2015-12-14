       subroutine myscint(paddle,bad_sc)
*

*
*
*    OUTPUT: bad_sc, if 1 paddle is bad, if 0 baddle is O.K.
*

******************************************************************         
      implicit none
*
      integer paddle,bad_sc
*


      bad_sc =0

c	print *,paddle	 
         if (paddle.eq.145) bad_sc=1

  
        
         if (paddle.eq.245) bad_sc=1


     
         if (paddle.eq.324) bad_sc=1
         if (paddle.eq.337) bad_sc=1
         if (paddle.eq.338) bad_sc=1
         if (paddle.eq.342) bad_sc=1
         if (paddle.eq.345) bad_sc=1
         if (paddle.eq.346) bad_sc=1
    
      
         if (paddle.eq.434) bad_sc=1

    

         if (paddle.eq.542) bad_sc=1

         
         if (paddle.eq.636) bad_sc=1
         if (paddle.eq.637) bad_sc=1
         if (paddle.eq.644) bad_sc=1

  


      return
      end      
