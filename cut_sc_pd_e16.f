       subroutine myscint(paddle,bad_sc,at_mod)
*

*
*
*    OUTPUT: bad_sc, if 1 paddle is bad, if 0 baddle is O.K.
*

******************************************************************         
      implicit none
*
      integer paddle,bad_sc
      integer at_mod
*
c     [02-18-16] To check if Fortran:logical*1 is compatible with CPP:bool
c      if (paddle.eq.311) then
c        write(*,"(I1)") at_mod
c        if (at_mod.eq.1) then
c          write(*,*) "at_mod=1"
c        elseif (at_mod.eq.0) then
c          write(*,*) "at_mod=0"
c        else
c          write(*,*) "at_mod=unknown"
c        endif
c      endif

      bad_sc =0

c	print *,paddle	 
         if (paddle.eq.145) bad_sc=1

  
        
         if (paddle.eq.245) bad_sc=1
c        [02-18-16]Add at_mod paddles
         if (at_mod.eq.1) then
           if (paddle.eq.205) bad_sc=1
         end if

     
         if (paddle.eq.324) bad_sc=1
         if (paddle.eq.337) bad_sc=1
         if (paddle.eq.338) bad_sc=1
         if (paddle.eq.342) bad_sc=1
         if (paddle.eq.345) bad_sc=1
         if (paddle.eq.346) bad_sc=1
c        [02-18-16]Add at_mod paddles
         if (at_mod.eq.1) then
           if (paddle.eq.311) bad_sc=1
           if (paddle.eq.340) bad_sc=1
           if (paddle.eq.347) bad_sc=1
         end if
    
      
         if (paddle.eq.434) bad_sc=1

    

         if (paddle.eq.542) bad_sc=1
c        [02-18-16]Add at_mod paddles
         if (at_mod.eq.1) then
           if (paddle.eq.505) bad_sc=1
           if (paddle.eq.520) bad_sc=1
           if (paddle.eq.532) bad_sc=1
           if (paddle.eq.540) bad_sc=1
         end if


         
         if (paddle.eq.636) bad_sc=1
         if (paddle.eq.637) bad_sc=1
         if (paddle.eq.644) bad_sc=1
c       [02-18-16]Add at_mod paddles
         if (at_mod.eq.1) then
           if (paddle.eq.632) bad_sc=1
           if (paddle.eq.633) bad_sc=1
           if (paddle.eq.634) bad_sc=1
           if (paddle.eq.635) bad_sc=1
           if (paddle.eq.638) bad_sc=1
           if (paddle.eq.639) bad_sc=1
         end if

c       [02-18-16] debug
c        if (paddle.eq.311) then
c          write(*,*) "bad_sc:"
c          write(*,"(L1)") bad_sc
c        endif  


      return
      end      
