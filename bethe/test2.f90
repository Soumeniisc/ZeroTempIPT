! it include the global and print w, that was created in another subroutine but showing error.
program IPT
use grid
include "global"


integer i

Wmax  = 20.0
Wmin  = -20.0
call omega()
	do i=-nbin/2,nbin/2-1
         	print*, w(i)  
      	enddo


end program IPT
