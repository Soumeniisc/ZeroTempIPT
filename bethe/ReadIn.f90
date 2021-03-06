! it created on 19/05/18
! it contain  ReadParams(), readSigma() and omega(),  function fermic(x,beta)
! last modified 21/05/18

subroutine ReadParams()

	use Global
	character(len=100) :: buffer, label
	!character(len=30) :: delim
  	integer :: pos
	integer, parameter :: fh = 15
  	integer :: ios = 0
  	integer :: line = 0
	
	! set the default values of parameters 
	Wmax  = 20.0
	Wmin  = -20.0
	
	t = 1.0
	t2 = 0.0
	delta = 1.0
	U = 0.0
	mu = 0.0	
	sub_lattice = 0

	kktransformation = 2
	benchmark = 1

	eta = 0.001
	MinNoLoop = 2
	MaxNoLoop = 2
	ToleranceDos = 0.0001	! Tolerance for convergence; may be increased if greater accuracy is required

	! ios is negative if an end of record condition is encountered or if
 	! an endfile condition was detected.  It is positive if an error was
  	! detected.  ios is zero otherwise.

    ! set the parameter value through PARAMS_IPT
	open(fh, file='PARAMS_IPT')
  	do while (ios == 0)

    	read(fh, '(A)', iostat=ios) buffer
     	if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, " ")
        label = buffer(1:pos)
        buffer = buffer(pos+1:)
		!print*, label, buffer
        select case (label)
        case ('t:')
           read(buffer, *, iostat=ios) t
           print *, 't:					     ', t
        case ('t2:')
           read(buffer, *, iostat=ios) t2
           print *, 't2:					 ', t2
		case ('delta:')
           read(buffer, *, iostat=ios) delta
           print *, 'delta:					 ', delta
		case ('U:')
           read(buffer, *, iostat=ios) U
           print *, 'U:					     ', U
        case ('mu:')
           read(buffer, *, iostat=ios) mu
           print *, 'mu:					 ', mu
		case ('sub_lattice:')
           read(buffer, *, iostat=ios) sub_lattice
           print *, 'sub_lattice:			 ', sub_lattice
		case ('kktransformation:')
           read(buffer, *, iostat=ios) kktransformation
           print *, 'kktransformation:		  ', kktransformation
		case ('benchmark:')
           read(buffer, *, iostat=ios) benchmark
           print *, 'benchmark:				  ', benchmark
		case ('eta:')
           read(buffer, *, iostat=ios) eta
           print *, 'eta:					  ', eta
        case ('MinNoLoop:')
           read(buffer, *, iostat=ios) MinNoLoop
           print *, 'MinNoLoop:				  ', MinNoLoop
		case ('MaxNoLoop:')
           read(buffer, *, iostat=ios) MaxNoLoop
           print *, 'MaxNoLoop:			      ', MaxNoLoop
		case ('ToleranceDos:')
           read(buffer, *, iostat=ios) ToleranceDos
           print *, 'ToleranceDos:			  ', ToleranceDos
        case default
			
           print *, 'Skipping invalid label at line of ', line
        end select
     end if
  end do

  	call omega()

end subroutine ReadParams

subroutine ReadSigma()
	
	use global
	implicit none

	integer i
	double precision omega, ReSigmaUp,ImSigmaUp,ReSigmaDn,ImSigmaDn

	!reading SigmaA
	open(unit=1, file="SigA.out", status="unknown")
	do i=-nbin/2,nbin/2-1
         	read(1,*) omega,ReSigmaUp,ImSigmaUp,ReSigmaDn,ImSigmaDn
		sigma_a(1,i) = dcmplx(ReSigmaUp,ImSigmaUp)
		sigma_a(2,i) = dcmplx(ReSigmaDn,ImSigmaDn)
      	enddo 	
	close(1)

	!reading SigmaB
	open(unit=2, file="SigB.out", status="unknown")
	do i=-nbin/2,nbin/2-1
         	read(2,*) omega,ReSigmaUp,ImSigmaUp,ReSigmaDn,ImSigmaDn
		sigma_b(1,i) = dcmplx(ReSigmaUp,ImSigmaUp)
		sigma_b(2,i) = dcmplx(ReSigmaDn,ImSigmaDn)
      	enddo 
	close(2)
	
	
	if (benchmark == 1) then
	! printing SigA.out to SigA.out.out to make sure that data is correctly read 
	open(unit=3,file='SigA.out.out',status='unknown')
	open(unit=4,file='SigB.out.out',status='unknown')
	do i=-nbin/2,nbin/2-1
         	write(3,*) w(i), real(sigma_a(1,i)), imag(sigma_a(1,i)), real(sigma_a(2,i)), imag(sigma_a(2,i))
		write(4,*) w(i), real(sigma_b(1,i)), imag(sigma_b(1,i)), real(sigma_b(2,i)), imag(sigma_b(2,i))
      	enddo
	close(3)
	close(4)
	end if !if (benchmark == 1)
	
	!print*, "self energy is read"
end

subroutine omega()
	!include "global"
	!its creates uniform frequency grid.

   	use Global
	integer i
    double precision beta 
	
	beta = 0.d0
   	binsize=(Wmax-Wmin)/dfloat(nbin)
	
   	do i=-nbin/2,nbin/2-1
      	w(i)=(dfloat(i)+0.5d0)*binsize  
		fermi(i) =  fermic(w(i),beta)
    enddo
	if (benchmark == 1) then
		print*, "look  at grid and fermi function in grid_fermi.dat file. printed from grid.f90/omega()"
		open(unit=238,file='grid_fermi.dat',status='unknown')
		do i=-nbin/2,nbin/2-1
			write(238,*)i,w(i),fermi(i)
		end do
		close(238)
	end if

end


double precision function fermic(x,beta)
	double precision x,beta
!	fermi function which avoids underflow/overflow.  If beta=0,
!	it is assumed that T=0 is meant!

	if(beta.eq.0.0D0) then
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0
	  else if(x.eq.0.0D0) then
	    fermic=0.5D0
	  else
	    fermic=0.0D0
	  end if
	else
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0/(dexp(beta*x)+1.0D0)
	  else
	    fermic=dexp(-beta*x)/(dexp(-beta*x)+1.0D0)
	  end if
	end if
	return
end

! this is a test function to check that function call work fine 
function area(a,b)
	real area
	real, intent(in):: a,b
	print*, "function test"
	area =  a*b
end function area


