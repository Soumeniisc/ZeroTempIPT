! here we are creating a module for the global veriables and its defined in main program and will be used in all subroutines, including the main program.

! /home/soumen/INSTALLPACKAGES/mpich2/bin/mpif90 Global_mod.f90 grid.f90 ReadIn.f90 convolution.f90 ksum.f90 test3.f90

subroutine Initialisation()
	use Global
	implicit none
	integer i
	!real ab, area

	call ReadParams()


	!test function call
	!ab = area(3.,4.)
        !print*,"area", ab

	
	!writing in a file omega
	if (benchmark == 1) then
	open(unit=3,file='SigA.out',status='unknown')
	open(unit=4,file='SigB.out',status='unknown')
	do i=-nbin/2,nbin/2-1
         	!write(3,*) w(i), w(i)+1, w(i)+2, w(i)+3, w(i)+4
		write(3,*) w(i), 0.0, 0.0, 0.0, 0.0
		write(4,*) w(i), 0.0, 0.0, 0.0, 0.0
      	enddo
	close(3)
	close(4)
	end if !if (benchmark == 1)
	
	!guess self energy or ReadSigma()
	call ReadSigma()
	!do self consistency
	call FindGfBetheDosIntegration()
	
	
end subroutine Initialisation


subroutine FindGfBetheDosIntegration()
	
	! it calculate the green's function of IHM on bethe lattice uisng dos integration it work fine when eta=0.001, nbin=20000, Nenergy=3000
	! arti's eta=0.005, nbin=16000,	
	! created on 19/05/18
	! last modified 19/05/18
	
	!NOTE: Green's function calculated(green) and assigned to corresponding sublattice green_{\alpha}
	
	use Global

        implicit none
        complex*16 gammaA, gammaB, z1
        integer io,i,wo
        real*8 de
        include 'Glob_cons'


	
	!here we are using w grid for Bethe dos but very time consuming
	!if (benchmark == 1) then
	!  open(unit=237,file='bethedos.dat',status='unknown')
	!  do i=-nbin/2,nbin/2-1
	!	if (abs(w(i)) .lt. 2*t ) then
	!		write(237,'(f12.6,e25.12)')w(i), sqrt( (2*t)**2 - w(i)**2)/(2*pi*t**2)
	!	else
	!		write(237,'(f12.6,e25.12)')w(i), sqrt( (2*t)**2 - w(i)**2)/(2*pi*t**2)
	!	end if
	!  end do
        !   close(237)
	!end if !if (benchmark == 1)

	!write the dos for Bethe lattice here you see we useing smaller no of Grid Nenergy=100, compared to "w" grid of size=16000
	de = 4*t/Nenergy
	do wo = 1,Nenergy
		energy(wo)=-2*t + de*(wo-1)
	end do
	if (benchmark == 1) then
		open(unit=237,file='bethedos.dat',status='unknown')
		do wo = 1,Nenergy
	 		write(237,'(f12.6,e25.12)')energy(wo), sqrt((2*t)**2 - energy(wo)**2)/(2*pi*t**2)
		end do
        close(237)
	end if !if (benchmark == 1)
	
	
    if (sub_lattice == 0) then ! sub_lattice == 0 means A sublattice
    	do io = 1,Norbs
        	do i = -nbin/2,nbin/2-1
            	z1 = w(i) + ii*eta
            	gammaA = z1 + mu - delta - sigma_a(io,i)
	   	 		gammaB = z1 + mu + delta - sigma_b(io,i)

            	! do for dos integration
            	green(io,i) = 0.0
	    		!do wo = -nbin/2,nbin/2-1
	    		do wo = 1,Nenergy
					!green_a(io,i)= green_a(io,i) + binsize*gammaB*sqrt((2*t)**2 - w(wo)**2)/((gammaA*gammaB - w(wo)**2)*2*pi*t**2 )
					green(io,i)= green(io,i) + de*gammaB*sqrt((2*t)**2 - energy(wo)**2)/((gammaA*gammaB - energy(wo)**2)*2*pi*t**2 )
            	end do ! w0
				green_a(io,i) = green(io,i) ! assigning to sublattic

        	end do  !i
       	end do !io

	else ! for B sublattice  dos is sqrt(4t^2-e^2)/2*pi*t^2
    	do io=1,Norbs
        	do i=-nbin/2,nbin/2-1
            	z1 = w(i) + ii*eta
            	gammaA = z1 + mu - delta - sigma_a(io,i)
	    		gammaB = z1 + mu + delta - sigma_b(io,i)

            	! do for dos integration
            	green(io,i) = 0.0
	    		!do wo = -nbin/2,nbin/2-1
	    		do wo = 1,Nenergy
					!green_b(io,i)= green_b(io,i) + binsize*gammaA*sqrt((2*t)**2 - w(wo)**2)/((gammaA*gammaB - w(wo)**2)*2*pi*t**2 )
					green(io,i)= green(io,i) + de*gammaA*sqrt((2*t)**2 - energy(wo)**2)/((gammaA*gammaB - energy(wo)**2)*2*pi*t**2 )
					green(io,i) = green_b(io,i)
            	end do ! w0
				green_b(io,i) = green(io,i) ! assigning to sublattic

          	end do !i
       	end do   !io
	end if !if (sub_lattice==0)

	if (benchmark == 1) then
		open(unit=238,file='GfA.out',status='unknown')
	   	open(unit=239,file='GfB.out',status='unknown')
	   	do i=-nbin/2,nbin/2-1
			!write(*,FMT="(5(F20.15,1x))"), w(i), real(green_a(1,i)), imag(green_a(1,i)), real(green_a(2,i)), imag(green_a(2,i))		
			write(238,FMT="(5(F20.15,1x))") w(i), real(green_a(1,i)), imag(green_a(1,i)), real(green_a(2,i)), imag(green_a(2,i))		
			write(239,FMT="(5(F20.15,1x))") w(i), real(green_b(1,i)), imag(green_b(1,i)), real(green_b(2,i)), imag(green_b(2,i))		
	   	end do
       	close(238)
	   	close(239)
	end if !if (benchmark == 1)

	if (sub_lattice == 0) then ! sub_lattice == 0 means A sublattice
		do io = 1,Norbs
			do i = -nbin/2,nbin/2-1
				green0(io,i) = 1.d0 + (sigma_a(io,i)- HFa(io))*green_a(io,i)  
				green0(io,i) = green_a(io,i) / green0(io,i)
				green0_a(io,i) = green0(io,i)
			end do
		end do
	else ! B sublattice calculation
		do io = 1,Norbs
			do i = -nbin/2,nbin/2-1
				green0(io,i) = 1.d0 + (sigma_b(io,i)- HFb(io))*green_b(io,i)  
				green0(io,i) = green_b(io,i) / green0(io,i)
				green0_b(io,i) = green0(io,i)
			end do
		end do
	end if !(sub_lattice == 0)
	
	!calculating host currected green's function
	
	
	!TODO set HFa and HFb
    
     return
end subroutine FindGfBetheDosIntegration


subroutine WriteOutput()
	use Global
   	implicit none
	real*8 r1,r2,r3,r4,r5,r6
	integer info,io, nloop
   	integer i,j
	
	
	open(unit=3,file='GfA.out',status='unknown')
	open(unit=4,file='Sig.out',status='unknown')

	do i=-nbin/2,nbin/2-1
    	write(3,*) w(i), real(green(1,i)), imag(green(1,i)), real(green(2,i)), imag(green(2,i))
		write(4,*) w(i), real(Sigma(1,i)), imag(Sigma(1,i)), real(Sigma(2,i)), imag(Sigma(2,i))
   	enddo
	close(3)
	close(4)	
end subroutine WriteOutput

