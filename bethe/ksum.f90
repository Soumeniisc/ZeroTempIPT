!TODO FindGfBetheAnalytic() from analytic expression of Gf for Bethe lattice Dos integration 
!TODO make a clear understanding about eta factor , dos energy grid(de), green's function grid size(binsize).

subroutine FindGfBetheDosIntegration() !(delta,BetheDos)=>(green_a&b, green0_a&b, na&b, n0a&b) non interating half filled case
	
	! it calculate the noninteracting green's function of IHM at B site on bethe lattice uisng dos integration it work fine when eta=0.001, nbin=20000, Nenergy=3000
	! arti's eta=0.005, nbin=16000,	
	! created on 19/05/18
	! last modified 21/05/18
	
	!NOTE: Green's function calculated(green) and assigned to corresponding sublattice green_{\alpha}
	
	use Global

        implicit none
        complex*16 gammaA, gammaB, z1
        integer io,i,wo
        real*8 de
		complex*16:: G
        include 'Glob_cons'

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
	
	
	
    sub_lattice =0  ! sub_lattice == 0 means A sublattice
    	do io = 1,Norbs
        	do i = -nbin/2,nbin/2-1
            	z1 = w(i) + ii*eta
            	gammaA = z1  - delta !- sigma_a(io,i)
	   	 		gammaB = z1  + delta !- sigma_b(io,i)

            	! do for dos integration
            	green(io,i) = 0.0
	    		!do wo = -nbin/2,nbin/2-1
	    		do wo = 1,Nenergy
					!green_a(io,i)= green_a(io,i) + binsize*gammaB*sqrt((2*t)**2 - w(wo)**2)/((gammaA*gammaB - w(wo)**2)*2*pi*t**2 )
					green(io,i)= green(io,i) + de*gammaB*sqrt((2*t)**2 - energy(wo)**2)/((gammaA*gammaB - energy(wo)**2)*2*pi*t**2 )
            	end do ! w0
				green_a(io,i) = green(io,i) ! assigning to sublattic
        	end do  !i
			HF(io) = 0.d0
       	end do !io

	!calculating Host green function for B sublattice
	do io = 1, Norbs
		do i=-nbin/2,nbin/2-1	
			z1 = w(i) + ii*eta	
			G = z1  + delta - t**2*green_a(io,i)			
			green0(io,i) = 1.d0/G
			green0_b(io,i) = green0(io,i)
		end do 
	end do

	!call HostGreen()   !(HF,mu,delta,sub_lattice,green) => (green0,n0) dpnding on sublat, sgn of delta will be chosn
	call occupancies() !green => n
	call occupancies0() !green0 => n0
	print*,"non interacting A sublattice n(1), n(2):", n(1), n(2)
	do io = 1, Norbs
		na(io) = n(io) 
		n0a(io) = n0(io) 
	end do


	sub_lattice =1 ! for B sublattice  
    	do io=1,Norbs
        	do i=-nbin/2,nbin/2-1
            	z1 = w(i) + ii*eta
            	gammaA = z1  - delta !- sigma_a(io,i)
	    		gammaB = z1  + delta !- sigma_b(io,i)

            	! do for dos integration
            	green(io,i) = 0.0
	    		!do wo = -nbin/2,nbin/2-1
	    		do wo = 1,Nenergy
					!green_b(io,i)= green_b(io,i) + binsize*gammaA*sqrt((2*t)**2 - w(wo)**2)/((gammaA*gammaB - w(wo)**2)*2*pi*t**2 )
					green(io,i)= green(io,i) + de*gammaA*sqrt((2*t)**2 - energy(wo)**2)/((gammaA*gammaB - energy(wo)**2)*2*pi*t**2 )
					!green(io,i) = green_b(io,i)
            	end do ! w0
				green_b(io,i) = green(io,i) ! assigning to sublattic
				HF(io) = 0.d0
          	end do !i
       	end do   !io

	!calculating Host green function for A sublattice  w - delta - t^2GB
	do io = 1, Norbs
		do i=-nbin/2,nbin/2-1
			z1 = w(i) + ii*eta
			G = z1  - delta- t**2*green_b(io,i)
			green0(io,i) = 1.d0/G
			green0_a(io,i) = green0(io,i)
		end do 
	end do

	!call HostGreen()   !(HF,mu,delta,sub_lattice,green) => (green0,n0) dpnding on sublat, sgn of delta will be chosn
	call occupancies() !green => n
	call occupancies0() !green0 => n0
	print*,"non interacting B sublattice n(1), n(2):", n(1), n(2)
	do io = 1, Norbs
		nb(io) = n(io) 
		n0b(io) = n0(io) 
	end do

	if (benchmark == 1) then
		open(unit=238,file='GfA.out.noninteracting',status='unknown')
	   	open(unit=239,file='GfB.out.noninteracting',status='unknown')
		open(unit=240,file='Gf0A.out.noninteracting',status='unknown')
	   	open(unit=241,file='Gf0B.out.noninteracting',status='unknown')
	   	do i=-nbin/2,nbin/2-1
			!write(*,FMT="(5(F20.15,1x))"), w(i), real(green_a(1,i)), imag(green_a(1,i)), real(green_a(2,i)), imag(green_a(2,i))		
			write(238,FMT="(5(F20.15,1x))") w(i), real(green_a(1,i)), imag(green_a(1,i)), real(green_a(2,i)), imag(green_a(2,i))		
			write(239,FMT="(5(F20.15,1x))") w(i), real(green_b(1,i)), imag(green_b(1,i)), real(green_b(2,i)), imag(green_b(2,i))	
			write(240,FMT="(5(F20.15,1x))") w(i), real(green0_a(1,i)), imag(green0_a(1,i)), real(green0_a(2,i)), imag(green0_a(2,i))		
			write(241,FMT="(5(F20.15,1x))") w(i), real(green0_b(1,i)), imag(green0_b(1,i)), real(green0_b(2,i)), imag(green0_b(2,i))	
	   	end do
       	close(238)
	   	close(239)
	end if !if (benchmark == 1)

	 
	
     return
end subroutine FindGfBetheDosIntegration

subroutine LocalGreen() !(green0, Sigma2) => (green,n HF) 
	use Global
	implicit none
	integer io, jo, i, info
	complex*16:: G

	do io = 1, Norbs
		do i=-nbin/2,nbin/2-1
			G = 1.d0 - Sigma2(io,i)*green0(io,i)
			green(io,i) = green0(io,i)/G
		end do 
	end do
	
	call occupancies()

	call Hartee()
end subroutine LocalGreen

subroutine HostGreen() !(HF,mu,delta,sub_lattice,green)=> (green0,n0)

	use Global
	implicit none
	integer io, jo, i, info
	complex*16:: G
	
	if (sub_lattice==0) then
		do io = 1, Norbs
			do i=-nbin/2,nbin/2-1
				G = w(i)  - delta - HF(io) - t**2*green(io,i)
				green0(io,i) = 1.d0/G
				!write(250,FMT="(5(F20.15,1x))") w(i), real(green0(1,i)), imag(green0(1,i)), real(green0(2,i)), imag(green0(2,i))
			end do 
		end do

		if(benchmark==1) then
			open(unit=250,file='Gf0A.out.loop',status='unknown')
			do i=-nbin/2,nbin/2-1
				write(250,FMT="(5(F20.15,1x))") w(i), real(green0(1,i)), imag(green0(1,i)), real(green0(2,i)), imag(green0(2,i))
			end do
		close(250)
		end if
	 
		
	else
		
		do io = 1, Norbs
			do i=-nbin/2,nbin/2-1
				G = w(i) + mu + delta - HF(io) - t**2*green(io,i)
				green0(io,i) = 1.d0/G
				write(251,FMT="(5(F20.15,1x))") w(i), real(green0(1,i)), imag(green0(1,i)), real(green0(2,i)), imag(green0(2,i))
			end do 
		end do
		if(benchmark==1) then
			open(unit=251,file='Gf0B.out.loop',status='unknown')
			do i=-nbin/2,nbin/2-1
				write(251,FMT="(5(F20.15,1x))") w(i), real(green0(1,i)), imag(green0(1,i)), real(green0(2,i)), imag(green0(2,i))
			end do
		close(251)
		end if
	end if
	
	call occupancies0()
end subroutine HostGreen


subroutine occupancies() !green => n

	!when delta = 1.0, t=1.0, U = 0.0 na = 0.116382562203721 nb = 0.883617437796280 which agree with the our earlier calculation

	!NOTE: occupancies (n, n0) calculated from green and green0 but not yet assigned to corresponding sublattice n_{\alpha} and n0_{\alpha}. 
	use Global
	implicit none
 	real*8 r1,r2
	integer i,j,io,jo,ie,flag_occ
	include 'Glob_cons'

	do jo = 1,Norbs ! Orbitals
    	r2 = 0.d0
      	r1 = 0.d0
	 	do ie = -nbin/2,nbin/2-1
            if(w(ie) .lt. 0.0) r2 = r2 + binsize*(-dimag(green(jo,ie))/pi)
            r1 = r1 + binsize*(-dimag(green(jo,ie))/pi)
	    end do
        n(jo) = (r2/r1)



        if(dabs(r1-1.d0).gt.0.02d0) then
        	!if(task_id==0) write(6,*) 'WARNING: Gf not normalized-',jo,r1
			write(6,*) 'WARNING: Gf not normalized-',jo,r1
       	end if	 
    	
	enddo   ! Orbitals
end subroutine occupancies

subroutine occupancies0() !green0 => n0

	!when delta = 1.0, t=1.0, U = 0.0 na = 0.116382562203721 nb = 0.883617437796280 which agree with the our earlier calculation

	!NOTE: occupancies (n, n0) calculated from green and green0 but not yet assigned to corresponding sublattice n_{alpha} and n0_{alpha}.
 
	use Global
	implicit none
 	real*8 r1,r2
	integer i,j,io,jo,ie,flag_occ
	include 'Glob_cons'

	do jo = 1,Norbs
	    r1 = 0.d0
		r2 = 0.d0
    	do ie = -nbin/2,nbin/2-1
	    	if(w(ie) .lt. 0.0) r2 = r2 + binsize*(-dimag(green0(jo,ie))/pi)   
	     	r1 = r1 + binsize*(-dimag(green0(jo,ie))/pi)
       	end do
       	n0(jo) = (r2/r1)
      	if(dabs(r1-1.d0).gt.0.02d0) then
       		!if(task_id==0) write(6,*) 'WARNING: G0 not normalized-',jo,r2
			write(6,*) 'WARNING: G0 not normalized-',jo,r2

      	end if	     	
	end do
   
end subroutine occupancies0    



subroutine Hartee()

	use Global
	implicit none
	integer io, jo
	

	do io=1,Norbs    !1, 2
     	do jo=1,Norbs  !2, 1
           	if(jo.ne.io) then 
				HF(io) = U*n(jo)
	 		end if !if(jo.ne.io) then
     	end do ! End orbital index jo
	end do   ! End orbital index io
	
	
end subroutine Hartee

