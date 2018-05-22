!/home/soumen/INSTALLPACKAGES/mpich2/bin/mpif90 Global_mod.f90 ReadIn.f90 convolution.f90 ksum.f90 Main_bethe.f90
program Bethe

	use Global
	implicit none
	real*8 r1,r2,r3,r4,r5,r6
	integer info,io, nloop
   	integer i,j
	real*8 rho,convA,convB,tolf,lval,Lutt_val,mu_0tmp,nf_tot_tmp
   	real*8 ::  muin(1),t5,t4,tmp_muc
   	real*8,allocatable ::prhoA(:,:),prhoB(:,:)
	logical check
	include 'Glob_cons'
	
	allocate(prhoA(Norbs,-nbin/2:nbin/2-1),stat=info)
	allocate(prhoB(Norbs,-nbin/2:nbin/2-1),stat=info)

	nloop=1

	! set the parameters values by hand
	call ReadParams()
	open(unit=25,file="dataU=0.5",status="new")
	call FindGfBetheDosIntegration() ! !(delta,BetheDos)=>(green_a&b, green0_a&b, na&b, n0a&b) non interating half filled case
	
	do while( nloop.le.MinNoLoop.or.((convA.gt.ToleranceDos).and.(convB.gt.ToleranceDos).and.(nloop.le.MaxNoLoop)) ) 
		
		print*, "loop no-----------------------------------------:", nloop
		do j=1,Norbs                                    !Orbitals
	 		do i=-nbin/2,nbin/2-1
	     		prhoA(j,i)=-dimag(green_a(j,i))/pi ! previous iteration dos, we will be converging this quantity
	   		enddo
		enddo	   
	
		! A sublattice Local Green Function
		sub_lattice = 0
		
		call InitialisingA()	!green=green_b, n=n_a, HF=HFa 
		print*, "1. initialised for A sub"
		!call SaveFile("green_after_A_intialization",green)
		call HostGreen() 		!(HF,mu,delta,sub_lattice,green) => (green0,n0) dpnding on sublat, sgn of delta will be chosn	
		print*, "2. HostGreen for A sub"	
        call SelfEnergy	    	!(n,no,green0)=> Sigma2
		print*, "3. SelfEnergy for A sub"	
		call SaveFile("Sigma2_after_A_intialization",Sigma2)
        call LocalGreen() 		!(green0, Sigma2) => (green,n, HF) 
		print*, "4. local green for A sub"
		call SaveData()			!(if sub_lattice=alpha) green_alpha=green, green0_alpha=green0, n_alpha=n, n0_alpha=n0)
		print*, "4. saving green=> green_a for A sub"
	  	r1=0.d0
	  	r2=0.d0
       	do j=1,Norbs			                      !Orbitals
	     	do i=-nbin/2,nbin/2-1
	       		rho=-dimag(green(j,i))/pi
	       		r1=r1+dabs(rho-prhoA(j,i))*binsize
	       		r2=r2+rho*binsize
	    	 end do
      	enddo                                    ! Orbitals
	  	convA=r1/r2
		
		
		do j=1,Norbs                                    !Orbitals
	 		do i=-nbin/2,nbin/2-1
	     		prhoB(j,i)=-dimag(green_b(j,i))/pi ! previous iteration dos, we will be converging this quantity
	   		enddo
		enddo	

		! B sublattice Local Green Function
		sub_lattice = 1
		call InitialisingB() 	!green=green_a, n=n_b, HF=HFb 
		call HostGreen() 	 	!(HF,mu,delta,sub_lattice,green) => (green0,n0) dpnding on sublat, sgn of delta will be chosn
		call SelfEnergy	     	!(n,no,green0)=> Sigma2
		call LocalGreen()	 	!green0, Sigma2) => (green,n, HF) 
		call Savedata()		 	!(if sub_lattice=alpha) green_alpha=green, green0_alpha=green0, n_alpha=n, n0_alpha=n0)
		
		!save the data to B sublattice G, GO, na nb

	  	r1=0.d0
	  	r2=0.d0
       	do j=1,Norbs			                      !Orbitals
	     	do i=-nbin/2,nbin/2-1
	       		rho=-dimag(green_b(j,i))/pi
	       		r1=r1+dabs(rho-prhoB(j,i))*binsize
	       		r2=r2+rho*binsize
	    	 end do
      	enddo                                    ! Orbitals
	  	convB=r1/r2
         
     
     	write(6,*)&
      	'**************************************************'
		write(6,*) "nloop,convA, convB, na(1), nb(1)"
       	write(6,*) nloop,convA, convB, na(1), nb(1)
		write(6,*)&      
       	'**************************************************'
		call WriteOutput(nloop)
		write(25,FMT="(9(F10.6,1x))") nloop, na(1),na(2),nb(1),nb(2),n0a(1),n0b(2),n0b(1),n0b(2)

       	nloop=nloop+1

   	end do			! END NLOOP

	deallocate(prhoA)
	deallocate(prhoB)
	close(25)
end program Bethe

subroutine SaveData() !(if=>sub_lattice=alpha) green_alpha=green, green0_alpha=green0, n_alpha=n, n0_alpha=n0)	
	use Global
	implicit none
	real*8 r1,r2,r3,r4,r5,r6
	integer info,io, nloop
   	integer i,j
	
	if (sub_lattice == 0) then
		if( benchmark==1) print*, "saving data for A sublattice using SaveData() in Main_bethe.f90"
		do io = 1,Norbs
			do i = -nbin/2,nbin/2-1
				green_a(io,i) 	= green(io,i)
				green0_a(io,i)  = green0(io,i)
				Sigma2a			= Sigma2(io,i)
				na(io)			= n(io)
				n0a(io)			= n0(io)
				HFa(io)		    = HF(io)
				!if( benchmark==1) then
				!	if (i.gt.-10 .and. i.lt.10 ) print*, i,green_a(io,i), green(io,i)
				!end if
			end do
		end do
	else
		if( benchmark==1) print*, "saving data for B sublattice using SaveData() in Main_bethe.f90"
		do io = 1,Norbs
			do i = -nbin/2,nbin/2-1
				green_b(io,i) 	= green(io,i)
				green0_b(io,i)  = green0(io,i)
				Sigma2b			= Sigma2(io,i)
				nb(io)			= n(io)
				n0b(io)			= n0(io)
				HFb(io)		    = HF(io)
				!if( benchmark==1) then
				!	if (i.gt.-10 .and. i.lt.10 ) print*, i,green_b(io,i), green(io,i)
				!end if
			end do
		end do
	end if
		
end subroutine SaveData

subroutine InitialisingA() !green=green_b, n=n_a, HF=HFa
	! call this before calculaitng self energy because selfenergy need n,n0,HF0,Delta=t^2*GB
	use Global
	implicit none
	real*8 r1,r2,r3,r4,r5,r6
	integer info,io, nloop
   	integer i,j
	
	do io = 1,Norbs
		do i = -nbin/2,nbin/2-1
			green(io,i) 	= green_b(io,i)				
			n(io)			= na(io)
			HF(io)		    = HFa(io)
		end do
	end do
end subroutine InitialisingA

subroutine InitialisingB() !!green=green_a, n=n_b, HF=HFb
	! call this before calculaitng self energy because selfenergy need n,n0,HF0, Delta = t^2*GA
	use Global
	implicit none
	real*8 r1,r2,r3,r4,r5,r6
	integer info,io, nloop
   	integer i,j
	
	do io = 1,Norbs
		do i = -nbin/2,nbin/2-1
			green(io,i) 	= green_a(io,i)				
			n(io)			= nb(io)
			HF(io)		    = HFb(io)
		end do
	end do
end subroutine InitialisingB

subroutine SaveFile(FileName,DataSave)
	use Global
	implicit none
	character(len=30),intent(in) :: FileName
	complex*16, intent(in)::	 DataSave(Norbs,-nbin/2:nbin/2-1)
	integer i
	
	open(unit=250,file=FileName,status='unknown')
	do i=-nbin/2,nbin/2-1
		write(250,FMT="(5(F20.15,1x))") w(i), real(DataSave(1,i)), imag(DataSave(1,i)), real(DataSave(2,i)), imag(DataSave(2,i))
	end do
	close(250)
end subroutine SaveFile	

subroutine WriteOutput(loop)
	use Global
	implicit none
	integer,intent(in) :: loop
	integer i
	character(len=1024) :: filename
    character(len=1024) :: format_string1,format_string2
   

    
        if (loop < 10) then
            format_string2 = "(A8,I1)" ! SigA.out,Gf0A.out
			format_string1 = "(A7,I1)" ! GfA.out
        else
            format_string2 = "(A8,I2)" ! SigA.out,Gf0A.out
			format_string1 = "(A7,I2)" ! GfA.out
        endif

    write (filename,format_string1) "GfA.out", loop
	open(unit=250,file=filename,status='unknown')

	write (filename,format_string1) "GfB.out", loop
	open(unit=251,file=filename,status='unknown')

	write (filename,format_string2) "Gf0A.out", loop
	open(unit=252,file=filename,status='unknown')

	write (filename,format_string2) "Gf0B.out", loop
	open(unit=253,file=filename,status='unknown')

	write (filename,format_string2) "SigA.out", loop
	open(unit=254,file=filename,status='unknown')

	write (filename,format_string2) "SigB.out", loop
	open(unit=255,file=filename,status='unknown')

	do i=-nbin/2,nbin/2-1
		write(250,FMT="(5(F20.15,1x))") w(i), real(green_a(1,i)), imag(green_a(1,i)), real(green_a(2,i)), imag(green_a(2,i))
		write(251,FMT="(5(F20.15,1x))") w(i), real(green_b(1,i)), imag(green_b(1,i)), real(green_b(2,i)), imag(green_b(2,i))
		write(252,FMT="(5(F20.15,1x))") w(i), real(green0_a(1,i)), imag(green0_a(1,i)), real(green0_a(2,i)), imag(green0_a(2,i))
		write(253,FMT="(5(F20.15,1x))") w(i), real(green0_b(1,i)), imag(green0_b(1,i)), real(green0_b(2,i)), imag(green0_b(2,i))
		write(254,FMT="(5(F20.15,1x))") w(i), real(Sigma2a(1,i)), imag(Sigma2a(1,i)), real(Sigma2a(2,i)), imag(Sigma2a(2,i))
		write(255,FMT="(5(F20.15,1x))") w(i), real(Sigma2b(1,i)), imag(Sigma2b(1,i)), real(Sigma2b(2,i)), imag(Sigma2b(2,i))
	end do
	close(250)
	close(251)
	close(252)
	close(253)
	close(254)
	close(255)

end subroutine WriteOutput
