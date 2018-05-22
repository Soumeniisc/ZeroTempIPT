subroutine SelfEnergy !(n,no,green0)=> Sigma2
	
	!It takes Hartee corrected Host Green's function(\mathcal{G}_0), "n" and "n_0"(later two needed only for A facotr calculation)
	!It calculate Hartee corrected SelfEnergy with "A" factor. It uses other two external subroutines EvalualteChi and convolution.
	!  
  	use Global
	implicit none
	integer io, jo, i, info
	real*8 :: A(1:Norbs)

	real*8,allocatable :: Rho1(:),Rho2(:),RhoSigma1(:),RhoSigma2(:),RhoSigma(:),ReSigma(:)
	allocate(Rho1(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of Rho1 failed'
    allocate(Rho2(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of Rho2 failed'
    allocate(RhoSigma1(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of Rhosigma1 failed'
	allocate(RhoSigma2(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of Rhosigma2 failed'
	allocate(RhoSigma(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of Rhosigma failed'
	allocate(ReSigma(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of Resigma failed'
	
	
	call EvaluateChi ! evaluate chi1 and chi2

   	Sigma2=zero
  	do io=1,Norbs    !1, 2

     	do jo=1,Norbs  !2, 1

           	if(jo.ne.io) then !here spin mixing happen
				
				!TODO check spin mixing is correct or not. Do analytical calculation for that
            	do i=-nbin/2,nbin/2-1
	      			Rho1(i)=(-dimag(green0(io,-i))/pi)*fermi(i)
	      			Rho2(i)=chi1(jo,i)
	    		end do
  	    		call convolution(Rho1,Rho2,RhoSigma1) 
   
             	do i=-nbin/2,nbin/2-1
	     			Rho1(i)=(-dimag(green0(io,-i))/pi)*fermi(-i)
	     			Rho2(i)=chi2(jo,i)
	    		end do
	    		call convolution(Rho1,Rho2,RhoSigma2)
             	do i=-nbin/2,nbin/2-1
	     			RhoSigma(i)=RhoSigma1(i) + RhoSigma2(i)
	   			end do
            	!write(6,*) 'CONVOLUTION DONE, NOW KKTRANSF'
				
				if (kktransformation==1 ) then 
					print*, "kktransf1 is used to get real part of self energy"
            		call kktransf(RhoSigma,ReSigma) ! my code giving infinty at two end points
				else
					print*, "kktransf2 is used to get real part of self energy"
					call kktransf2(RhoSigma,ReSigma) ! used by dasari
				end if
            


            	do i = -nbin/2,nbin/2-1
             		Sigma2(io,i)= U**2 * (ReSigma(i)-pi*ii*RhoSigma(i))   
           		end do

           	end if !if(jo.ne.io) then

     	end do ! End orbital index jo

	end do   ! End orbital index io
        
	!   Now calculate the n, A and B factors
	
	do io=1,Norbs    !1, 2
     	do jo=1,Norbs  !2, 1
           	if(jo.ne.io) then 
				A(io) = n(jo)*(1.d0 - n(jo))/(n0(jo)*(1.d0 - n0(jo)))
	 		end if !if(jo.ne.io) then
     	end do ! End orbital index jo
	end do   ! End orbital index io
	
	

	do io = 1, Norbs
		do i=-nbin/2,nbin/2-1
			Sigma2(io,i) = A(io)*Sigma2(io,i)
		end do 
	end do
	
	! Print Hartee corrected Sigma
	if (HarteeCorrectedSig==1) then
		open(unit=3,file='HarteeCorrectedSig.out',status='unknown')
		do i=-nbin/2,nbin/2-1
			write(3,*) w(i), Sigma2(io,i)
		  	enddo
		close(3)
	end if

end subroutine SelfEnergy

subroutine EvaluateChi

    use Global
	real*8 r1,r2,r3,r4,r5,r6
    integer i,j,info
	real*8,allocatable :: rho1(:),rho2(:),chi(:)
	include 'Glob_cons'

    allocate(rho1(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of rho1 failed'
    allocate(rho2(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of rho2 failed'
    allocate(chi(-nbin/2:nbin/2-1),stat=info)
    if(info.ne.0) write(6,*) 'Allocation of chi failed'



    do j=1,Norbs
		do i=-nbin/2,nbin/2-1
	    	rho1(i)=-fermi(i)*dimag(green0(j,i))/pi
	    	rho2(i)=-fermi(-i)*dimag(green0(j,i))/pi
		end do
        call convolution(rho1,rho2,chi)
        do i=-nbin/2,nbin/2-1
           	chi1(j,i)=chi(i)
        end do
  	enddo !j=1,Norbs

    !if(task_id==0) write(6,*) 'Chi1 done'
  	do j=1,Norbs
	  	do i=-nbin/2,nbin/2-1
	    	rho1(i)=-fermi(-i)*dimag(green0(j,i))/pi
	    	rho2(i)=-fermi(i)*dimag(green0(j,i))/pi
	  	end do
        call convolution(rho1,rho2,chi)
        !if(task_id==0) write(6,*) 'Chi2 done,orb=',j
       	do i=-nbin/2,nbin/2-1
        	chi2(j,i)=chi(i)
       	end do
 	enddo

    return
end subroutine EvaluateChi

subroutine convolution(rho1,rho2,rho3)
!			      /	
!	Calculate   rho3(w) = | dw' rho1(w') rho2(w+w')
!			      /
	! check Arti's code what she used. this part of the code need regorous benchmarking https://en.wikipedia.org/wiki/Convolution

	use Global
	real*8 r1
	integer i,j
    integer jmin,jmax
    real*8 :: rho1(-nbin/2:nbin/2-1),rho2(-nbin/2:nbin/2-1), rho3(-nbin/2:nbin/2-1)
       

    do i = -nbin/2,nbin/2-1
      	jmin = max(-nbin/2-i, -nbin/2)
        jmax = min(nbin/2-i-1,nbin/2-1)
        r1 = 0.d0
       	do j = jmin,jmax
           	ji = j+i 
            r1=r1+rho1(j)*rho2(ji)*binsize
       	end do
        rho3(i)=r1
  	end do !i = -nbin/2,nbin/2-1

	if(benchmark==2) then
		print*, "jmin and jmax in convolution(f(j),f(i+j)) for nbin=20"
		do i = -10,10-1 ! here i have manualy considered nbin/2 = 10, one cane check jmin and jmax value
      		jmin = max(-10-i, -10)
        	jmax = min(10-i-1,10-1)
			print*, "i:  ",i,"jmin:  ", jmin, "jmax:  ", jmax
		end do
	end if !if(benchmark==1)

end subroutine convolution

subroutine kktransf(RhoSigma_,ReSigma_)
	use Global
    implicit none
	real*8 r1,r2,r3,r4,r5,r6
    integer i,j,info
	integer llim,ulim,Ntmp

	!real*8 :: RhoSigma(-nbin/2:nbin/2-1),rlsg(-Nm:Nm+1), rlsigma(-nbin/2:nbin/2-1),wo(-nbin/2:nbin/2-1)
	real*8 :: RhoSigma_(-nbin/2:nbin/2-1), ReSigma_(-nbin/2:nbin/2-1)
    real*8 :: woj,dspj,resgi,spi,dspi,dwi,woji1,woji2

!		      oo		
!		       /			
!	rlsigma(w)= -P | dw' rhosigma(w')
!		       /     -----------	
!	             -oo       w' - w


!		KRAMERS-KRONIG, Raja used different algorith

	do i = -nbin/2,nbin/2-1
		resgi = 0.d0
		do j = -nbin/2,nbin/2-1
			if (i .ne. j) then
				resgi = resgi - (RhoSigma_(i)-RhoSigma_(j))/(w(i) - w(j))*binsize
			else 
				resgi = resgi - RhoSigma_(j)*dlog((w(nbin/2-1) - w(j))/(w(j)-w(-nbin/2)))
			end if
		end do
		ReSigma_(i) = resgi
	end do
	
	
	 return
	 end

subroutine kktransf2(RhoSigma_,ReSigma_)

	use Global
    implicit none
	real*8 r1,r2,r3,r4,r5,r6
    integer i,j,info
	integer llim,ulim,Ntmp

	!real*8 :: RhoSigma(-nbin/2:nbin/2-1),rlsg(-Nm:Nm+1), rlsigma(-nbin/2:nbin/2-1),wo(-nbin/2:nbin/2-1)
	real*8 :: RhoSigma_(-nbin/2:nbin/2-1), rlsg(-nbin/2:nbin/2), ReSigma_(-nbin/2:nbin/2-1), wo(-nbin/2:nbin/2-1)
    real*8 :: woj,dspj,resgj,spi,dspi,dwi,woji1,woji2

!		      oo		
!		       /			
!	rlsigma(w)= -P | dw' rhosigma(w')
!		       /     -----------	
!	             -oo       w' - w


!		KRAMERS-KRONIG, algorith used by Dasari

	llim = -nbin/2 + 1
	ulim = nbin/2 - 1		! calculate ReSig(w) for all w

	do j = llim,ulim

!	interlacing grid

		wo(j) = .5D0*(w(j-1)+w(j))
		woj   = wo(j)
		dspj  = RhoSigma_(j) - Rhosigma_(j-1)
		resgj = 0.0D0
	
		do i = -nbin/2,j-2
	  		spi   = RhoSigma_(i)
	  		dspi  = RhoSigma_(i+1) - RhoSigma_(i)
	  		dwi   = w(i+1) - w(i)
	  		woji1 = w(i) - woj
	  		woji2 = w(i+1) - woj
	  		r1    = dlog(woji2/woji1)

	  		resgj = resgj-(spi*r1 + dspi )
	  		resgj = resgj-(dspi/dwi)*(woj -w(i))*r1
	  	end do

!	 skip the interval (j-1) to j
 	
		do i = j,nbin/2-1-1

	  		spi   = RhoSigma_(i)
	  		dspi  = RhoSigma_(i+1) - RhoSigma_(i)
	  		dwi   = w(i+1)-w(i)
	  		woji1 = w(i)-woj
	  		woji2 = w(i+1)-woj
	  		r1    = dlog(woji2/woji1)

	  		resgj=resgj-(spi*r1 + dspi )
	  		resgj=resgj-(dspi/dwi)*(woj -w(i))*r1

	  	end do ! w(i) integration end here  for real part at w(j). here i run from -nbin/2 to nbin/2-2

	  	resgj=resgj - dspj
	  	rlsg(j)=resgj

	end do

	rlsg(ulim+1)=rlsg(ulim) ! here might come memorry alocation error.
	rlsg(llim-1)=rlsg(llim)
	
	do i=llim-1,ulim
	ReSigma_(i)=0.5d0*(rlsg(i)+rlsg(i+1))
	end do


	return
end subroutine kktransf2



subroutine TestKKTransformation()
	use Global
    implicit none
	integer i
	real*8 :: RhoSigma_(-nbin/2:nbin/2-1), ReSigma_(-nbin/2:nbin/2-1)

	if (sub_lattice==0) then
		print*, "calculating real part of green's function from up spin dos of A sublattice. look at GfA.out.Kramers"
		do i = -nbin/2,nbin/2-1
			RhoSigma_(i) = -imag(green_a(1,i))/pi
		end do

		if(kktransformation==1) then
			call kktransf(RhoSigma_,ReSigma_)
			open(unit=3,file='GfA.out.Kramers',status='unknown')
		else
			call kktransf2(RhoSigma_,ReSigma_)
			open(unit=3,file='GfA.out.Kramers2',status='unknown')
		end if

		do i=-nbin/2,nbin/2-1
			write(3,*) w(i), ReSigma_(i), RhoSigma_(i)		
   		enddo
		close(3)

	else

		print*, "calculating real part of green's function from up spin dos of B sublattice. look at GfB.out.Kramers"
		do i = -nbin/2,nbin/2-1
			RhoSigma_(i) = -imag(green_b(1,i))/pi
		end do
	
		if(kktransformation==1) then
			call kktransf(RhoSigma_,ReSigma_)
			open(unit=3,file='GfB.out.Kramers',status='unknown')
		else
			call kktransf2(RhoSigma_,ReSigma_)
			open(unit=3,file='GfB.out.Kramers2',status='unknown')
		end if

		do i=-nbin/2,nbin/2-1
			write(3,*) w(i), ReSigma_(i), RhoSigma_(i)		
   		enddo
		close(3)

	end if
	
	

end subroutine TestKKTransformation
	
	

	




