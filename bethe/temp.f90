
	
	subroutine FindGfBethe()
	
	! it calculate the green's function of IHM on bethe lattice
	use Global

        implicit none
        complex*16 gammaA, gammaB, z1
        integer io,i,jo,wo,pw,ie,info
        real*8 m2,m4,r1,r2norm,de
        include 'Glob_cons'

	if (benchmark == 1) then
	   open(unit=237,file='bethedos.dat',status='unknown')
	   do i=-nbin/2,nbin/2-1
		if (abs(w(i)) .lt. 2*t ) then
			write(237,'(f12.6,e25.12)')w(i), sqrt( (2*t)**2 - w(i)**2)/(2*pi*t**2)
		else
			write(237,'(f12.6,e25.12)')w(i), sqrt( (2*t)**2 - w(i)**2)/(2*pi*t**2)
		end if
	   end do
           close(237)
	end if !if (benchmark == 1)
	
	
        if (sub_lattice == 0) then ! sub_lattice == 0 means A sublattice
         do io = 1,Norbs
          do i = -nbin/2,nbin/2-1
            z1 = w(i) + ii*eta
            gammaA = z1 + mu -delta - sigma_a(io,i)
	    gammaB = z1 + mu + delta -sigma_b(io,i)

            ! do for dos integration
            green_a(io,i) = 0.0
	    do wo = -nbin/2,nbin/2-1
			green_a(io,i)= green_a(io,i) + binsize*gammaB*sqrt((2*t)**2 - w(w0)**2)/((gammaA*gammaB - w(w0)**2)*2*pi*t**2 )
            end do ! w0

          end do  !i
         end do !io

	 else ! for B sublattice  dos is sqrt(4t^2-e^2)/2*pi*t^2
         do io=1,Norbs
          do i=-nbin/2,nbin/2-1
            z1 = w(i) + ii*eta
            gammaA = z1 + mu -delta - sigma_a(io,i)
	    gammaB = z1 + mu + delta -sigma_b(io,i)

            ! do for dos integration
            green_b(io,i) = 0.0
	    do wo = -nbin/2,nbin/2-1
			green_b(io,i)= green_b(io,i) + binsize*gammaA*sqrt((2*t)**2 - w(w0)**2)/((gammaA*gammaB - w(w0)**2)*2*pi*t**2 )
            end do ! w0

          end do !i
        end do   !io
       
       end if !if (sub_lattice==0)
 
     
        
       do io=1,Norbs
        do i=-N,N
        Gfscript(io,i)=(1.d0+(sigma(io,i)+ep_f(io)-mu_c+mu0)*Gf(io,i))
        Gfscript(io,i)=Gf(io,i)/Gfscript(io,i)
        if(dimag(Gfscript(io,i)).gt.0.d0)  &
            Gfscript(io,i)=dcmplx(dreal(Gfscript(io,i)),-1.D-9)
        end do !i
      end do !io

      return
      end
	
