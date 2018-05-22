module Global
	implicit none
	save
   	integer nbin, Nenergy
   	parameter(nbin=20000)
   	parameter(Nenergy=3000)

   	double precision 	Wmax, Wmin, binsize
   	double precision	 w(-nbin/2:nbin/2-1),energy(1:Nenergy), fermi(-nbin/2:nbin/2-1)
	
	!Hamiltonian Parameters
	double precision mu, delta, t, t2, U, eta
	integer    sub_lattice, Norbs
	parameter(Norbs=2)

	!Global constants
	complex*16::ii,zero
	real*8  :: pi

	! iterations and convergence
	integer MinNoLoop, MaxNoLoop
	double precision ToleranceDos

	! switches to on and off firent functionality
	integer benchmark           ! whether you want to bechmark code
	integer kktransformation    ! = 1 for staightforward calculation, =2 for little more corrected calculation 
	integer HarteeCorrectedSig  ! whether to print HarteeCorrectedSig in SelfEnergy Subroutine

	!obervables
	double precision n(1:Norbs), n0(1:Norbs), HF(1:Norbs), na(1:Norbs), nb(1:Norbs), n0a(1:Norbs), n0b(1:Norbs), HFa(1:Norbs), HFb(1:Norbs)


   	!Green's function and self energy
	complex*16 green0(Norbs,-nbin/2:nbin/2-1), sigma0(Norbs,-nbin/2:nbin/2-1)
	complex*16 green0_a(Norbs,-nbin/2:nbin/2-1), sigma0_a(Norbs,-nbin/2:nbin/2-1)
	complex*16 green0_b(Norbs,-nbin/2:nbin/2-1), sigma0_b(Norbs,-nbin/2:nbin/2-1)
	complex*16 green(Norbs,-nbin/2:nbin/2-1)
	complex*16 green_a(Norbs,-nbin/2:nbin/2-1), sigma_a(Norbs,-nbin/2:nbin/2-1)
 	complex*16 green_b(Norbs,-nbin/2:nbin/2-1), sigma_b(Norbs,-nbin/2:nbin/2-1)
	complex*16 chi1(Norbs,-nbin/2:nbin/2-1), chi2(Norbs,-nbin/2:nbin/2-1)
	complex*16 Sigma2(Norbs,-nbin/2:nbin/2-1), Sigma2a(Norbs,-nbin/2:nbin/2-1), Sigma2b(Norbs,-nbin/2:nbin/2-1)! Hartee self energy

end module Global
