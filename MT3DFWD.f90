!=======================================================================
!   ifort -I/opt/intel/mkl/include MARE3DCSEM.f90 -lmkl_rt -L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 -o mt3d
!   GENERAL REMARKS
!
!   This code is freely available under the following conditions:
!
!   1) The code is to be used only for non-commercial purposes.
!   2) No changes and modifications to the code without prior permission of the developer.
!   3) No forwarding the code to a third party without prior permission of the developer.
!
!   This code is an finite element modeling program with unstructured mesh for solving the 3D MT problem in anisotropic media using A-phi potentials.
!
!   Dr. Yixin Ye
!   School of Geophysics   
!   East China University of Technology 
!   418 Guanglan Road, nanchang 330013, China 
!   yxye@ecit.cn
!  
!======================================================================! 
!======================================================================!
!======================================================================! 
!======================================================================! 
!                    MT3Dani
! 
! Adaptive Finite Element Modeling based on Unstructured Meshes
!
!           Finite Element Forward Code           
!           For 3D MT Methods in anisotropic media            
!
!   MT Forward Code: Written by Yixin Ye (yixinye321@126.com)    
!
!    Version 1.1 April 1, 2020 
!  
!======================================================================! 
!
!  Acknowledgments:
!       
!  The unstructured tetrahedral meshes used in the FE method were created 
!  using the TetGen mesh generator (Si, 2007)
!
!  The resulting sparse linear system of equations is solved by using the 
!  symmetric successive over-relaxation preconditioned conjugate gradient 
!  (SSOR-PCG) method
!
!======================================================================! 
!==========================================================! FEM_PARAMS
!======================================================================! 
 
  	module  fem_params

   	integer,         parameter :: double = 8
!
! Constants:
!
	real(double),    parameter :: PI  = 3.141592653589793d0
	real(double),    parameter :: PI2  = 2.d0*pi      
	real(double),    parameter :: MU0 = PI * 4.d-7
	complex(double), parameter :: IC  = (0.d0,1.d0)   
!
! QMR solver parameters:
!
	real(double), parameter    :: qmrmaxerr  = 1.d-15   ! reduction in residual
	integer, parameter         :: qmrmaxiter = 5000     ! # qmr iterations
!
! Adaptive refinement error estimator parameters
!
	real(double), parameter    :: qhmaxerr     = 1.d-35 ! L2 projection QMR solve reduction in residual
	integer, parameter         :: qhmaxiter    = 5000   ! max qmr iterations for L2 projection
        real(double),parameter     :: dualqmrmxerr = 1.d-8  !QMR max error for dual solution       
!
!  Convergence Test field minimum field strength (in wavenumber domain)
!  
	real(double),    parameter :: ecutoff = 1d-20
	real(double),    parameter :: hcutoff = 1d-18
!
! Minimum area limit for grid triangles
!
        real(double), parameter    :: sknlimit   = 5.
        real(8), parameter         :: minelearea = 1. ! minimum area in m^2, to ensure no run away refinement near singularities
     
	end module fem_params
!
!======================================================================!
!=============================================================! DATINFO
!======================================================================!
!
	module  datinfo
	use     fem_params
!
! Model and site file in RUNFILE:  
!
	character*50               :: elementfile,nodefile  ! name of element file , name of node file 
!
! Name strings used for the element file, node file:
!
	character*50               :: modelroot,adaptroot,currentname,lastname
!
! TetGen program call string:
!
	character(256)             :: tristr  
! 
! Adaptive mesh options in RUNFILE
!
	integer                    :: strtmodnum  ! starting model number 
	integer                    :: maxadapt    ! vector of 1 or 0 for each frequency  
	real(double)               :: pctrefine
	integer                    :: nsmiter  = 2  ! number of smoothing operations to use
	integer,parameter          :: weighterr = 1 ! 1=yes, 0=no to "weight error" using duality 
        integer,parameter          :: errortype = 3 ! both !E,H,EH: 1,2,3.  field to use for error estimate
        real(double)               :: tolerance     ! tolerance for adaptive solution convergence
        integer,DIMENSION(:)       :: adapt         ! vector of 1 or 0 for each frequency  
	allocatable                :: adapt  
!
! Frequencies to use for refinement, listed in RUNFILE: 
!
        integer                    :: nfrqrefine  ! number of frequencies to consider during refinement interation
        real(8), dimension(:)      :: frqrefine   ! vector of frequencies to consider during refinement interation
        allocatable                :: frqrefine	
!        
! Frequencies to compute solution in RUNFILE
!
	integer                    :: nfreq    ! number of requencies to compute solution
	real(double), dimension(:) :: freqvec  ! frequency vector
	allocatable                :: freqvec
!
! Site location info in RUNFILE:
!	 
        integer		           :: nsite  ! number of sites
	real(double), dimension(:) :: xsite,ysite,zsite,siteslope  ! x- y- z- coordinates of receiver sites
	allocatable                :: xsite,ysite,zsite,siteslope           
!
!   Resistivity block parameters:
!
        INTEGER                    :: NRO    ! number of resistivity blocks
        REAL(double), DIMENSION(:) :: roxx,royy,rozz,alphx,alphy,alphz  
        allocatable                :: roxx,royy,rozz,alphx,alphy,alphz

        REAL(double), DIMENSION(:,:) :: SG  ! conductivity, resistivity
        allocatable                  :: sg 
	
	end module datinfo
!      
!======================================================================!
!================================================================! MESH
!======================================================================!
! This is the mesh info, this stuff is the same for each call to Mt3Dfw
!	
	module  mesh	
	use fem_params
!
! Used by: input									     !
!
	integer		           :: nnod,nele,iskip,e
	integer                    :: nnodlast,skip  ! size of "previous" grid
	integer, dimension(:)      :: nodlab,irg, nelein	
	integer, dimension(:,:)    :: emap, nodele, indinod1
	real(double), dimension(:) :: xn,yn,zn         
	integer,      dimension(:) :: bcnod
	
	allocatable                :: bcnod,nodele,nelein,indinod1
	allocatable                :: nodlab,emap,xn,yn,zn,irg
	
	end module mesh
!      
!======================================================================! 
!===============================================================! FEMOUT
!======================================================================! 
! Output solutions of FEM modeling
!
	module          femout
	use             fem_params
!	
! Total fields at (isites,ifreq):
!
	complex(double),dimension(:,:)   :: ex1, ey1, ez1, hx1, hy1, hz1
	complex(double),dimension(:,:)   :: ex2, ey2, ez2, hx2, hy2, hz2
        complex(double),dimension(:,:)   :: zxy, zyx
        real(double),dimension(:,:)      :: rhoxy, rhoyx, phaxy, phayx
!
	allocatable   :: ex1,ey1,ez1,hx1,hy1,hz1 ! deallocated in driver
	allocatable   :: ex2,ey2,ez2,hx2,hy2,hz2 ! deallocated in driver
        allocatable   :: zxy, zyx, rhoxy, rhoyx, phaxy, phayx
                  	
	end module femout    
! 
!======================================================================! 
!=============================================================! LHS_MOD
!======================================================================! 
!
        module      lhs_mod
!
! Stores the global finite element lhs matrix in compressed sparse row
! format.  
!     
        complex(8), dimension(:), allocatable :: val
        integer,    dimension(:), allocatable :: col, irw
        integer                               :: nnz 
  
        complex(8), dimension(:),allocatable  ::  vtempex,vtempey,vtempez,vtemp
   
        end module  lhs_mod      
! 
!======================================================================! 
!==========================================================! CSEMARRAYS
!======================================================================! 
!	
	module  csemarrays
	
	integer                    :: nlfreq
	real(8), dimension(:)      :: localfreq
	real(8), dimension(:)      :: errnrm_eh,errnrm, area
!	real(8), dimension(:)      :: egradnrm_eh,egradnrm  !egradnrmeh,
!	complex(8), dimension(:,:) :: esmgradx, esmgrady, esmgradz 

        complex(8), dimension(:)   :: rhs1 !rhs1 is the array for storing the A-phi potentials solutions induced by the x-direcional sources 
        complex(8), dimension(:)   :: rhs2 !rhs2 is the array for storing the A-phi potentials solutions induced by the y-direcional sources

        complex(8), dimension(:)   :: rhs
        integer, parameter         :: nrhs=1        !indicate the number of the right-hand side	
	
	allocatable  :: localfreq  !esmgradx, esmgrady, esmgradz,
	allocatable  :: errnrm_eh,errnrm,area  !, egradnrm_eh,egradnrm  !egradnrmeh,
	allocatable  :: rhs1, rhs2, rhs !  

	end module csemarrays	
!
!======================================================================! 
!==========================================================! MARE3DCSEM
!======================================================================! 
!	
	program MARE3DCSEM
!
!  Driver routine for finite element forward modeling
!  based on unstructured mesh for 3D MT in anisotropic media
!
!  Written by Yixin Ye
!======================================================================! 	
	use fem_params
	use mesh
	use datinfo
	use femout
        use csemarrays
	implicit none
!
! working variables for driver routine
!
	integer        	        :: iadapt  !,i,j,ilay
        character(256)          :: str
	logical        	        :: lrefine,lconv
        real(8)                 :: Time0, Time1
!
! Set the time (s) when the computation starts,
! and any FORTRAN90/95 compiler should support the function cpu_time(). 
!
        call cpu_time(Time0)
!
!  Read in RUNFILE
!
	call readrunfile
!
!  Allocate storage for solution matrices, open log files:
!
        call allocate_solution_arrays  
!             
	   nnodlast = 0 ! stores size of last grid as error check in case
	                ! Triangle errors and doesn't refine grid
	   iadapt = 0
!
! Initialize values at each loop of Tx:
!
           lrefine = .false.
	   lconv   = .false. 
! 
!  Set logical flag for refinement
!
	   if (nfrqrefine.gt.0) lrefine=.true.
!
!  Loop over number of mesh adaptations (refinements)
!
	   do while (iadapt.le.maxadapt+1)
              iadapt = iadapt+1
!
!  If last iteration, don't refine mesh
!
	      if (iadapt.eq.maxadapt+1) lrefine=.false.        
!
!  Update element and node file names to current iTx and iadapt
!
              call updatefilenames(iadapt)
!
!  Read in mesh NODEFILE and ELEMENTFILE
!
	      call readmesh
                       
	      call gen_nodlab  ! nodlabs are for internal and boundary node labels
!
!  Compute forward MT(x,y,z) response 
!
!              write(6,*)'!!'
	      call mt3dfw(lrefine)
!              write(6,*)'!!!'
!
	      call write_mt_response
!
!  Test for convergence, if converged stop, else refine mesh
!
	      if (lrefine) then
		 call testconvergence(lrefine,iadapt,lconv)
		 if (lconv) then   !! compute for remaining frequencies
		    lrefine = .false.
                    if(nlfreq.ne.nfreq)then
	                call mt3dfw(lrefine)
	                call write_mt_response
                    endif                    
		 endif

	      endif
!
! If xyz solution computed, output CSEM fields at receivers to a file 
!		        		
	      if (lrefine) then    		
		 str = trim(tristr)//' '//trim(currentname)                              
		 call system(str)
	      endif			
!
! Move on if not refining
!
	      if (.not.lrefine)  then
		 iadapt = maxadapt+2 ! break out of "do while" loop
	      endif
!    
! Deallocate locale arrays:
!
              deallocate (rhs1,rhs2,rhs)  !
!
!  Deallocate last mesh
!
	      deallocate (xn,yn,zn,nodlab,bcnod,emap,irg)

           enddo  ! loop of iadapt      
!
! Deallocate arrays:
!
	deallocate ( freqvec, frqrefine )	
	deallocate ( ex1,hx1,ey1,hy1,ez1,hz1 )
	deallocate ( ex2,hx2,ey2,hy2,ez2,hz2 )
        deallocate ( zxy, zyx, rhoxy, rhoyx, phaxy, phayx )

! Set the time (s) when the computation ends,
! and any FORTRAN90/95 compiler should support the function cpu_time(). 
!
        call cpu_time(Time1)
!
! output the time this computation costs
!
        print ''
        print '("Total time for computations: ",f10.4," minites.")',(Time1 - Time0)/60.d0
	
	end program MARE3DCSEM
!
!======================================================================!
!================================================! calculate_mt_response
!======================================================================!
!
	subroutine calculate_mt_response(ifreq)	
!
! calculate the MT impetance responses 
! and the apparent resistivities and pahases at the EM sites
!
	use femout
	use datinfo
	use csemarrays	
	implicit none 
		
	integer              :: isite,ifreq
        real(double)  ::omega
        common /kxw/omega
!
! calculate mt solutions at all sites, all frequencies
!
		do isite = 1,nsite
                   zxy(isite,ifreq)=(ex2(isite,ifreq)*hx1(isite,ifreq)-ex1(isite,ifreq)*hx2(isite,ifreq))/(hx1(isite,ifreq)*hy2(isite,ifreq)-hx2(isite,ifreq)*hy1(isite,ifreq))
                   zyx(isite,ifreq)=(ey1(isite,ifreq)*hy2(isite,ifreq)-ey2(isite,ifreq)*hy1(isite,ifreq))/(hx1(isite,ifreq)*hy2(isite,ifreq)-hx2(isite,ifreq)*hy1(isite,ifreq))
                   rhoxy(isite,ifreq)=(abs(zxy(isite,ifreq)))**2/(omega*mu0)
                   rhoyx(isite,ifreq)=(abs(zyx(isite,ifreq)))**2/(omega*mu0)
                   phaxy(isite,ifreq)=atan(dimag(zxy(isite,ifreq))/dreal(zxy(isite,ifreq)))*180.d0/Pi
                   phayx(isite,ifreq)=atan(dimag(zyx(isite,ifreq))/dreal(zyx(isite,ifreq)))*180.d0/Pi
 							   
		enddo
	
	end subroutine calculate_mt_response
!
!======================================================================!
!===================================================! WRITE_mt_RESPONSE
!======================================================================!
!
	subroutine write_mt_response	
!
! Writes out the apparent resistivities and phases responses
! at the EM sites
!
	use femout
	use datinfo
	use csemarrays	
	implicit none 
		
	integer              :: isite,ifreq
        real(double)  ::omega
        common /kxw/omega
!
! Open .mt files for mt responses
!
	open(unit=20,file = trim(currentname)//'.mt',status='replace')
!
! Write out solutions at all sites, all frequencies
!
        write(20,*) nsite, nlfreq

        write(20,fmt = '(3a10, 2x, 5(a16,2x))') 'X','Y','Z', 'frequency(Hz)','resistivity(Rxy)','resistivity(Ryx)','phase(Pxy)','phase(Pyx)'
	          
        do ifreq=1,nlfreq
  
		do isite = 1,nsite

		   write(20,'(3f10.2, 2x, 5(f16.3,2x))') xsite(isite), ysite(isite), zsite(isite),  localfreq(ifreq),  &     
		   rhoxy(isite,ifreq), rhoyx(isite,ifreq),  phaxy(isite,ifreq), phayx(isite,ifreq)
 							   
		enddo
	enddo

	close(20)
	
	deallocate (localfreq) 
	
	end subroutine write_mt_response		
!
!
!======================================================================!
!=============================================! allocate_solution_arrays
!======================================================================!
!
	subroutine allocate_solution_arrays
!
!  Setup storage for the solution arrays 
!
!======================================================================!
        use datinfo
	use femout
!	
! Field matrices in spatial domain	
!
	allocate(ex1(nsite,nfreq),ey1(nsite,nfreq),ez1(nsite,nfreq)) 
	allocate(hx1(nsite,nfreq),hy1(nsite,nfreq),hz1(nsite,nfreq))
	allocate(ex2(nsite,nfreq),ey2(nsite,nfreq),ez2(nsite,nfreq)) 
	allocate(hx2(nsite,nfreq),hy2(nsite,nfreq),hz2(nsite,nfreq))

        allocate(zxy(nsite,nfreq),zyx(nsite,nfreq))
        allocate(rhoxy(nsite,nfreq),rhoyx(nsite,nfreq),phaxy(nsite,nfreq),phayx(nsite,nfreq))
	
        ex1 = 0; ey1 = 0; ez1 = 0; hx1 = 0; hy1 = 0; hz1 = 0
        ex2 = 0; ey2 = 0; ez2 = 0; hx2 = 0; hy2 = 0; hz2 = 0
        zxy = 0; zyx = 0
        rhoxy=0; rhoyx=0; phaxy=0; phayx=0
	
	end subroutine allocate_solution_arrays
!
!======================================================================!
!=====================================================! TESTCONVERGENCE
!======================================================================!	
!	
	subroutine testconvergence(lrefine,iadapt,lconv)
!
! Convergence metric is to take relative differences between present and 
! previous solutions.  Reads in x,y,z solutions from files.
!  
!======================================================================! 	

	use fem_params
	use datinfo
	use mesh
	use csemarrays	
	implicit none
	
	logical      :: lconv,lrefine
        integer      :: iadapt,isite,nlsite,ifreq
	real(8)      :: junk,lfreq
	real(8)      :: e1r,e1i,h1r,h1i,e0r,e0i,h0r,h0i,r1
	integer      :: icthxs,ictexs	
        real(8)      :: maxdexs,maxdhxs,errdhxs,errdexs,rmsexs,rmshxs
        real(8)      :: tol,err
        character(120)  ::stringline
!
! Initialize
!
	lconv   = .false.
!
! See if convergence test is possible and desired
!
	if ( .not.lrefine ) then
	    return
	elseif ((strtmodnum.eq.1).and.(iadapt.eq.1))  then  ! new refinement run
!
! Open a file that will record the convergence all iterations
!	
	    open(unit=11,file=trim(modelroot)//'.conv',status='REPLACE')
	    write(11,100)
	    close(11)
	    write(*,*) modelroot,'debug test'
	    write(*,*) 'press any key to continue!'
	    read(*,*)
	    return
	endif
!
! Initialize some things:
!
	maxdexs = 0d0
	maxdhxs = 0d0
	errdexs = 0d0
	errdhxs = 0d0	
	ictexs = 0
	icthxs = 0
!
! Open .kxsites file for current and previous grids:
!
	open(unit=20,file = trim(currentname)//'.mt',status='OLD')
	read(20,*) nlsite, nlfreq
        read(20,'(a)') stringline
!
	open(unit=21,file = trim(lastname)//'.mt',status='OLD')
	read(21,*) nlsite, nlfreq
        read(21,'(a)') stringline
!
! Loop through all wavenumbers and all sites
!
        do ifreq=1,nlfreq    	    

            do isite = 1,nlsite
		   read(20,*) junk,junk,junk,junk, e1r,e1i,h1r,h1i
		   read(21,*) junk,junk,junk,junk, e0r,e0i,h0r,h0i       
!
! Compute relative diffrences
! 
			if (e1r.gt.ecutoff) then
			    err = dabs (e1r-e0r) / e1r
			    if (err.gt.maxdexs) maxdexs = err
			    errdexs = errdexs+err
			    ictexs = ictexs+1	  	 
			endif

			if (e1i.gt.ecutoff) then
			    err = dabs(e1i-e0i) /e1i
			    if (err.gt.maxdexs) maxdexs = err
			    errdexs = errdexs+err
			    ictexs = ictexs+1	  	 
			endif
             			   		 
	    enddo ! loop over sites

        enddo

	close(20)
	close(21)

        if(ictexs.eq.0)ictexs=1
!        if(icthxs.eq.0)icthxs=1
        rmsexs = errdexs/ictexs
!	rmshxs = dsqrt(errdhxs/icthxs)
!
! Write out another line to the .conv convergence file
!
	open(unit=11,file=trim(modelroot)//'.conv',position='APPEND') 
	write(11,101) iadapt,nnod,rmsexs*100,maxdexs*100
	close(11)
!
! See if requested tolerance obtained:
!
	tol = tolerance/100d0
        if ( (maxdexs.le.tol) ) then 
		lconv = .true.
	endif
!
! Make an announcement to standard output
!	
	write(6,*) ' '
	write(6,*) 'Convergence Test using relative differences: '
	write(6,*) ''
	write(6,102) 'rxy&ryx:','RMS[%]','Max[%]'	    
	write(6,103) rmsexs*100,maxdexs*100
	write(6,*) ''           
        write(6,fmt ='(a,f8.3)') ' Requested tolerance (%): ',tolerance

        if (lconv) then
            write(6,*) 'Solution has converged !'
        else
       	    write(6,*) 'Solution has not converged yet... '
        endif
	return
	
100 format(' Grid,   # Nodes,  mdiff(App.res),   Max(App.res) [%]')	
101 format(i5,1x,i10,2x,f10.1,2x,f10.1)	
102 format(a10,a10,2x,a10)      
103 format(10x,f10.1,2x,f10.1)	
!
	
        end subroutine testconvergence	
!
!======================================================================! 
!===========================================================!  MT3DFW
!======================================================================! 
!
        subroutine mt3dfw(lrefine)
!
!======================================================================!
!                                                                      !
!  Program to solve the coupled partial differential equations for the !  
!  magnetotellurics in arbitrialy anisotropic medias.                  !
!                                                                      !                                                     
!                                                                      !
!======================================================================!
 
	use            fem_params
	use            mesh
	use            datinfo
	use            femout
	use            csemarrays
	use            lhs_mod
	
	implicit none
	
	logical              :: lrefine
	integer              :: ifreq, i, ierr
        real(double)         :: omega,Time0,Time1
        complex(double), allocatable, dimension(:)  :: diag, Sm
!
        common /kxw/omega 
  
!
!        write(*,*)'**'  
        allocate ( diag(4*nnod),Sm(4*nnod) )
        			
!
! Allocate local arrays
!
	call allocatecsem3dfw(lrefine)  
!         write(*,*)'****'
!
!  Compute the FE solutions computions for requested frequencys.
! 
	do ifreq=1,nlfreq
!           call cpu_time(Time0)
!
!  Angular frequency:
            write(*,'(a24,i6)') 'Frequency #:',ifreq    
            omega=2.d0*pi*localfreq(ifreq)		
!
!  Generate boundary condition flag array
!
	    call gen_bcnod ! since this is perturbed in errorestimator, reset here for everywavenumber

            call gen_nodele
!
!  Construct the left hand side (LHS) and right hand side (RHS) of finite element linear system for the MT problem in anisotropic media, 
!  induced by the x-directional sources                   
!
            call gen_diago_lu_ex
!
! Solve the linear system Ax=b using the SSOR_PCG program:
!
            call SSOR_PCG_CSR(4*nnod,NNZ,val,rhs1,col,irw)
!
! Solve the linear system Ax=b using the pardiso program:
!  
!            rhs=rhs1
!            call solver_pardiso
!            rhs1=rhs
	    if (lrefine) then		
		call errorindicator_land(rhs1)   !
		errnrm = errnrm  + errnrm_eh  ! add to previous frequency's error
            endif
!
!  Construct the left hand side (LHS) and right hand side (RHS) of finite element linear system for the MT problem in anisotropic media, 
!  induced by the y-directional sources   
!
            call gen_diago_lu_ey

            call SSOR_PCG_CSR(4*nnod,NNZ,val,rhs2,col,irw) 
!
! Solve the linear system Ax=b using the pardiso program:
! 
!            rhs=rhs2
!            call solver_pardiso    
!            rhs2=rhs           
!
!  If refinement iteration, compute error estimator 
!           
             if (lrefine) then	
		call errorindicator_land(rhs2)    !
		errnrm = errnrm  + errnrm_eh  ! add to previous frequency's error
!
! then output total normed error for each element
!               		
                if (ifreq.eq.nlfreq) then ! sum up Exs and Hxs errors       
		    call write_errnrm(errnrm,area)             
	        endif

	      endif
!
! Compute the EM fields from the A-phi potentials,Lastly, store solutions to  ex1, ex2, ey1, ey2, etc. arrays at the receivers
!
                 call solve_EM_fields_ex(ifreq)
                 call solve_EM_fields_ey(ifreq)   
!
! calculate the apparent resistivity and phase responses from the EM fields
!
                 call calculate_mt_response(ifreq)   
!
!            call cpu_time(Time1)
!            print '("Computational time: ",f10.1," Seconds")',(Time1 - Time0)
!		
        enddo   !move to next frequency 
!
! Remove storage for error estimator variables if refinement iteration
!
   	if (lrefine) then 
	    deallocate ( errnrm_eh,errnrm,area )  !,egradnrmeh,egradnrm_eh,egradnrm
   	endif

        deallocate (vtempex,vtempey,vtempez,vtemp,val,col,irw)  !vtempex,vtempey,vtempez,

        end subroutine mt3dfw

!======================================================================! 
!==============================================================! pardiso
!======================================================================!    
!        subroutine solver_pardiso
!
! Solve the linear system Ax=b using the mkl pardiso:
!
!        use mesh
!        use lhs_mod
!        use csemarrays
!   
!        implicit none
!        integer(8) :: ipt(64)
!        integer :: i, idum, maxfct, mnum, mtype, phase, error, msglvl
!        integer :: iparm(64)
!        integer :: mkl_get_max_threads
!        external   mkl_get_max_threads
!        complex(8) :: ddum
!        complex(8), dimension(:), allocatable:: x
!    
!          allocate (x(4*nnod))
!          maxfct=1
!          mnum=1
!          nrhs=1

!          do i=1,64
!             iparm(i)=0
!             ipt(i)=0
!          enddo
!      iparm(1) = 1 ! no solver default
!      iparm(2) = 2 ! fill-in reordering from METIS
!      iparm(3) = 1 ! numbers of processors
!      iparm(4) = 0 ! no iterative-direct algorithm
!      iparm(5) = 0 ! no user fill-in reducing permutation
!      iparm(6) = 0 ! =0 solution on the first n compoments of x
!      iparm(7) = 0 ! not in use
!      iparm(8) = 10 ! numbers of iterative refinement steps
!      iparm(9) = 0 ! not in use
!      iparm(10) = 8 !   perturbe the pivot elements with 1E-8
!      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
!      iparm(12) = 0 ! not in use
!      iparm(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
!      iparm(14) = 0 ! Output: number of perturbed pivots
!      iparm(15) = 0 ! not in use
!      iparm(16) = 0 ! not in use
!      iparm(17) = 0 ! not in use
!      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
!      iparm(19) = -1 ! Output: Mflops for LU factorization
!      iparm(20) = 0 ! Output: Numbers of CG Iterations
!       iparm(27) = 1 !  checks whether column
!       iparm(21)=  2 !robust
!       iparm(21)=  0 
!preprocessing method based on symmetric weighted
!matchings and 1x1 and 2x2 Bunch and Kaufman pivoting
!will be used in the factorization process.
!indices are sorted in increasing order within each row.
!      iparm(60) = 2 ! out of core (when compuation memory is insufficent)
    ! iparm(60) = 1
! 
!      error = 0 ! initialize error flag
!      msglvl = 1 ! print statistical information
!      mtype = 6 ! symmetric, complex             
!     
! First, factorize the matrix. The factors are stored in factor() handle. 
!
!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization
!      phase = 11 ! only reordering and symbolic factorization

!      call pardiso( ipt, maxfct, mnum, mtype, phase, 4*nnod, val, irw, col, idum, &
!			        & nrhs, iparm, msglvl, ddum, ddum, error) 
!      WRITE(*,*) 'Reordering completed ... '
!      IF (error .NE. 0) THEN
!         WRITE(*,*) 'The following ERROR was detected: ', error
!         STOP
!      END IF
!      WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!      WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
!      print *,'iparm(15) is ',iparm(15)
!                                  
!.. Factorization.
!      phase = 22 ! only factorization
!      call pardiso( ipt, maxfct, mnum, mtype, phase, 4*nnod, val, irw, col, idum, &
!			        & nrhs, iparm, msglvl, ddum, ddum, error) 
!      WRITE(*,*) 'Factorization completed ... '
!      IF (error .NE. 0) THEN
!         WRITE(*,*) 'The following ERROR was detected: ', error
!         STOP
!      ENDIF
! Second, solve the linear system using the pardiso. 
!
!.. Back substitution and iterative refinement
!           iparm(8) = 2 ! max numbers of iterative refinement steps
!           phase = 33 ! solve the linear system
!	   call pardiso( ipt, maxfct, mnum, mtype, phase, 4*nnod, val, irw, col, idum, &
!			        & nrhs, iparm, msglvl, rhs, x, error) 
!
!           write(*,*)'Solve completed...'
!           if(error.ne.0)then
!              write(*,*)'The following ERROR was detected:',error
!              stop
!           endif

!           rhs=x
!                                
!           phase=-1  ! release internal memory!
!	   call pardiso( ipt, maxfct, mnum, mtype, phase, 4*nnod, ddum, idum,idum, idum, &
!			        & nrhs, iparm, msglvl, ddum, ddum, error) 
!   do i = 1 , 2*nnod
!	write(*,'(a,i,e16.8,e16.8)') 'i=',i,rhs(i)
!   enddo
!           deallocate(X)
!
!        end subroutine solver_pardiso
!       
!======================================================================! 
!======================================================! solve_EM_fields
!======================================================================!
!    
        subroutine solve_EM_fields_ex(ifreq)
!
! solve the EM fields from the A-phi potentials induced by the x-directional sources 
! using the moving least square method.
!
        use mesh
        use lhs_mod
        use csemarrays
        use datinfo
        use femout
!  
        implicit none
!
        INTEGER,PARAMETER ::LWORK=200
        integer :: i, iele,isite,ifreq, nnzV, innz, ineednode, inod1, inod(4),rank,inf,nrhs0
        integer, dimension(200) :: cols
        real(8) :: hd, w, x1, y1, z1,xp,yp,zp
        real(8), dimension(200):: d, xl, yl, zl
        real(8), dimension(:,:), allocatable:: A,B
        real(8) :: work(LWORK)
        real(8)    :: s(4),rwork(20)
        real(double)  ::omega
        common /kxw/omega


         nrhs0=8

        do isite=1,nsite
           xp = xsite(isite)
           yp = ysite(isite)
           zp = zsite(isite)
        do inod1=1,nnod
           
           if (xn(inod1).eq.xp.and.yn(inod1).eq.yp.and.zn(inod1).eq.zp) then
              nnzV = 0
              cols = 0
              x1=xn(inod1)
              y1=yn(inod1)
              z1=zn(inod1)

              do iele =1, nelein(inod1)          
                 e = nodele(inod1,iele)     ! get element number

                 if(irg(e).ne.1)then
                 inod(1) = emap(1,e)
                 inod(2) = emap(2,e)
                 inod(3) = emap(3,e)
                 inod(4) = emap(4,e)

                 do i=1,4
                    ineednode = 1 

                    do innz=1,nnzV
                       if (cols(innz).eq.inod(i)) ineednode = 0
                    enddo  
     
                    if (ineednode.eq.1) then
                       nnzV = nnzV+1
                       cols(nnzV) = inod(i)
                    endif
                 enddo

                 endif

              enddo
           
              allocate (A(nnzv,4),B(nnzv,nrhs0))

              do i=1,nnzv
                 xl(i)=xn(cols(i))
                 yl(i)=yn(cols(i))
                 zl(i)=zn(cols(i))
                 d(i) =(xl(i)-x1)**2+(yl(i)-y1)**2+(zl(i)-z1)**2 
              enddo
          
              hd=maxval(d(1:nnzv))
              do i=1,nnzv
                 w=exp(-d(i)/hd)
                 A(i,1)=w*xl(i)
                 A(i,2)=w*yl(i)
                 A(i,3)=w*zl(i)
                 A(i,4)=w
                 B(i,1)  =w*dreal(rhs1(4*cols(i)))
                 B(i,2)  =w*dimag(rhs1(4*cols(i)))              
                 B(i,3)  =w*dreal(rhs1(4*cols(i)-3))
                 B(i,4)  =w*dimag(rhs1(4*cols(i)-3))
                 B(i,5)  =w*dreal(rhs1(4*cols(i)-2))
                 B(i,6)  =w*dimag(rhs1(4*cols(i)-2))
                 B(i,7)  =w*dreal(rhs1(4*cols(i)-1))
                 B(i,8)  =w*dimag(rhs1(4*cols(i)-1))

              enddo

              call dgelss(nnzv,4,nrhs0,A,nnzv,B,nnzv,s,0.0d0,rank,work,LWORK,inf)
               
              ex1(isite,ifreq)=ic*omega*(rhs1(4*inod1-3)+B(1,1)+ic*B(1,2))
              ey1(isite,ifreq)=ic*omega*(rhs1(4*inod1-2)+B(2,1)+ic*B(2,2))
              ez1(isite,ifreq)=ic*omega*(rhs1(4*inod1-1)+B(3,1)+ic*B(3,2))

              hx1(isite,ifreq)= (B(2,7)+ic*B(2,8)-B(3,5)-ic*B(3,6))/mu0
              hy1(isite,ifreq)= (B(3,3)+ic*B(3,4)-B(1,7)-ic*B(1,8))/mu0
              hz1(isite,ifreq)= (B(1,5)+ic*B(1,6)-B(2,3)-ic*B(2,4))/mu0

              deallocate(A,B)

              exit  

           endif

        enddo 

        enddo
                                     
        end subroutine solve_EM_fields_ex
!======================================================================! 
!======================================================! solve_EM_fields
!======================================================================!
!    
        subroutine solve_EM_fields_ey(ifreq)
!
! solve the EM fields from the A-phi potentials induced by the y-directional sources 
! using the moving least square method.
!
        use mesh
        use lhs_mod
        use csemarrays
        use datinfo
        use femout    
        implicit none
!
        INTEGER,PARAMETER ::LWORK=200
        integer :: i, iele, nnzV,isite,ifreq, innz, ineednode, inod1, inod(4),rank,inf,nrhs0
        integer, dimension(200) :: cols
        real(8) :: hd, w, x1, y1, z1,xp,yp,zp
        real(8), dimension(200):: d, xl, yl, zl
        real(8), dimension(:,:), allocatable:: A,B
        real(8) :: work(LWORK)
        real(8)    :: s(4),rwork(20)
        real(double)  ::omega
        common /kxw/omega

! 

         nrhs0=8
        do isite=1,nsite
           xp = xsite(isite)
           yp = ysite(isite)
           zp = zsite(isite)
        
        do inod1=1,nnod
           
           if (xn(inod1).eq.xp.and.yn(inod1).eq.yp.and.zn(inod1).eq.zp) then

              nnzV = 0
              cols = 0
              x1=xn(inod1)
              y1=yn(inod1)
              z1=zn(inod1)

              do iele =1, nelein(inod1)          
                 e = nodele(inod1,iele)     ! get element number

                 if(irg(e).ne.1)then
                 inod(1) = emap(1,e)
                 inod(2) = emap(2,e)
                 inod(3) = emap(3,e)
                 inod(4) = emap(4,e)

                 do i=1,4
                    ineednode = 1 

                    do innz=1,nnzV
                       if (cols(innz).eq.inod(i)) ineednode = 0
                    enddo  
     
                    if (ineednode.eq.1) then
                       nnzV = nnzV+1
                       cols(nnzV) = inod(i)
                    endif
                 enddo

                 endif

              enddo
           
              allocate (A(nnzv,4),B(nnzv,nrhs0))

              do i=1,nnzv
                 xl(i)=xn(cols(i))
                 yl(i)=yn(cols(i))
                 zl(i)=zn(cols(i))
                 d(i) =(xl(i)-x1)**2+(yl(i)-y1)**2+(zl(i)-z1)**2 
              enddo
          
              hd=maxval(d(1:nnzv))
              do i=1,nnzv
                 w=exp(-d(i)/hd)
                 A(i,1)=w*xl(i)
                 A(i,2)=w*yl(i)
                 A(i,3)=w*zl(i)
                 A(i,4)=w
                 B(i,1)  =w*dreal(rhs2(4*cols(i)))
                 B(i,2)  =w*dimag(rhs2(4*cols(i)))              
                 B(i,3)  =w*dreal(rhs2(4*cols(i)-3))
                 B(i,4)  =w*dimag(rhs2(4*cols(i)-3))
                 B(i,5)  =w*dreal(rhs2(4*cols(i)-2))
                 B(i,6)  =w*dimag(rhs2(4*cols(i)-2))
                 B(i,7)  =w*dreal(rhs2(4*cols(i)-1))
                 B(i,8)  =w*dimag(rhs2(4*cols(i)-1))

              enddo

              call dgelss(nnzv,4,nrhs0,A,nnzv,B,nnzv,s,0.0d0,rank,work,LWORK,inf)
               
              ex2(isite,ifreq)=ic*omega*(rhs2(4*inod1-3)+B(1,1)+ic*B(1,2))
              ey2(isite,ifreq)=ic*omega*(rhs2(4*inod1-2)+B(2,1)+ic*B(2,2))
              ez2(isite,ifreq)=ic*omega*(rhs2(4*inod1-1)+B(3,1)+ic*B(3,2))

              hx2(isite,ifreq)= (B(2,7)+ic*B(2,8)-B(3,5)-ic*B(3,6))/mu0
              hy2(isite,ifreq)= (B(3,3)+ic*B(3,4)-B(1,7)-ic*B(1,8))/mu0
              hz2(isite,ifreq)= (B(1,5)+ic*B(1,6)-B(2,3)-ic*B(2,4))/mu0

              deallocate(A,B)  

              exit

           endif

        enddo 

        enddo

        deallocate( nodele, nelein,indinod1)
                                     
        end subroutine solve_EM_fields_ey
!    
!======================================================================! 
!===================================================!  ALLOCATECSEM3DFW	
!======================================================================! 
!	
	subroutine allocatecsem3dfw(lrefine)	
	
	use csemarrays
	use datinfo
	use mesh
	use lhs_mod	
	implicit none	
	logical        :: lrefine 
!
        if (lrefine) then 
!
	    allocate ( errnrm_eh(nele),errnrm(nele),area(nele)) 
!            allocate ( egradnrm_eh(nele),egradnrm(nele) )  !,egradnrmeh(4*nele)
	    errnrm = 0	 
!            egradnrm = 0
!           egradnrmeh = 0
!
! Set local frequencies to use for refinement
!
            nlfreq = nfrqrefine
            allocate (localfreq(nlfreq))
            localfreq = frqrefine
	
	else  ! full solution requested
!
! Set local frequencies to use for final computation
!
            nlfreq = nfreq
            allocate (localfreq(nlfreq))
            localfreq = freqvec 				
   	endif				
! 	   	
!  Set aside storage.
!       
        allocate ( rhs1(4*nnod),rhs(4*nnod)) !
        allocate ( rhs2(4*nnod)) !
	allocate ( val(50*4*nnod), col(50*4*nnod), irw(4*nnod+1))
        allocate ( vtempEX(4*nnod), vtempEY(4*nnod), vtempEZ(4*nnod), vtemp(4*nnod))  

	end subroutine allocatecsem3dfw
!
!======================================================================! 
!===========================================================!  SUMERROR	
!======================================================================! 
!	
!	subroutine sumerror
	
!	use datinfo
!	use mesh
!	use csemarrays	
!	implicit none
           	 
!	if (errortype.eq.1) then   !use M potential
!           do e=1,nele    	
!              errnrm(e) = MAX(errnrmeh(4*e-3),errnrmeh(4*e-2),errnrmeh(4*e-1))!
!	   enddo!
!	elseif (errortype.eq.2) then  !use E potential
!	   do e=1,nele    	
!	      errnrm(e) = errnrmeh(4*e)
!	   enddo
!	else   !Use both E and M potentials
!	   do e=1,nele    	
!	      errnrm(e) = MAX(errnrmeh(4*e-3),errnrmeh(4*e-2),errnrmeh(4*e-1),errnrmeh(4*e)) 
!	   enddo
!	endif

!	call write_errnrm(errnrm,area)
	
!	end subroutine sumerror 
!	      
!   	 
!======================================================================! 
!==========================================================! GEN_NODLAB
!======================================================================! 

        subroutine gen_nodlab

!======================================================================!
!  Program to figure out where the mesh boundaries are.
!
! Strategy: 
! * To use with any triangle emap and node array, regardless of whether
! boundary markers are defined for the nodes.  This algorithm could easily
! be extended to add in assign nodlab with the seafloor boundary
!
! 0. set nodlab = 0
! 1. define y min and max, z min and max
! 2. loop over elements, then over nodes
! 3. cases: 
!   A. if on top boundary (usually top of air), set nodlab = 1 
!   B. elseif on bottom boundary (z=zmax), set nodlab = 3 
!   C. elseif node is on left or right boundary (y = ymin or ymax)
!          set nodlab=2                    
!   D. else node is on interior of mesh (y and z don't equal y or z min or max)
!          set nodlab = 4 
!
!======================================================================!
        use mesh
        use fem_params
        IMPLICIT NONE

        real(double)  :: xmin,xmax,ymin,ymax,zmin,zmax 
        integer       :: i
        
! 0. set nodlab = 0
        do i=1,nnod
            nodlab(i) = 0
        end do 
! 1. define y min and max, z min and max
        xmin =   1.d99
        xmax =  -1.d99
        ymin =   1.d99
        ymax =  -1.d99
        zmin =   1.d99
        zmax =  -1.d99

        do i=1,nnod
            if (xn(i).lt.xmin) xmin = xn(i)
            if (xn(i).gt.xmax) xmax = xn(i)
            if (yn(i).lt.ymin) ymin = yn(i)
            if (yn(i).gt.ymax) ymax = yn(i)
            if (zn(i).lt.zmin) zmin = zn(i)
            if (zn(i).gt.zmax) zmax = zn(i)
        enddo

! 2. loop over elements, then over 3 nodes

!        do e = 1,nele
        do i = 1,nnod
! 3. cases: 
!   A. if on top boundary (z=zmin), set nodlab = 1
            if ( zn(i).eq.zmin ) then
                   nodlab(i) = 1 
!   B. elseif on bottom boundary (z=zmax), set nodlab = 3 
            elseif ( zn(i).eq.zmax ) then
                   nodlab(i) = 3
!   C. elseif node is on left or right boundary (y = ymin or ymax)
!          set nodlab=2              
            elseif ((yn(i).eq.ymax).or.(yn(i).eq.ymin)) then             
                   nodlab(i) = 2
!    C. elseif node is on front or back boundary (x = xmin or xmax)
            elseif ((xn(i).eq.xmax).or.(xn(i).eq.xmin)) then
                   nodlab(i) = 5
!   D. else node is on interior of mesh (x, y and z don't equal x, y or z min or max)
!          set nodlab = 4           
            else  
                   nodlab(i) = 4 
            endif
        enddo
!        end do   
!  End of gen_nodlab subrountine
!
        end subroutine gen_nodlab      
!
!======================================================================! 
!===========================================================! GEN_BCNOD
!======================================================================! 

        subroutine gen_bcnod
!======================================================================!
!                                                                      !
!  Routine to remap the vector of node labels 'nodlab' to a new vector !
!  of labels 'bcnod' whose values are appropriate.This routine expects !
!  the following scheme for nodlab:                                    !
!                                                                      !
!     nodlab = 1  ...  node on top of model (top of air usually)       !
!     nodlab = 2  ...  node on left or right boundary                  !
!     nodlab = 3  ...  node on bottom bounrary                         !
!     nodlab = 4  ...  node in EARTH interior region                   ! 
!                                                                      !
!  The routine returns the following in bcnod:                         !
!                                                                      !
!     nodlab          bcnod                                            !
!          1            1                                              !
!          2            1                                              !
!          3            1                                              !
!          4            0                                              !
!                                                                      !
!  Values within bcnod indicate the fate of the node they represent:   !
!                                                                      !
!    0  ...  interior node where FE solution is obtained.              !
!    1  ...  boundary node where Dirichlet condition is applied.       !
!                                                                      !
!======================================================================!

        use           fem_params
        use           mesh
        IMPLICIT NONE
        integer                  :: i
      
        bcnod = 0
!
!  Remap the node labels 
!
        do i=1,nnod
           if (nodlab(i).ne.4) bcnod(i) = 1
        enddo

        end subroutine gen_bcnod
!
!======================================================================! 
!===============================================================! gen_k1
!======================================================================! 
!
        subroutine gen_k1(xe,ye,ze,del,ommu,sige,k1)
!======================================================================!
!                                                                      ! 
!  Routine to compute the 16x16 element matrix K1 on a tetrahedral     !
!  element.                                                            !
!                                                                      ! 
!  written by yixin ye, 2016.1.6                                       !
!                                                                      ! 
!======================================================================! 
        use mesh
        IMPLICIT NONE
        complex(8) :: k0,ommu
        complex(8) :: k1(16,16)
        real(8), dimension(4) :: xe,ye,ze
        real(8), dimension(4) :: a,b,c
        real(8) :: M1(4,4),M2(4,4),M3(4,4),M4(4,4)
        real(8) :: M5(4,4),M6(4,4),M7(4,4),M8(4,4)
        real(8) :: del,ck1,ck2,ck3,a1,a2,b1,b2,area,sige(6)

        k1=(0.d0,0.d0)
!
!  Computer the coefficients of linear interpolation in an element.
!
        a(1) = ye(4)*(ze(3)-ze(2)) + ye(2)*(ze(4)-ze(3)) + ye(3)*(ze(2)-ze(4))
        a(2) = ye(4)*(ze(1)-ze(3)) + ye(1)*(ze(3)-ze(4)) + ye(3)*(ze(4)-ze(1))
        a(3) = ye(4)*(ze(2)-ze(1)) + ye(1)*(ze(4)-ze(2)) + ye(2)*(ze(1)-ze(4))
        a(4) = ye(3)*(ze(1)-ze(2)) + ye(1)*(ze(2)-ze(3)) + ye(2)*(ze(3)-ze(1))

        b(1) = xe(3)*(ze(4)-ze(2)) + xe(4)*(ze(2)-ze(3)) + xe(2)*(ze(3)-ze(4))
        b(2) = xe(3)*(ze(1)-ze(4)) + xe(4)*(ze(3)-ze(1)) + xe(1)*(ze(4)-ze(3))
        b(3) = xe(2)*(ze(4)-ze(1)) + xe(4)*(ze(1)-ze(2)) + xe(1)*(ze(2)-ze(4))
        b(4) = xe(2)*(ze(1)-ze(3)) + xe(3)*(ze(2)-ze(1)) + xe(1)*(ze(3)-ze(2))

        c(1) = xe(4)*(ye(3)-ye(2)) + xe(2)*(ye(4)-ye(3)) + xe(3)*(ye(2)-ye(4))
        c(2) = xe(4)*(ye(1)-ye(3)) + xe(1)*(ye(3)-ye(4)) + xe(3)*(ye(4)-ye(1))
        c(3) = xe(4)*(ye(2)-ye(1)) + xe(1)*(ye(4)-ye(2)) + xe(2)*(ye(1)-ye(4))
        c(4) = xe(3)*(ye(1)-ye(2)) + xe(1)*(ye(2)-ye(3)) + xe(2)*(ye(3)-ye(1))
!
!  Compute the volume of the current element
!
        del = dabs((xe(2)-xe(1))*(ye(3)*ze(4)-ze(3)*ye(4)) + (xe(1)- xe(3))*(ye(2)*ze(4)-ze(2)*ye(4)) + (xe(4)-xe(1))*(ye(2)*ze(3)-ze(2)*ye(3)) &
             &  +(xe(3)- xe(2))*(ye(1)*ze(4)-ze(1)*ye(4)) + (xe(2)- xe(4))*(ye(1)*ze(3)-ze(1)*ye(3)) + (xe(4)- xe(3))*(ye(1)*ze(2)-ze(1)*ye(2)))/6.d0

        M1(:,1) = (/a(1)*a(1), a(1)*a(2), a(1)*a(3), a(1)*a(4)/)
        M1(:,2) = (/a(2)*a(1), a(2)*a(2), a(2)*a(3), a(2)*a(4)/)
        M1(:,3) = (/a(3)*a(1), a(3)*a(2), a(3)*a(3), a(3)*a(4)/)
        M1(:,4) = (/a(4)*a(1), a(4)*a(2), a(4)*a(3), a(4)*a(4)/) !check

        M2(:,1) = (/b(1)*b(1), b(1)*b(2), b(1)*b(3), b(1)*b(4)/)
        M2(:,2) = (/b(2)*b(1), b(2)*b(2), b(2)*b(3), b(2)*b(4)/)
        M2(:,3) = (/b(3)*b(1), b(3)*b(2), b(3)*b(3), b(3)*b(4)/)
        M2(:,4) = (/b(4)*b(1), b(4)*b(2), b(4)*b(3), b(4)*b(4)/) !check

        M3(:,1) = (/c(1)*c(1), c(1)*c(2), c(1)*c(3), c(1)*c(4)/)
        M3(:,2) = (/c(2)*c(1), c(2)*c(2), c(2)*c(3), c(2)*c(4)/)
        M3(:,3) = (/c(3)*c(1), c(3)*c(2), c(3)*c(3), c(3)*c(4)/)
        M3(:,4) = (/c(4)*c(1), c(4)*c(2), c(4)*c(3), c(4)*c(4)/) !check

        M4(:,1) = (/2, 1, 1, 1/)
        M4(:,2) = (/1, 2, 1, 1/)
        M4(:,3) = (/1, 1, 2, 1/)
        M4(:,4) = (/1, 1, 1, 2/) !check

        M5(:,1) = (/a(1), a(1), a(1), a(1)/)
        M5(:,2) = (/a(2), a(2), a(2), a(2)/)
        M5(:,3) = (/a(3), a(3), a(3), a(3)/)
        M5(:,4) = (/a(4), a(4), a(4), a(4)/) !check

        M6(:,1) = (/b(1), b(1), b(1), b(1)/)
        M6(:,2) = (/b(2), b(2), b(2), b(2)/)
        M6(:,3) = (/b(3), b(3), b(3), b(3)/)
        M6(:,4) = (/b(4), b(4), b(4), b(4)/) !check 

        M7(:,1) = (/c(1), c(1), c(1), c(1)/)
        M7(:,2) = (/c(2), c(2), c(2), c(2)/)
        M7(:,3) = (/c(3), c(3), c(3), c(3)/)
        M7(:,4) = (/c(4), c(4), c(4), c(4)/) !check

        M8(1,1) = sige(1)*a(1)*a(1)+sige(2)*a(1)*b(1)+sige(3)*a(1)*c(1)+sige(2)*b(1)*a(1)+sige(4)*b(1)*b(1)+sige(5)*b(1)*c(1)+sige(3)*c(1)*a(1)+sige(5)*c(1)*b(1)+sige(6)*c(1)*c(1) 
        M8(1,2) = sige(1)*a(1)*a(2)+sige(2)*a(1)*b(2)+sige(3)*a(1)*c(2)+sige(2)*b(1)*a(2)+sige(4)*b(1)*b(2)+sige(5)*b(1)*c(2)+sige(3)*c(1)*a(2)+sige(5)*c(1)*b(2)+sige(6)*c(1)*c(2)
        M8(1,3) = sige(1)*a(1)*a(3)+sige(2)*a(1)*b(3)+sige(3)*a(1)*c(3)+sige(2)*b(1)*a(3)+sige(4)*b(1)*b(3)+sige(5)*b(1)*c(3)+sige(3)*c(1)*a(3)+sige(5)*c(1)*b(3)+sige(6)*c(1)*c(3)
        M8(1,4) = sige(1)*a(1)*a(4)+sige(2)*a(1)*b(4)+sige(3)*a(1)*c(4)+sige(2)*b(1)*a(4)+sige(4)*b(1)*b(4)+sige(5)*b(1)*c(4)+sige(3)*c(1)*a(4)+sige(5)*c(1)*b(4)+sige(6)*c(1)*c(4)
        M8(2,1) = M8(1,2)
        M8(2,2) = sige(1)*a(2)*a(2)+sige(2)*a(2)*b(2)+sige(3)*a(2)*c(2)+sige(2)*b(2)*a(2)+sige(4)*b(2)*b(2)+sige(5)*b(2)*c(2)+sige(3)*c(2)*a(2)+sige(5)*c(2)*b(2)+sige(6)*c(2)*c(2)
        M8(2,3) = sige(1)*a(2)*a(3)+sige(2)*a(2)*b(3)+sige(3)*a(2)*c(3)+sige(2)*b(2)*a(3)+sige(4)*b(2)*b(3)+sige(5)*b(2)*c(3)+sige(3)*c(2)*a(3)+sige(5)*c(2)*b(3)+sige(6)*c(2)*c(3)
        M8(2,4) = sige(1)*a(2)*a(4)+sige(2)*a(2)*b(4)+sige(3)*a(2)*c(4)+sige(2)*b(2)*a(4)+sige(4)*b(2)*b(4)+sige(5)*b(2)*c(4)+sige(3)*c(2)*a(4)+sige(5)*c(2)*b(4)+sige(6)*c(2)*c(4)
        M8(3,1) = M8(1,3)
        M8(3,2) = M8(2,3)
        M8(3,3) = sige(1)*a(3)*a(3)+sige(2)*a(3)*b(3)+sige(3)*a(3)*c(3)+sige(2)*b(3)*a(3)+sige(4)*b(3)*b(3)+sige(5)*b(3)*c(3)+sige(3)*c(3)*a(3)+sige(5)*c(3)*b(3)+sige(6)*c(3)*c(3)
        M8(3,4) = sige(1)*a(3)*a(4)+sige(2)*a(3)*b(4)+sige(3)*a(3)*c(4)+sige(2)*b(3)*a(4)+sige(4)*b(3)*b(4)+sige(5)*b(3)*c(4)+sige(3)*c(3)*a(4)+sige(5)*c(3)*b(4)+sige(6)*c(3)*c(4)
        M8(4,1) = M8(1,4)
        M8(4,2) = M8(2,4)
        M8(4,3) = M8(3,4)
        M8(4,4) = sige(1)*a(4)*a(4)+sige(2)*a(4)*b(4)+sige(3)*a(4)*c(4)+sige(2)*b(4)*a(4)+sige(4)*b(4)*b(4)+sige(5)*b(4)*c(4)+sige(3)*c(4)*a(4)+sige(5)*c(4)*b(4)+sige(6)*c(4)*c(4)
         
        ck1 = 1.d0/36.d0/del
        ck2 = del/20.d0
        ck3 = 1.d0/24.d0

        k1(1:4,1:4)    = -ck1*(M1+M2+M3)+ommu*sige(1)*ck2*M4 
        k1(5:8,5:8)    = -ck1*(M1+M2+M3)+ommu*sige(4)*ck2*M4
        k1(9:12,9:12)  = -ck1*(M1+M2+M3)+ommu*sige(6)*ck2*M4

        k1(1:4,5:8)    = ommu*sige(2)*ck2*M4
        k1(1:4,9:12)   = ommu*sige(3)*ck2*M4
        k1(5:8,9:12)   = ommu*sige(5)*ck2*M4

        k1(5:8,1:4)    = k1(1:4,5:8)
        k1(9:12,1:4)   = k1(1:4,9:12)
        k1(9:12,5:8)   = k1(5:8,9:12)          

        k1(1:4,13:16)  = ommu*ck3*(sige(1)*M5+sige(2)*M6+sige(3)*M7)
        k1(5:8,13:16)  = ommu*ck3*(sige(2)*M5+sige(4)*M6+sige(5)*M7)
        k1(9:12,13:16) = ommu*ck3*(sige(3)*M5+sige(5)*M6+sige(6)*M7)

        k1(13:16,1:4)   = ommu*ck3*transpose(sige(1)*M5+sige(2)*M6+sige(3)*M7)  
        k1(13:16,5:8)   = ommu*ck3*transpose(sige(2)*M5+sige(4)*M6+sige(5)*M7) 
        k1(13:16,9:12)  = ommu*ck3*transpose(sige(3)*M5+sige(5)*M6+sige(6)*M7)

        k1(13:16,13:16) = ommu*ck1*M8        
        return

        end subroutine gen_k1   

!======================================================================!
!===========================================================! gen_nodele
!======================================================================!

        subroutine gen_nodele
        use mesh

        implicit none
        integer       :: inod(4)

        allocate( nodele(nnod,120), nelein(nnod), indinod1(nnod,120))
!
! Store array list of elements containing each node
!
        nelein = 0
        nodele = 0
        indinod1 = 0

        do e = 1,nele   
           inod(1) = emap(1,e)
           inod(2) = emap(2,e)
           inod(3) = emap(3,e)
           inod(4) = emap(4,e)

           nelein(inod(1)) = nelein(inod(1)) + 1
           nelein(inod(2)) = nelein(inod(2)) + 1
           nelein(inod(3)) = nelein(inod(3)) + 1
           nelein(inod(4)) = nelein(inod(4)) + 1
        
           nodele(inod(1),nelein(inod(1)))   = e ! record this element
           nodele(inod(2),nelein(inod(2)))   = e ! record this element
           nodele(inod(3),nelein(inod(3)))   = e ! record this element
           nodele(inod(4),nelein(inod(4)))   = e ! record this element
         
           indinod1(inod(1),nelein(inod(1))) = 1 ! store which node of emap(1:3,e) is inod1
           indinod1(inod(2),nelein(inod(2))) = 2 ! store which node of emap(1:3,e) is inod1
           indinod1(inod(3),nelein(inod(3))) = 3 ! store which node of emap(1:3,e) is inod1  
           indinod1(inod(4),nelein(inod(4))) = 4 ! store which node of emap(1:3,e) is inod1          
        enddo

!        open(10,file='nelein.dat')

!        do i=1,nnod
!           write(10,*)nelein(i)
!        enddo

!        close(10)     
        end subroutine  gen_nodele 
!
!
!======================================================================!
!=========================================================! gen_diago_lu
!======================================================================!
!
        subroutine gen_diago_lu_ex
!
! Compressed Row storage format:
!
!  val - nnz x 1 array of non-zero values of matrix, stored in order of 
!        the rows from 1 to nrows
!  col - nnz x 1 array of column number of each non-zero entry 
!        (doesn't have to be in order...)
!  irw - nnod+1 x 1 array, ith entry is points to first entry of the ith row in val.
!        last entry is equal to nnz.
!
! Modules used:
!
        use fem_params
        use lhs_mod
        use mesh
        use datinfo
        use csemarrays

        implicit none
    
        integer        :: i,i1,jj,n,inod(4),inod1,eps(8)
        real(8)        :: xe(4),ye(4),ze(4),del    
        integer, dimension(400) :: cols
        integer ::iele, nnzV, innz, ineednode
        complex(8)    :: ommu,k0,ke1(16,16)
        real(double)  ::omega,sige(6)
        common /kxw/omega
    
        ommu=ic*omega*mu0
        val  = (0.d0,0.d0)
        col = 0
        irw = 0
        nnz = 0   
        rhs1 = (0.d0,0.d0) 
!        rhs2 = (0.d0,0.d0)
!        rhs  = (0.d0,0.d0)
! 
! Loop over nodes, 
!
!        open(11,file='nnz.dat')

        do inod1 = 1,nnod ! note there are 4*nnod unknowns (Ax Ay Az and phi)         
           vtempEX = (0.d0,0.d0)
           vtempEY = (0.d0,0.d0)
           vtempEZ = (0.d0,0.d0)
           vtemp   = (0.d0,0.d0)

           nnzV = 0
           cols = 0
!                     
! Loop over elements containing inod1
!  
         if (nodlab(inod1).ne.1.and.nodlab(inod1).ne.3) then   ! inod1 that doesnot lied on the top or bottom boundary

            do iele =1, nelein(inod1)           
              e = nodele(inod1,iele)      ! get element number         
              i = indinod1(inod1,iele)
              inod(1) = emap(1,e)
              inod(2) = emap(2,e)
              inod(3) = emap(3,e)
              inod(4) = emap(4,e)
              jj=irg(e)
!
!  Yank out the (y,z) coordinates for the current element 
!        
              do n=1,4
                 xe(n)=xn(inod(n))
                 ye(n)=yn(inod(n))
                 ze(n)=zn(inod(n))
              enddo

              sige=sg(jj,:)

              call gen_K1(xe,ye,ze,del,ommu,sige,Ke1) 
!	                     
!  Wipe out the row and column in the mass and stiffness 
!  matrices for those elements which lie in the air
!  of the FE mesh since a phi value is zero there.
!
              do n=1,4
                 if(nodlab(inod(n)).eq.2.or.nodlab(inod(n)).eq.1.or.nodlab(inod(n)).eq.5) then  ! node lied on the lef-right boundary
                     ke1(:,n+4)=(0.d0,0.d0)
                     ke1(:,n+8)=(0.d0,0.d0)
                     ke1(:,n+12)=(0.d0,0.d0)
                 elseif (nodlab(inod(n)).eq.3) then     ! node lied on the bottom boundary              
                     ke1(:,n)=(0.d0,0.d0)
                     ke1(:,n+4)=(0.d0,0.d0)
                     ke1(:,n+8)=(0.d0,0.d0)
                     ke1(:,n+12)=(0.d0,0.d0)
                 endif
              enddo                   
!
                 if(jj.eq.1)then !node doesnot lied on the top boundary, but in the air 
                    do n=1,4
                       ke1(:,n+12)=(0.d0,0.d0)
                    enddo
                 endif  
!                                                    
!! 
!! Append parts corresponding to inod1's E and H vectors to a
!! temporary vector of nnod*2 length:        
!!                                    
              do n=1,4
                 vtempex(4*inod(n)-3) = vtempex(4*inod(n)-3) + ke1(i,n)           
                 vtempex(4*inod(n)-2) = vtempex(4*inod(n)-2) + ke1(i,n+4)      
                 vtempex(4*inod(n)-1) = vtempex(4*inod(n)-1) + ke1(i,n+8)    
                 vtempex(4*inod(n)  ) = vtempex(4*inod(n)  ) + ke1(i,n+12) 

                                                                      
                 vtempey(4*inod(n)-3) = vtempey(4*inod(n)-3) + ke1(i+4,n)           
                 vtempey(4*inod(n)-2) = vtempey(4*inod(n)-2) + ke1(i+4,n+4)      
                 vtempey(4*inod(n)-1) = vtempey(4*inod(n)-1) + ke1(i+4,n+8)    
                 vtempey(4*inod(n)  ) = vtempey(4*inod(n)  ) + ke1(i+4,n+12)
                                                                     
                 vtempez(4*inod(n)-3) = vtempez(4*inod(n)-3) + ke1(i+8,n)           
                 vtempez(4*inod(n)-2) = vtempez(4*inod(n)-2) + ke1(i+8,n+4)      
                 vtempez(4*inod(n)-1) = vtempez(4*inod(n)-1) + ke1(i+8,n+8)   
                 vtempez(4*inod(n)  ) = vtempez(4*inod(n)  ) + ke1(i+8,n+12)

                 if(jj.ne.1)then ! as the element lied in the air, the phi value is set to zero
                    vtemp(4*inod(n)-3)   = vtemp(4*inod(n)-3)   + ke1(i+12,n)           
                    vtemp(4*inod(n)-2)   = vtemp(4*inod(n)-2)   + ke1(i+12,n+4)      
                    vtemp(4*inod(n)-1)   = vtemp(4*inod(n)-1)   + ke1(i+12,n+8)    
                    vtemp(4*inod(n)  )   = vtemp(4*inod(n)  )   + ke1(i+12,n+12)
                 endif  

              enddo           
!
! Store list of unique columns used
!
              do n=1,4
                 ineednode = 1 
                 do innz=1,nnzV
                    if (cols(innz).eq.inod(n)) ineednode = 0
                 enddo
       
                 if (ineednode.eq.1) then
                    nnzV = nnzV+1
                    cols(nnzV) = inod(n)
                 endif
              enddo  
       
           enddo           ! loop over elements

           if(vtemp  (4*inod1  ).eq.0.d0)then  !if inod1 is lied in the air space
                 vtemp  (4*inod1  ) = (1.d0,0.d0)
           endif 

           if(nodlab(inod1).eq.2.or.nodlab(inod1).eq.5) then  ! node lied on the top boundary
                 vtempey(4*inod1-2) = (1.d0,0.d0)
                 vtempez(4*inod1-1) = (1.d0,0.d0)
                 vtemp  (4*inod1  ) = (1.d0,0.d0)
           endif

        endif          
!        
! Apply boundary conditions for inod1 (note this is not needed since inprod_stored
! applies this as a boundary condition
!       
        if (nodlab(inod1).eq.1.or.nodlab(inod1).eq.3) then  ! if inod1 is lied on the top surface
           vtempex(4*inod1-3) = (1.d0,0.d0)
           vtempey(4*inod1-2) = (1.d0,0.d0)
           vtempez(4*inod1-1) = (1.d0,0.d0)
           vtemp  (4*inod1  ) = (1.d0,0.d0)
           nnzV=1
           cols(nnzV)=inod1
           if (nodlab(inod1).eq.1)then
            rhs1(4*inod1-3)    = (1.d0,0.d0)   ! ! Apply boundary conditions for inod1 lied on the top surface Ax
!           rhs2(4*inod1-2)    = (1.d0,0.d0)   ! ! Apply boundary conditions for inod1 lied on the top surface Ay
           endif
        endif
!
! Store  non zero elements in a compressed sparse row matrix:
!
        irw(4*inod1-3)   = nnz+1

           nnz=nnz+1
           val(nnz)  = vtempex(4*inod1-3)
           col(nnz)  = 4*inod1-3

        if(nodlab(inod1).ne.1.and.nodlab(inod1).ne.3)then
           do i1 = 1,nnzV 
              if (nodlab(cols(i1)).ne.1) then   !
                 if (vtempex(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1-3)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif 
                 if (vtempex(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1-3)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif   
                 if (vtempex(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1-3)) then     !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1)-1)
                    col(nnz) = 4*cols(i1)-1
                 endif
                 if (vtempex(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1-3)) then     !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
                 rhs1(4*inod1-3)    = rhs1(4*inod1-3)-1.d0*vtempex(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
!                 rhs2(4*inod1-3)    = rhs2(4*inod1-3)-1.d0*vtempex(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface
              endif

           enddo
        endif

        irw(4*inod1-2)   = nnz+1

           nnz=nnz+1 
           val(nnz)  = vtempey(4*inod1-2)
           col(nnz)  = 4*inod1-2

        if(nodlab(inod1).eq.4)then ! 
           do i1 = 1,nnzV 
              if (nodlab(cols(i1)).ne.1) then  !
                 if (vtempey(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1-2)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif  
                 if (vtempey(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1-2)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif
                 if (vtempey(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1-2)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1)-1)
                    col(nnz) =4*cols(i1)-1
                 endif
                 if (vtempey(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1-2)) then     !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
                 rhs1(4*inod1-2)    = rhs1(4*inod1-2)-1.d0*vtempey(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
!                 rhs2(4*inod1-2)    = rhs2(4*inod1-2)-1.d0*vtempey(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface
              endif

           enddo
        endif

        irw(4*inod1-1)   = nnz+1

           nnz=nnz+1
           val(nnz)  = vtempez(4*inod1-1)
           col(nnz)  = 4*inod1-1

        if(nodlab(inod1).eq.4)then ! 
           do i1 = 1,nnzV  
              if (nodlab(cols(i1)).ne.1)then    !
                 if (vtempez(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1-1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif 
                 if (vtempez(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1-1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif  
                 if (vtempez(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1-1)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1)-1)
                    col(nnz) =4*cols(i1)-1
                 endif
                 if (vtempez(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1-1)) then    !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
                 rhs1(4*inod1-1)    = rhs1(4*inod1-1)-1.d0*vtempez(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
!                 rhs2(4*inod1-1)    = rhs2(4*inod1-1)-1.d0*vtempez(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface
              endif

           enddo
        endif

        irw(4*inod1)   = nnz+1

           nnz=nnz+1
           val(nnz)  = vtemp(4*inod1)
           col(nnz)  = 4*inod1

        if(nodlab(inod1).eq.4)then  !
           do i1 = 1,nnzV 
              if (nodlab(cols(i1)).ne.1) then    ! 
                 if (vtemp(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif 
                 if (vtemp(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif 
                 if (vtemp(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1)) then   !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1)-1)
                    col(nnz) =4*cols(i1)-1
                 endif
                 if (vtemp(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1)) then    !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
                 rhs1(4*inod1)    = rhs1(4*inod1)-1.d0*vtemp(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
!                 rhs2(4*inod1)    = rhs2(4*inod1)-1.d0*vtemp(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface                 
              endif

           enddo
        endif

!	write(11,*) nnz 
    
        enddo        ! loop over nodes


!       close(11)
!
! Finally record the total number of non-zeros:
!
        irw(4*nnod+1) = nnz+1
! 
!        open(12,file='xxx.dat')
!        do i=1,nnz
!           write(12,*)val(i),col(i)
!        enddo

!        do i=1,4*nnod+1
!           write(12,*)irw(i)
!        enddo

!        close(12)  
!
! Deallocate arrays:
!
!        deallocate (indinod1)
    
        end subroutine gen_diago_lu_ex      
!
!======================================================================!
!=========================================================! gen_diago_lu
!======================================================================!
!
        subroutine gen_diago_lu_ey
!
! Compressed Row storage format:
!
!  val - nnz x 1 array of non-zero values of matrix, stored in order of 
!        the rows from 1 to nrows
!  col - nnz x 1 array of column number of each non-zero entry 
!        (doesn't have to be in order...)
!  irw - nnod+1 x 1 array, ith entry is points to first entry of the ith row in val.
!        last entry is equal to nnz.
!
! Modules used:
!
        use fem_params
        use lhs_mod
        use mesh
        use datinfo
        use csemarrays

        implicit none
    
        integer        :: i,i1,jj,n,inod(4),inod1,eps(8)
        real(8)        :: xe(4),ye(4),ze(4),del    
        integer, dimension(200) :: cols
        integer ::iele, nnzV, innz, ineednode
        complex(8)    :: ommu,k0,ke1(16,16)
        real(double)  ::omega,sige(6)
        common /kxw/omega

    
        ommu=ic*omega*mu0
        val  = (0.d0,0.d0)
        col = 0
        irw = 0
        nnz = 0   
!        rhs1 = (0.d0,0.d0) 
        rhs2 = (0.d0,0.d0)
! 
! Loop over nodes, 
!
!        open(11,file='nnz.dat')

        do inod1 = 1,nnod ! note there are 4*nnod unknowns (Ax Ay Az and phi)         
           vtempEX = (0.d0,0.d0)
           vtempEY = (0.d0,0.d0)
           vtempEZ = (0.d0,0.d0)
           vtemp   = (0.d0,0.d0)

           nnzV = 0
           cols = 0
!                     
! Loop over elements containing inod1
!  
         if (nodlab(inod1).ne.1.and.nodlab(inod1).ne.3) then   ! inod1 that doesnot lied on the top boundary

            do iele =1, nelein(inod1)           
              e = nodele(inod1,iele)      ! get element number         
              i = indinod1(inod1,iele)
              inod(1) = emap(1,e)
              inod(2) = emap(2,e)
              inod(3) = emap(3,e)
              inod(4) = emap(4,e)
              jj=irg(e)
!
!  Yank out the (y,z) coordinates for the current element 
!        
              do n=1,4
                 xe(n)=xn(inod(n))
                 ye(n)=yn(inod(n))
                 ze(n)=zn(inod(n))
              enddo

              sige=sg(jj,:)

              call gen_K1(xe,ye,ze,del,ommu,sige,Ke1) 
!	                     
!  Wipe out the row and column in the mass and stiffness 
!  matrices for those elements which lie in the air
!  of the FE mesh since a phi value is zero there.
!
              do n=1,4
                 if(nodlab(inod(n)).eq.5.or.nodlab(inod(n)).eq.2.or.nodlab(inod(n)).eq.1) then  ! node lied on the top boundary
                     ke1(:,n)=(0.d0,0.d0)
                     ke1(:,n+8)=(0.d0,0.d0)
                     ke1(:,n+12)=(0.d0,0.d0)
                 elseif (nodlab(inod(n)).eq.3) then                   
                     ke1(:,n)=(0.d0,0.d0)
                     ke1(:,n+4)=(0.d0,0.d0)
                     ke1(:,n+8)=(0.d0,0.d0)
                     ke1(:,n+12)=(0.d0,0.d0)
                 endif
              enddo

                 if(jj.eq.1)then !node doesnot lied on the top boundary, but in the air 
                    do n=1,4
                    ke1(:,n+12)=(0.d0,0.d0)
                    enddo
                 endif    
!                                  
!! 
!! Append parts corresponding to inod1's E and H vectors to a
!! temporary vector of nnod*2 length:        
!!                                    
              do n=1,4
                 vtempex(4*inod(n)-3) = vtempex(4*inod(n)-3) + ke1(i,n)           
                 vtempex(4*inod(n)-2) = vtempex(4*inod(n)-2) + ke1(i,n+4)      
                 vtempex(4*inod(n)-1) = vtempex(4*inod(n)-1) + ke1(i,n+8)    
                 vtempex(4*inod(n)  ) = vtempex(4*inod(n)  ) + ke1(i,n+12) 
                                                                      
                 vtempey(4*inod(n)-3) = vtempey(4*inod(n)-3) + ke1(i+4,n)           
                 vtempey(4*inod(n)-2) = vtempey(4*inod(n)-2) + ke1(i+4,n+4)      
                 vtempey(4*inod(n)-1) = vtempey(4*inod(n)-1) + ke1(i+4,n+8)    
                 vtempey(4*inod(n)  ) = vtempey(4*inod(n)  ) + ke1(i+4,n+12)
                                                                     
                 vtempez(4*inod(n)-3) = vtempez(4*inod(n)-3) + ke1(i+8,n)           
                 vtempez(4*inod(n)-2) = vtempez(4*inod(n)-2) + ke1(i+8,n+4)      
                 vtempez(4*inod(n)-1) = vtempez(4*inod(n)-1) + ke1(i+8,n+8)   
                 vtempez(4*inod(n)  ) = vtempez(4*inod(n)  ) + ke1(i+8,n+12)

                 if(jj.ne.1)then ! as the element lied in the air, the phi value is set to zero
                    vtemp(4*inod(n)-3)   = vtemp(4*inod(n)-3)   + ke1(i+12,n)           
                    vtemp(4*inod(n)-2)   = vtemp(4*inod(n)-2)   + ke1(i+12,n+4)      
                    vtemp(4*inod(n)-1)   = vtemp(4*inod(n)-1)   + ke1(i+12,n+8)    
                    vtemp(4*inod(n)  )   = vtemp(4*inod(n)  )   + ke1(i+12,n+12)
                 endif  

              enddo           
!
! Store list of unique columns used
!
              do n=1,4
                 ineednode = 1 
                 do innz=1,nnzV
                    if (cols(innz).eq.inod(n)) ineednode = 0
                 enddo
       
                 if (ineednode.eq.1) then
                    nnzV = nnzV+1
                    cols(nnzV) = inod(n)
                 endif
              enddo  
       
           enddo           ! loop over elements


           if(vtemp  (4*inod1  ).eq.0.d0)then  !if inod1 is lied in the air space
              vtemp  (4*inod1  ) = (1.d0,0.d0)
           endif 

           if(nodlab(inod1).eq.5.or.nodlab(inod1).eq.2) then  ! node lied on the top boundary
                 vtempex(4*inod1-3) = (1.d0,0.d0)
                 vtempez(4*inod1-1) = (1.d0,0.d0)
                 vtemp  (4*inod1  ) = (1.d0,0.d0)
           endif


        endif          
!        
! Apply boundary conditions for inod1 (note this is not needed since inprod_stored
! applies this as a boundary condition
!       
        if (nodlab(inod1).eq.1.or.nodlab(inod1).eq.3) then  ! if inod1 is lied on the top surface
           vtempex(4*inod1-3) = (1.d0,0.d0)
           vtempey(4*inod1-2) = (1.d0,0.d0)
           vtempez(4*inod1-1) = (1.d0,0.d0)
           vtemp  (4*inod1  ) = (1.d0,0.d0)
           nnzV=1
           cols(nnzV)=inod1
           if (nodlab(inod1).eq.1)then
            rhs2(4*inod1-2)    = (1.d0,0.d0)   ! ! Apply boundary conditions for inod1 lied on the top surface Ay
           endif
        endif
!
! Store  non zero elements in a compressed sparse row matrix:
!
        irw(4*inod1-3)   = nnz+1

        nnz=nnz+1
        val(nnz)  = vtempex(4*inod1-3)
        col(nnz)  = 4*inod1-3

        if(nodlab(inod1).eq.4)then  !
           do i1 = 1,nnzV 
              if (nodlab(cols(i1)).ne.1)then   !
                 if (vtempex(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1-3)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif 
                 if (vtempex(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1-3)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif   
                 if (vtempex(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1-3)) then     !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1)-1)
                    col(nnz) = 4*cols(i1)-1
                 endif
                 if (vtempex(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1-3)) then     !
                    nnz = nnz+1
                    val(nnz) = vtempex(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
!                 rhs1(4*inod1-3)    = rhs1(4*inod1-3)-1.d0/(ic*omega)*vtempex(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
                 rhs2(4*inod1-3)    = rhs2(4*inod1-3)-1.d0*vtempex(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface
              endif

           enddo
        endif

        irw(4*inod1-2)   = nnz+1

        nnz=nnz+1 
        val(nnz)  = vtempey(4*inod1-2)
        col(nnz)  = 4*inod1-2

        if(nodlab(inod1).ne.1.and.nodlab(inod1).ne.3)then
           do i1 = 1,nnzV 
              if (nodlab(cols(i1)).ne.1)then   !
                 if (vtempey(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1-2)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif  
                 if (vtempey(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1-2)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif
                 if (vtempey(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1-2)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1)-1)
                    col(nnz) =4*cols(i1)-1
                 endif
                 if (vtempey(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1-2)) then     !
                    nnz = nnz+1
                    val(nnz) = vtempey(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
!                 rhs1(4*inod1-2)    = rhs1(4*inod1-2)-1.d0/(ic*omega)*vtempey(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
                 rhs2(4*inod1-2)    = rhs2(4*inod1-2)-1.d0*vtempey(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface
              endif

           enddo
        endif

        irw(4*inod1-1)   = nnz+1

        nnz=nnz+1
        val(nnz)  = vtempez(4*inod1-1)
        col(nnz)  = 4*inod1-1

        if(nodlab(inod1).eq.4)then   !
           do i1 = 1,nnzV  
              if (nodlab(cols(i1)).ne.1)then   !
                 if (vtempez(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1-1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif 
                 if (vtempez(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1-1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif  
                 if (vtempez(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1-1)) then   !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1)-1)
                    col(nnz) =4*cols(i1)-1
                 endif
                 if (vtempez(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1-1)) then    !
                    nnz = nnz+1
                    val(nnz) = vtempez(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
!                 rhs1(4*inod1-1)    = rhs1(4*inod1-1)-1.d0/(ic*omega)*vtempez(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
                 rhs2(4*inod1-1)    = rhs2(4*inod1-1)-1.d0*vtempez(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface
              endif

           enddo
        endif

        irw(4*inod1)   = nnz+1

        nnz=nnz+1
        val(nnz)  = vtemp(4*inod1)
        col(nnz)  = 4*inod1

        if(nodlab(inod1).eq.4)then
           do i1 = 1,nnzV 
              if (nodlab(cols(i1)).ne.1)then   !
                 if (vtemp(4*cols(i1)-3).ne.0d0.and.(4*cols(i1)-3).gt.(4*inod1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1)-3)
                    col(nnz) =4*cols(i1)-3
                 endif 
                 if (vtemp(4*cols(i1)-2).ne.0d0.and.(4*cols(i1)-2).gt.(4*inod1)) then  !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1)-2)
                    col(nnz) =4*cols(i1)-2
                 endif 
                 if (vtemp(4*cols(i1)-1).ne.0d0.and.(4*cols(i1)-1).gt.(4*inod1)) then   !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1)-1)
                    col(nnz) =4*cols(i1)-1
                 endif
                 if (vtemp(4*cols(i1)).ne.0d0.and.(4*cols(i1)).gt.(4*inod1)) then    !
                    nnz = nnz+1
                    val(nnz) = vtemp(4*cols(i1))
                    col(nnz) = 4*cols(i1)
                 endif
              endif

              if (nodlab(cols(i1)).eq.1)then                    
!                 rhs1(4*inod1)    = rhs1(4*inod1)-1.d0/(ic*omega)*vtemp(4*cols(i1)-3)  ! ! Apply boundary conditions for inod1 lied on the top surface
                 rhs2(4*inod1)    = rhs2(4*inod1)-1.d0*vtemp(4*cols(i1)-2)  ! ! Apply boundary conditions for inod1 lied on the top surface                 
              endif

           enddo
        endif

!	write(11,*) nnz 
    
        enddo        ! loop over nodes

!       close(11)
!
! Finally record the total number of non-zeros:
!
        irw(4*nnod+1) = nnz+1
! 
!        open(12,file='xxx.dat')
!        do i=1,nnz
!           write(12,*)val(i),col(i)
!        enddo

!        do i=1,4*nnod+1
!           write(12,*)irw(i)
!        enddo

!        close(12)  
!
! Deallocate arrays:
!
    
        end subroutine gen_diago_lu_ey 
!      
!======================================================================!
!=========================================================! readrunfile
!======================================================================!
!
        subroutine readrunfile
!
!  Read in the RUNFILE for starting the 3D MT Modeling
!
!======================================================================!
        use    	datinfo        
        implicit none
        integer   i,itemp,lerr,j
        character(56)  vers,ctemp
        integer, dimension(:), allocatable 	::	temp
        real(8)       :: rho1,lalphx,lalphy,lalphz
        real(8)       :: sigc(3,3),sigt(3,3),rx(3,3),ry(3,3),rz(3,3),rxyz(3,3)
!
!  Open RUNFILE and read in modeling info:
!
	open (unit=10,file='RUNFILE',status='old',iostat=lerr)
	if (lerr .ne. 0) then
		write(*,*) ' Error opening RUNFILE'
		stop
	end if

        read(unit=10,fmt='(18x,a)')    vers
        if  (vers(1:16) .ne. 'MT3Dforward')  then
    	    write(*,*) 'Error: Runfile format ',vers
    	    write(*,*) 'not supported, stopping.'
    	    write(*,*) 'Try using MT3Dforward'
    	    stop
        endif

        write(6,*) '======== Reading in "RUNFILE" ',&
                    & '============================================'
        write(*,*) ' '

 	read(unit=10,fmt='(18x,a)')    	tristr
	write(6,fmt='(a24,a32)') 'TetGen Cmd:', trim(tristr)

	read(unit=10,fmt='(18x,a)')    	modelroot
	write(6,fmt='(a24,a32)') 'Modelroot:', trim(modelroot)

	read(unit=10,fmt='(18x,i4)')    strtmodnum
        write(6,fmt='(a24,i32)') 'Starting Model #:', (strtmodnum)

	read(unit=10,fmt='(18x,i4)')  maxadapt
	write(6,fmt='(a24,i32)') 'Max # refinements:',maxadapt

	read(unit=10,fmt='(18x,a56)')ctemp
	read(ctemp,*) pctrefine
        write(6,fmt='(a24,f32.3)') 'Percent to refine:',pctrefine
        pctrefine = pctrefine/100d0

        read(unit=10,fmt='(18x,a56)')ctemp
	read(ctemp,*) tolerance
	write(*,fmt='(a24,f32.3)') 'Stopping Tolerance (%):',tolerance

	read(unit=10,fmt='(18x,a56)')ctemp
	read(ctemp,*)  nro
        write(6,fmt='(a24,i32)') '# resistivity blocks:',nro
        allocate(sg(nro,6),roxx(nro),royy(nro),rozz(nro),alphx(nro),alphy(nro),alphz(nro))
        write(6,fmt='(a12,5x,a12)') 'Block No:','Resistivity' 

        DO  I=1,NRO
            READ(10,*)j,roxx(i),royy(i),rozz(i),alphx(i),alphy(i),alphz(i)
            write(6,fmt='(i2,2x,6e12.3)')j,roxx(i),royy(i),rozz(i),alphx(i),alphy(i),alphz(i)

            lalphx=alphx(i)*pi/180.d0
            lalphy=alphy(i)*pi/180.d0
            lalphz=alphz(i)*pi/180.d0
           
            sigc(:,1)=(/1.d0/roxx(i), 0.d0, 0.d0/)
            sigc(:,2)=(/0.d0, 1.d0/royy(i), 0.d0/)
            sigc(:,3)=(/0.d0, 0.d0, 1.d0/rozz(i)/)

            rz(:,1)=(/dcos(lalphx),dsin(lalphx),0.d0/)
            rz(:,2)=(/-dsin(lalphx),dcos(lalphx),0.d0/)
            rz(:,3)=(/0.d0,0.d0,1.d0/)

            rx(:,1)=(/1.d0,0.d0,0.d0/)
            rx(:,2)=(/0.d0,dcos(lalphy),dsin(lalphy)/)
            rx(:,3)=(/0.d0,-dsin(lalphy),dcos(lalphy)/)

            ry(:,1)=(/dcos(lalphz),dsin(lalphz),0.d0/)
            ry(:,2)=(/-dsin(lalphz),dcos(lalphz),0.d0/)
            ry(:,3)=(/0.d0,0.d0,1.d0/)
      

            rxyz=matmul(matmul(rz,rx),ry)
            sigt=matmul(matmul(rxyz,sigc),transpose(rxyz))
            
            sg(i,1)=sigt(1,1)
            sg(i,2)=sigt(1,2)
            sg(i,3)=sigt(1,3)
            sg(i,4)=sigt(2,2)
            sg(i,5)=sigt(2,3)
            sg(i,6)=sigt(3,3)
!            write(6,fmt='(i2,2x,6e12.3)')i,sg(i,1),sg(i,2),sg(i,3),sg(i,4),sg(i,5),sg(i,6)
        enddo      
!
! Read in the requested frequencies:
!
        nfrqrefine = 0
        read(unit=10,fmt='(18x,a56)')ctemp
        read(ctemp,*) nfreq
        write(*,*)  ' '
        write(6,fmt='(a24,a32)') 'Frequencies:', '[Hz]:'
        allocate (freqvec(nfreq), temp(nfreq))

        do i  = 1,nfreq
            read(10,*) freqvec(i), temp(i)
            write(*,fmt='(i24,e32.4)') i, freqvec(i) 
            if (temp(i).eq.1)  nfrqrefine  = nfrqrefine+1        
        enddo

        allocate ( frqrefine(nfrqrefine) )
        nfrqrefine = 0
        write(6,fmt='(a24,a32)') 'Frequencies to refine:','[Hz]' 
   
        do i  = 1,nfreq    
            if (temp(i).eq.1) then
                nfrqrefine =nfrqrefine+1
                frqrefine(nfrqrefine) = freqvec(i)  
                write(6,fmt='(i24,e32.3)') nfrqrefine, frqrefine(nfrqrefine)
             endif
        enddo   
! 
!
! Read in the Site locations:
!
        read(unit=10,fmt='(18x,a56)')ctemp
        read(ctemp,*) nsite
        write(*,*) ' ' 
        write(*,fmt = '(a24,i6)') ' # Receivers:',nsite
        write(*,68) 'Receiver #', 'x[m]','y[m]','z[m]'
    
        allocate(xsite(nsite),ysite(nsite), zsite(nsite))
        do i = 1,nsite
           read(10,*) xsite(i),ysite(i),zsite(i)
           write(*,69) i, xsite(i),ysite(i),zsite(i)   
        enddo   
!
! Close RUNFILE
!
        close(10)
        write(*,*) ' '
        write(*,*) '====================== Forward computation starts ======================'  
	write(*,*) ' '

68      format(a24,3(1x,a12)) 
69      format(i24,3(1x,f12.1))

	end subroutine readrunfile
!
!======================================================================!
!=====================================================! UPDATEFILENAMES
!======================================================================!
!
	subroutine updatefilenames(iadapt)
!
! Upadates some strings to point to the current finite element grid.
!
!======================================================================!
	use datinfo     
	implicit none

	integer :: iTx,iadapt
	character(256) :: str,str1
	character(2)   :: cadapt,cTx,clast
!
!  Make character strings for current transmitter, frequency and mesh
!
	write(cadapt,'(i2)') iadapt-1+strtmodnum
	write(clast,'(i2)') iadapt-1+strtmodnum-1
!	write(cTx,'(i2)') iTx

	cadapt = adjustl(cadapt)
	clast = adjustl(clast)
!	cTx    = adjustl(cTx)

!        adaptroot   = trim(modelroot)//'.'//trim(cTx)
	currentname = trim(modelroot)//'.'//trim(cadapt)
        lastname    = trim(modelroot)//'.'//trim(clast)

	elementfile = trim(currentname)//'.ele'
	nodefile    = trim(currentname)//'.node'
!
! copy over initial starting model from triangle if first run:
!
 	if ((iadapt.eq.1).and.(strtmodnum.eq.1)) then
!
! copy file output files from triangle into new numbering scheme
!
 		str1 = 'cp '//trim(modelroot)//'.1'//'.ele '
 		str = trim(str1)//' '//elementfile
 		call system(str)
 		str1 = 'cp '//trim(modelroot)//'.1'//'.node '
 		str = trim(str1)//' '//nodefile
 		call system(str)
! 		str1 = 'cp '//trim(modelroot)//'.1'//'.poly '
! 		str = trim(str1)//' '//trim(currentname)//'.poly'               
! 		call system(str)
 	endif

	end subroutine updatefilenames
!
!======================================================================!
!============================================================! READMESH
!======================================================================!
!
        subroutine readmesh
!
!  Reads in 3D Tetrahedral grid in the format used by TetGen, the
!  open source C-code by Si, H., 2003. 
!   TETGEN: A 3D delaunay tetrahedral mesh generator. 
!   http://tetgen.berlios.de
!
	use        fem_params
	use        mesh
	use        datinfo
	implicit none
	integer       :: natt,n,err,i

	write(*,*) ' '
	write(6,*) '-------- Reading in Finite Element Grid Files ',&
                      & '-------------------'
        write(*,*) ' '
!
!  Read in NODE file
!
	open(10,file=trim(nodefile),status='old',iostat=err)

	if (err .ne. 0) then
		write(*,*) ' Error opening NODE file: '//trim(nodefile)
		stop
	end if

        write(6,fmt='(a24,a32)') 'Node Filename:',trim(nodefile)
        read (10,*,end=188, err=189) nnod,iskip,natt,iskip
        write(6,fmt='(a24,i32)') '# Nodes:',  nnod

        allocate (xn(nnod),yn(nnod),zn(nnod),nodlab(nnod),bcnod(nnod) )

        if (natt.eq.0) then
           do 12 n=1,nnod
               read(10,*) iskip, xn(n),yn(n),zn(n)
 12        continue
        else
           do 13 n=1,nnod
               read(10,*,end=188, err=189) iskip, xn(n),yn(n),zn(n),(skip,i=1,natt)
 13        continue
        endif

        close(10)
!
!  Read in ELEMENT file
!  Read in the element map where emap(i,e) contains the index of the
!  i'th node (i=1,3) which defines the e'th element (e=1,nele).
!
        open(unit=11,file=trim(ELEMENTFILE),status='old',iostat =err)

        if (err .ne. 0) then
            write(*,*) ' Error opening ELEMENT file: '//trim(ELEMENTFILE)
            stop
        end if

        write(6,fmt='(a24,a32)') 'Element Filename:',trim(ELEMENTFILE)
        read (11,*,end=178, err=179) nele,iskip,natt
        write(6,fmt='(a24,i32)') '# Elements:',  nele

        allocate (emap(4,nele), irg(nele) )

        do 20 e=1,nele
        read(11,*,end=178, err=179) iskip, (emap(i,e),i=1,4), irg(e)
 20     continue

        close(11)
        write(6,*) ' '

        return

178     write(*,*) ' ELEMENT file ended prematurely'
179     write(*,*) ' Error reading ELEMENT file'
188     write(*,*) ' NODE file ended prematurely'
189     write(*,*) ' Error reading NODE file'

        end subroutine readmesh  
!
!======================================================================! 
!======================================================================! 
!====================================================== errorindicator !
!======================================================================!
!
        subroutine errorindicator_land(rhs0)
!              							       !
!  Written by: Yuguo Li						       !
!              Scripps Institution of Oceanography 		       ! 
!  Modified by: Yixin Ye
!              East china University of Technology
!  Styled after errorindicator for marine 2D CSEM by KWK.                         
!======================================================================!   
        use             datinfo         
        use             fem_params
        use             mesh  
        use             csemarrays  ! inputs rhs and outputs SmQh grads, error estimates
        implicit none
    
        integer                            :: i,i2,idotpos
!       real(double),    dimension(:)      :: errnrmx,errnrmy,errnrmz
        complex(double), dimension(:)      :: dual,egradx,egrady,egradz,dualgx,dualgy,dualgz
        complex(double), dimension(:,:)    :: esmgradx,esmgrady,esmgradz,dualsmgx,dualsmgy,dualsmgz
        integer, dimension(:)              :: irele
        integer                            :: nregion
        complex(double), dimension(4*nnod) :: rhs0
!       logical			           :: lrefine ! flag for refinement iteration
    
        allocatable irele    !,errnrmx,errnrmy,errnrmz
        allocatable dual,egradx,egrady,egradz,esmgradx,esmgrady,esmgradz,dualgx,dualgy,dualgz,dualsmgx,dualsmgy,dualsmgz
!
!  Allocate local arrays:
!
        allocate (irele(nele))
        allocate (egradx(4*nele),egrady(4*nele),egradz(4*nele))
        allocate (esmgradx(4,4*nele), esmgrady(4,4*nele),esmgradz(4,4*nele)) 

        errnrm_eh    = 0d0
!        egradnrm_eh  = 0d0

        esmgradx  = (0d0,0d0)
        esmgrady  = (0d0,0d0)
        esmgradz  = (0d0,0d0)

        egradx  = (0d0,0d0)
        egrady  = (0d0,0d0)
        egradz  = (0d0,0d0) 
!
!  Get piecewise (PCWS) constant gradients of rhs1 for each element.  area of elements 
!  also returned.
!
        call getgradyz(rhs0,egradx,egrady,egradz,area)
!
!  Get regions of constant conductivity
!
        call getregions(nregion,irele)      
!
!  Compute L2 projection of gradient, smooth if necessary
!
        call SmQhregions_land(egradx,nregion,irele,esmgradx,area) 
        call SmQhregions_land(egrady,nregion,irele,esmgrady,area) 
        call SmQhregions_land(egradz,nregion,irele,esmgradz,area)  
!
!  Compute Error indicator using normed difference between 
!
        call l2norm_land(egradx,egrady,egradz,esmgradx,esmgrady,esmgradz,errnrm_eh,area)

!       open(16,file='element.error',status='REPLACE')
!       do e=1,nele
!          write(16,fmt = '(3(e15.8,2x))') errnrm_eh(3*e-2),errnrm_eh(3*e-1),errnrm_eh(3*e)
!       enddo
!       close(16)
!  
!       open(11,file='egrady.dat')
!	write(11,*) 'sd_n=',2*nele
!       do i = 1 , 2*nele
!	   write(11,'(a,i,3(e16.8,e16.8))') 'i=', i , egrady(i)
!       enddo
!       close(11)
!
!       open(11,file='egradz.dat')
!	write(11,*) 'sd_n=',2*nele
!       do i = 1 , 2*nele
!	   write(11,'(a,i,3(e16.8,e16.8))') 'i=', i , egradz(i)
!       enddo
!       close(11)
!
!       open(11,file='esmgrady.dat')
!	write(11,*) 'sd_n=',2*nele
!       do i = 1 , 2*nele
!	   write(11,'(a,i,3(e16.7,e16.7))') 'i=', i , esmgrady(1:3,i)
!       enddo
!       close(11)
!
!       open(11,file='esmgradz.dat')
!	write(11,*) 'sd_n=',2*nele
!       do i = 1 , 2*nele
!	   write(11,'(a,i,3(e16.7,e16.7))') 'i=', i , esmgradz(1:3,i)
!       enddo
!       close(11)
!
!  Compute weights for error by solving dual problem
!
        if (weighterr.eq.1) then ! compute dual weights

           allocate (dualgx(4*nele),dualgy(4*nele),dualgz(4*nele),dual(4*nnod))
           allocate (dualsmgx(4,4*nele), dualsmgy(4,4*nele),dualsmgz(4,4*nele)) 
!
!  Solving dual problem for error weighting function 
!
           call solve_dual(egradx,egrady,egradz,esmgradx,esmgrady,esmgradz,dual,area)
!
!          do i = 1 , 2*nnod
!	      write(*,'(a,i,e16.8,e16.8)') 'i=',i,dual(i)
!          enddo
!
!          open(11,file='dual.dat')
!	   write(11,*) 'sd_n=',2*nnod
!          do i = 1 , 2*nnod
!	      write(11,'(a,i,e16.8,e16.8)') 'i=',i,dual(i)
!          enddo
!          close(11)
!          stop
!   
!  Get piecewise (PCWS) constant gradients of rhs1 for each element.  area of elements 
!  also returned.
!
           call getgradyz(dual,dualgx,dualgy,dualgz,area)
!
!  Compute L2 projection of gradient, smooth if necessary
!
           call SmQhregions_land(dualgx,nregion,irele,dualsmgx,area)
           call SmQhregions_land(dualgy,nregion,irele,dualsmgy,area)   
           call SmQhregions_land(dualgz,nregion,irele,dualsmgz,area)
!             
	   call l2norm_dual_land(dualgx,dualgy,dualgz,dualsmgx,dualsmgy,dualsmgz,&
		                      errnrm_eh,area)   
 
        endif  ! if WEIGHTERR.eq.1
!
!  Deallocate arrays
!
        deallocate (irele)
        deallocate (egradx,egrady,egradz) !,errnrmx,errnrmy,errnrmz
        deallocate (esmgradx,esmgrady,esmgradz ) 
   
        if(weighterr.eq.1)then
           deallocate (dual,dualgx,dualgy,dualgz) 
           deallocate (dualsmgx,dualsmgy,dualsmgz)
        endif  

        end subroutine errorindicator_land
!    
!======================================================================! 
!======================================================! GEN_BCNOD_DUAL
!======================================================================! 
!
        subroutine gen_bcnod_dual

        use           fem_params
        use           mesh
        integer       :: i

        bcnod = 0
!
!  Remap the node labels according to the problem that is being
!  solved...
!    
        do i=1,nnod
           if (nodlab(i).ne.4) bcnod(i) = 1
        enddo
!
        end subroutine gen_bcnod_dual
!   
!======================================================================! 
!======================================================== write_errnrm !
!======================================================================!     
!
        subroutine write_errnrm(errnrm,area)      
!
! Flags elements for refinement, writes out the normed error estimate 
! and writes out the .area file used by Triangle for grid refinement
!
        use             mesh      
        use             datinfo         
        use             fem_params
        implicit none
    
        integer                            :: i,idotpos,nrefine
        real(double),    dimension(nele)   :: area,errnrm
        integer, dimension(:),allocatable  :: irele
        logical, dimension(:),allocatable  :: lrefined
        real(double)                       :: areamin, skndepth,sknarea

        allocate(irele(nele),lrefined(nele))   
!    
!  Sort normed error from min to max:
!
        call indexx(nele,errnrm,irele) 
!
!  Create new area constraint for a small % of largest error elements:
!
        nrefine = 0
        lrefined = .false.
        i = nele+1
        do while ((i.gt.1).and.(nrefine.lt.pctrefine*nele))
           i= i-1
           e = irele(i)
           if (area(e)*0.5d0.gt.minelearea) then
              nrefine = nrefine+1
              lrefined(e) = .true.
              area(e) = area(e)*0.5d0
           endif
        enddo
 
        if (nrefine.eq.0) then
    	   area(1:nele) = -1
        else
!
! make other area negative so no refinement:
!
    	   do e = 1,nele
      	      if (.not.lrefined(e)) area(e) = -1
    	   enddo
        endif

        deallocate (irele,lrefined)
!
! Save the total summed error for the current grid
!
        idotpos = scan((elementfile),achar(46),.true.)
        open(unit=13,file=trim(elementfile(1:idotpos))//'error',status='unknown') 
        do e=1,nele
           write(13,*) errnrm(e)  
        enddo
        close(13)

!
! save new area file:
!
        open(unit=13,file=trim(elementfile(1:idotpos))//'vol',status='unknown') 
        write(13,*) nele
        do e=1,nele
           write(13,fmt='(i8,2x,f20.3)' ) e,area(e)
        enddo
        close(13)   
        
        end subroutine write_errnrm     
!    
!======================================================================! 
!========================================================== SOLVE_DUAL !
!======================================================================! 
!
        subroutine solve_dual(egradx,egrady,egradz,esmgradx,esmgrady,esmgradz,dual,area)  
!======================================================================! 
!                                                                      !
!  Written by: Yuguo Li , Modified by: Yixin Ye                        ! 
!                                                                      ! 
!  Styled after solve_etws for 2D marine CSEM  by KWK.                  
!======================================================================!                                            

        use fem_params
        use mesh
        use datinfo
        use lhs_mod  
        implicit none
    
        integer                              :: i,isite,ep,jj,inod1,iele
        integer, dimension(4)                :: n
        real(8)                              :: wt,sige,xp,yp,zp
        real(double), dimension(4)           :: a,b,c,xe,ye,ze
        complex(double), dimension(4)        :: gx,gy,gz,gx1,gy1,gz1,gx2,gy2,gz2,gx3,gy3,gz3
        complex(double), dimension(4*nele)   :: egradx,egrady,egradz
        complex(double), dimension(4,4*nele) :: esmgradx,esmgrady,esmgradz
        complex(double), dimension(nnod*4)   :: dual
        real(double),dimension(nele) :: area
        real(8)     :: areamin, skndepth,sknarea
        real(double)  ::omega
        common /kxw/omega
!
!  Initialize the RHS1
!
        dual = (0.d0,0.d0)      
! 
!  Generate new boundary conditions 
!
        call gen_bcnod_dual
!
!  Make a list of the elements EM sites are located in:
!
        do isite=1,nsite

           call eget(xsite(isite),ysite(isite),zsite(isite),ep)

!           xp=xsite(isite)
!           yp=ysite(isite)
!           zp=zsite(isite)

!           do inod1=1,nnod           ! loop over nodes
!              if (xn(inod1).eq.xp.and.yn(inod1).eq.yp.and.zn(inod1).eq.zp) exit  ! get node number
!           enddo

!           do iele =1, nelein(inod1)          
!              ep = nodele(inod1,iele)     ! get element number
!
! Compute contribution to rhs1 from each element:
!
           do i=1,4
              n(i)  = emap(i,ep)
              xe(i) = xn(n(i))
              ye(i) = yn(n(i))
              ze(i) = zn(n(i))
           enddo
        
           a(1) = ye(4)*(ze(3)-ze(2)) + ye(2)*(ze(4)-ze(3)) + ye(3)*(ze(2)-ze(4))
           a(2) = ye(4)*(ze(1)-ze(3)) + ye(1)*(ze(3)-ze(4)) + ye(3)*(ze(4)-ze(1))
           a(3) = ye(4)*(ze(2)-ze(1)) + ye(1)*(ze(4)-ze(2)) + ye(2)*(ze(1)-ze(4))
           a(4) = ye(3)*(ze(1)-ze(2)) + ye(1)*(ze(2)-ze(3)) + ye(2)*(ze(3)-ze(1))
           b(1) = xe(3)*(ze(4)-ze(2)) + xe(4)*(ze(2)-ze(3)) + xe(2)*(ze(3)-ze(4))
           b(2) = xe(3)*(ze(1)-ze(4)) + xe(4)*(ze(3)-ze(1)) + xe(1)*(ze(4)-ze(3))
           b(3) = xe(2)*(ze(4)-ze(1)) + xe(4)*(ze(1)-ze(2)) + xe(1)*(ze(2)-ze(4))
           b(4) = xe(2)*(ze(1)-ze(3)) + xe(3)*(ze(2)-ze(1)) + xe(1)*(ze(3)-ze(2))
           c(1) = xe(4)*(ye(3)-ye(2)) + xe(2)*(ye(4)-ye(3)) + xe(3)*(ye(2)-ye(4))
           c(2) = xe(4)*(ye(1)-ye(3)) + xe(1)*(ye(3)-ye(4)) + xe(3)*(ye(4)-ye(1))
           c(3) = xe(4)*(ye(2)-ye(1)) + xe(1)*(ye(4)-ye(2)) + xe(2)*(ye(1)-ye(4))
           c(4) = xe(3)*(ye(1)-ye(2)) + xe(1)*(ye(2)-ye(3)) + xe(2)*(ye(3)-ye(1))
!
! Apply a penalty weighting binary term to very small elements...
!
!           jj=irg(ep)
!           sige=sg(jj,1)
!           skndepth = 500*dsqrt(2d0*PI / ( sige*omega ) ) ! skin depth in element e
!      	   areamin = 0.5*(skndepth/sknlimit)**2   ! approx min(area) limit parameter

!           if (area(ep).lt.areamin) then        
!              wt = (area(ep)/areamin)**2 ! quadratic down weighting
!           else
              wt = 1.d0
!           endif  
!       
! Approximate H1 seminorm:
!
           do i=1,4
              gx(i)  = wt*cdabs(esmgradx(i,4*ep-3)-egradx(4*ep-3) )
              gy(i)  = wt*cdabs(esmgrady(i,4*ep-3)-egrady(4*ep-3) )
              gz(i)  = wt*cdabs(esmgradz(i,4*ep-3)-egradz(4*ep-3) ) 
              gx1(i) = wt*cdabs(esmgradx(i,4*ep-2)-egradx(4*ep-2) )
              gy1(i) = wt*cdabs(esmgrady(i,4*ep-2)-egrady(4*ep-2) )
              gz1(i) = wt*cdabs(esmgradz(i,4*ep-2)-egradz(4*ep-2) )
              gx2(i) = wt*cdabs(esmgradx(i,4*ep-1)-egradx(4*ep-1) )
              gy2(i) = wt*cdabs(esmgrady(i,4*ep-1)-egrady(4*ep-1) )
              gz2(i) = wt*cdabs(esmgradz(i,4*ep-1)-egradz(4*ep-1) )
              gx3(i) = wt*cdabs(esmgradx(i,4*ep  )-egradx(4*ep  ) )
              gy3(i) = wt*cdabs(esmgrady(i,4*ep  )-egrady(4*ep  ) )
              gz3(i) = wt*cdabs(esmgradz(i,4*ep  )-egradz(4*ep  ) )
           enddo
!  Integral of (QhGradUh-Uh) * GradV, where v is the basis function 
  
           do i=1,4

              dual(4*n(i)-3) = dual(4*n(i)-3) +  &
              &  a(i)/24.d0*sum(gx(1:4)) + b(i)/24.d0*sum(gy(1:4)) + c(i)/24.d0*sum(gz(1:4))

              dual(4*n(i)-2) = dual(4*n(i)-2) +  &
              &  a(i)/24.d0*sum(gx1(1:4)) + b(i)/24.d0*sum(gy1(1:4)) + c(i)/24.d0*sum(gz1(1:4))

              dual(4*n(i)-1) = dual(4*n(i)-1) +  &
              &  a(i)/24.d0*sum(gx2(1:4)) + b(i)/24.d0*sum(gy2(1:4)) + c(i)/24.d0*sum(gz2(1:4))

              dual(4*n(i)) = dual(4*n(i)) +  &
              &  a(i)/24.d0*sum(gx3(1:4)) + b(i)/24.d0*sum(gy3(1:4)) + c(i)/24.d0*sum(gz3(1:4))
        
           enddo

!           enddo !loop over elements containing the node
        
        enddo  ! loop over sites
!
!  Set sourcing function to 0 on outer boundaries
!
        do i=1,nnod
           if (bcnod(i).ge.1) then
              dual(4*i-3) = 0d0
              dual(4*i-2) = 0d0
              dual(4*i-1) = 0d0
              dual(4*i) = 0d0
           endif
        enddo   
!
!  Solve dual problem for error weighting:
!
        call SSOR_PCG_CSR(4*nnod,NNZ,val,dual,col,irw)      

        end subroutine solve_dual
!      
!======================================================================! 
!================================================================! EGET
!======================================================================! 
!
        subroutine eget(xp,yp,zp,ep)

!======================================================================! 
!                                                                      !
!  Routine to extract the index of the element containing the point    !
!  (xp,yp,zp).                                                         !
!                                                                      !
!======================================================================! 
        use             fem_params
        use             mesh
        use             datinfo
        implicit none
        
        integer       :: i,jj,emin,ep,n(4)
        real(double)  :: xp,yp,zp,xe(4),ye(4),ze(4)
        real(double)  :: sigmin,sige,minxe,minye,minze,maxxe,maxye,maxze
        real(double)  :: v0,v1,v2,v3,v4,vmin
!
!  Loop over elements in the mesh
!
        sigmin = 1.d-99

        do e=1,nele
           jj=irg(e)
           sige=sg(jj,1)
!
!  Loop over nodes within the current element
!  
           do i=1,4
              n(i)=emap(i,e)
              xe(i) =   xn(n(i))
              ye(i) =   yn(n(i))
              ze(i) =   zn(n(i))
           enddo
!
!  Compute the element which contains the point xp, yp, zp
!
           maxxe=maxval(xe(1:4))
           minxe=minval(xe(1:4))
           maxye=maxval(ye(1:4))
           minye=minval(ye(1:4))
           maxze=maxval(ze(1:4))
           minze=minval(ze(1:4))

           if(xp.ge.minxe.and.xp.le.maxxe.and.yp.ge.minye.and.yp.le.maxye.and.zp.ge.minze.and.zp.le.maxze)then                  
                    v0=(xe(2)-xe(1))*(ye(3)*ze(4)-ze(3)*ye(4)) + (xe(1)- xe(3))*(ye(2)*ze(4)-ze(2)*ye(4)) + (xe(4)-xe(1))*(ye(2)*ze(3)-ze(2)*ye(3)) &
                     &  +(xe(3)- xe(2))*(ye(1)*ze(4)-ze(1)*ye(4)) + (xe(2)- xe(4))*(ye(1)*ze(3)-ze(1)*ye(3)) + (xe(4)- xe(3))*(ye(1)*ze(2)-ze(1)*ye(2)) 

                    v1=(xe(2)-xp)*(ye(3)*ze(4)-ze(3)*ye(4)) + (xp- xe(3))*(ye(2)*ze(4)-ze(2)*ye(4)) + (xe(4)-xp)*(ye(2)*ze(3)-ze(2)*ye(3)) &
                     &  +(xe(3)- xe(2))*(yp*ze(4)-zp*ye(4)) + (xe(2)- xe(4))*(yp*ze(3)-zp*ye(3)) + (xe(4)- xe(3))*(yp*ze(2)-zp*ye(2))    
                              
                    v2=(xp-xe(1))*(ye(3)*ze(4)-ze(3)*ye(4)) + (xe(1)- xe(3))*(yp*ze(4)-zp*ye(4)) + (xe(4)-xe(1))*(yp*ze(3)-zp*ye(3)) &
                     &  +(xe(3)- xp)*(ye(1)*ze(4)-ze(1)*ye(4)) + (xp- xe(4))*(ye(1)*ze(3)-ze(1)*ye(3)) + (xe(4)- xe(3))*(ye(1)*zp-ze(1)*yp)
                
                    v3=(xe(2)-xe(1))*(yp*ze(4)-zp*ye(4)) + (xe(1)- xp)*(ye(2)*ze(4)-ze(2)*ye(4)) + (xe(4)-xe(1))*(ye(2)*zp-ze(2)*yp) &
                     &  +(xp- xe(2))*(ye(1)*ze(4)-ze(1)*ye(4)) + (xe(2)- xe(4))*(ye(1)*zp-ze(1)*yp) + (xe(4)- xp)*(ye(1)*ze(2)-ze(1)*ye(2))
                     
                    v4=(xe(2)-xe(1))*(ye(3)*zp-ze(3)*yp) + (xe(1)- xe(3))*(ye(2)*zp-ze(2)*yp) + (xp-xe(1))*(ye(2)*ze(3)-ze(2)*ye(3)) &
                     &  +(xe(3)- xe(2))*(ye(1)*zp-ze(1)*yp) + (xe(2)- xp)*(ye(1)*ze(3)-ze(1)*ye(3)) + (xp- xe(3))*(ye(1)*ze(2)-ze(1)*ye(2))

                    if(v0*v4.ge.0.and.v0*v1.ge.0.and.v0*v2.ge.0.and.v0*v3.ge.0.and.sige.gt.sigmin)then            
                             sigmin=sige
                             ep = e
                    endif                                           	          
                                                      
           endif
!
!  Move on to the next element...
!
        enddo

        end subroutine eget
!        
!======================================================================! 
!========================================================= L2NORM_LAND !
!======================================================================! 
!
        subroutine l2norm_land(egradx,egrady,egradz,esmgradx,esmgrady,esmgradz,&
                           errnrm_eh,area)
!======================================================================! 
!                                                                      !
!  Written by: Yuguo Li                                                !
!              Scripps Institution of Oceanography                     !
!  Modified by: Yixin Ye                                               !
!              East china University of Technology                     ! 
!                                                                      !                       
!======================================================================!        
        use fem_params
        use mesh   
        use   datinfo
        implicit none
    
        integer                              :: i,kk,eleshft,jj
        complex(double), dimension(4*nele)   :: egradx, egrady, egradz ! element gradient
        complex(double), dimension(4,4*nele) :: esmgradx, esmgrady, esmgradz  
        complex(double), dimension(3,4)      :: gxyz
        real(double)                         :: wt(3),xi(5)
        real(double)                         :: s(3),f(3)
        real(double),dimension(4*nele)       :: errnrmx,errnrmy,errnrmz
        real(double),    dimension(nele)     :: area,errnrm_eh      
!
!  13 - Point Quadrature weights and points in area coordinates:
!
        DATA WT(1)/ 0.0131555555555D0/, WT(2)/0.0076222222222D0/,&
          WT(3)/ 0.0248888888888D0/

        DATA XI(1)/ 0.2500000000000D0/, XI(2)/0.785714285714286D0/,&
          XI(3)/ 0.0714285714285714D0/, XI(4)/0.399403576166799D0/,&
          XI(5)/ 0.100596423833201D0/
    
        eleshft = 3
      
        do kk = 1,4
           do e = 1,nele  
              s=0.d0          
              do i=1,4
                 gxyZ(1,I)  = esmgradx(i,4*e-eleshft)-egradx(4*e-eleshft)
                 gxyZ(2,I)  = esmgrady(i,4*e-eleshft)-egrady(4*e-eleshft)
                 gxYZ(3,I)  = esmgradz(i,4*e-eleshft)-egradz(4*e-eleshft)
              enddo
        
              DO I=1,3
                 F(1) = (CDABS(XI(1)*(GXYZ(I,1)+GXYZ(I,2)+GXYZ(I,3)+GXYZ(I,4))))**2

                 F(2) = (CDABS(XI(2)*GXYZ(I,1)+XI(3)*(GXYZ(I,2)+GXYZ(I,3)+GXYZ(I,4))))**2 &
                       +(CDABS(XI(2)*GXYZ(I,2)+XI(3)*(GXYZ(I,1)+GXYZ(I,3)+GXYZ(I,4))))**2 &
                       +(CDABS(XI(2)*GXYZ(I,3)+XI(3)*(GXYZ(I,1)+GXYZ(I,2)+GXYZ(I,4))))**2 &
                       +(CDABS(XI(2)*GXYZ(I,4)+XI(3)*(GXYZ(I,1)+GXYZ(I,2)+GXYZ(I,3))))**2 

                 F(3) = (CDABS(XI(4)*(GXYZ(I,1)+GXYZ(I,2))+XI(5)*(GXYZ(I,3)+GXYZ(I,4))))**2 &
                       +(CDABS(XI(4)*(GXYZ(I,1)+GXYZ(I,3))+XI(5)*(GXYZ(I,2)+GXYZ(I,4))))**2 &
                       +(CDABS(XI(4)*(GXYZ(I,1)+GXYZ(I,4))+XI(5)*(GXYZ(I,3)+GXYZ(I,2))))**2 &
                       +(CDABS(XI(4)*(GXYZ(I,3)+GXYZ(I,2))+XI(5)*(GXYZ(I,1)+GXYZ(I,4))))**2 &
                       +(CDABS(XI(4)*(GXYZ(I,4)+GXYZ(I,2))+XI(5)*(GXYZ(I,3)+GXYZ(I,1))))**2 &
                       +(CDABS(XI(4)*(GXYZ(I,3)+GXYZ(I,4))+XI(5)*(GXYZ(I,1)+GXYZ(I,2))))**2      
                 S(I) = AREA(E)*(SUM(WT(1:3)*F(1:3)))      
              ENDDO   
              errnrmx(4*e-eleshft) = dsqrt(s(1))           
              errnrmy(4*e-eleshft) = dsqrt(s(2))
              errnrmz(4*e-eleshft) = dsqrt(s(3))
           enddo 
           eleshft = eleshft-1
        enddo ! loop over kk
!
! Sum y and z terms to get total error
!
        do e=1,nele
           jj=irg(e)
           errnrm_eh(e) =  errnrmx(4*e-3)+ errnrmy(4*e-3)+ errnrmz(4*e-3)+ &
                           errnrmx(4*e-2)+ errnrmy(4*e-2)+ errnrmz(4*e-2)+ &
                           errnrmx(4*e-1)+ errnrmy(4*e-1)+ errnrmz(4*e-1)
!           errnrm_eh(e) = Max(errnrmx(4*e-3),errnrmy(4*e-3),errnrmz(4*e-3),errnrmx(4*e-2),errnrmy(4*e-2),errnrmz(4*e-2),errnrmx(4*e-1),errnrmy(4*e-1),errnrmz(4*e-1),errnrmx(4*e),errnrmy(4*e),errnrmz(4*e))
        enddo
    
        end subroutine l2norm_land       
!
!======================================================================! 
!============================================================== L2NORM !
!======================================================================! 
!
        subroutine l2norm_dual_land(dualgx,dualgy,dualgz,dualsmgradx,dualsmgrady,dualsmgradz, &
                    &   errnrm_eh,area)
!======================================================================! 
!                                                                      !
!  Written by: YUGUO LI, SIO, April, 2009                              !                     
!  Modified by: Yixin Ye                                               !
!               East china University of Technology                     !
!                                                                      !                       
!======================================================================!                   
        use   fem_params
        use   mesh
        use   datinfo             
        implicit none
    
        integer                              :: i,kk,eleshft,jj
        complex(double), dimension(4*nele)   :: dualgx ,dualgy , dualgz ! element gradient
        complex(double), dimension(4,4*nele) :: dualsmgradx, dualsmgrady, dualsmgradz  
        complex(double)                      :: gxyz(3,4),gx(4),gy(4),gz(4),eta
        real(double)                         :: wt(3),xi(5)
        real(double)                         :: s(3),f(3)
        real(double),dimension(4*nele)       :: errnrmx,errnrmy,errnrmz     
        complex(8)                           :: ommu
        real(double),    dimension(nele)     :: area,errnrm_eh  
        real(double)  ::omega
        common /kxw/omega
!
!  13 - Point Quadrature weights and points in area coordinates:
!
        DATA WT(1)/ 0.0131555555555D0/, WT(2)/0.0076222222222D0/,&
          WT(3)/ 0.0248888888888D0/

        DATA XI(1)/ 0.2500000000000D0/, XI(2)/0.785714285714286D0/,&
          XI(3)/ 0.0714285714285714D0/, XI(4)/0.399403576166799D0/,&
          XI(5)/ 0.100596423833201D0/

        ommu=ic*omega*mu0 
        eleshft = 3
        eta=-1.d0/ommu

        do kk = 1,4
           do e = 1,nele  
              s=0.d0
              do i=1,4
!                 gx(i)  = dualsmgradx(i,4*e-eleshft)-dualgx(4*e-eleshft)
!                 gy(i)  = dualsmgrady(i,4*e-eleshft)-dualgy(4*e-eleshft)
!                 gz(i)  = dualsmgradz(i,4*e-eleshft)-dualgz(4*e-eleshft)

                 GXYZ(1,I)=dualsmgradx(i,4*e-eleshft)-dualgx(4*e-eleshft)
                 GXYZ(2,I)=dualsmgrady(i,4*e-eleshft)-dualgy(4*e-eleshft)
                 GXYZ(3,I)=dualsmgradz(i,4*e-eleshft)-dualgz(4*e-eleshft)
              enddo
        
              DO I=1,3
                 F(1) = (CDABS(XI(1)*(GXYZ(I,1)+GXYZ(I,2)+GXYZ(I,3)+GXYZ(I,4))))**2

                 F(2) = (CDABS(XI(2)*GXYZ(I,1)+XI(3)*(GXYZ(I,2)+GXYZ(I,3)+GXYZ(I,4))))**2 &
                   +(CDABS(XI(2)*GXYZ(I,2)+XI(3)*(GXYZ(I,1)+GXYZ(I,3)+GXYZ(I,4))))**2 &
                   +(CDABS(XI(2)*GXYZ(I,3)+XI(3)*(GXYZ(I,1)+GXYZ(I,2)+GXYZ(I,4))))**2 &
                   +(CDABS(XI(2)*GXYZ(I,4)+XI(3)*(GXYZ(I,1)+GXYZ(I,2)+GXYZ(I,3))))**2 

                 F(3) = (CDABS(XI(4)*(GXYZ(I,1)+GXYZ(I,2))+XI(5)*(GXYZ(I,3)+GXYZ(I,4))))**2 &
                   +(CDABS(XI(4)*(GXYZ(I,1)+GXYZ(I,3))+XI(5)*(GXYZ(I,2)+GXYZ(I,4))))**2 &
                   +(CDABS(XI(4)*(GXYZ(I,1)+GXYZ(I,4))+XI(5)*(GXYZ(I,3)+GXYZ(I,2))))**2 &
                   +(CDABS(XI(4)*(GXYZ(I,3)+GXYZ(I,2))+XI(5)*(GXYZ(I,1)+GXYZ(I,4))))**2 &
                   +(CDABS(XI(4)*(GXYZ(I,4)+GXYZ(I,2))+XI(5)*(GXYZ(I,3)+GXYZ(I,1))))**2 &
                   +(CDABS(XI(4)*(GXYZ(I,3)+GXYZ(I,4))+XI(5)*(GXYZ(I,1)+GXYZ(I,2))))**2
    
                 S(I) = area(e)*(SUM(WT(1:3)*F(1:3)))  
              ENDDO 
              errnrmx(4*e-eleshft) = dsqrt(s(1))           
              errnrmy(4*e-eleshft) = dsqrt(s(2))
              errnrmz(4*e-eleshft) = dsqrt(s(3))    
           enddo 
           eleshft = eleshft-1
        enddo
!              
        do e=1,nele
           jj=irg(e)
           errnrm_eh(e) = errnrm_eh(e)*(errnrmx(4*e-3)+ errnrmy(4*e-3)+ errnrmz(4*e-3)+&
                                      errnrmx(4*e-2)+ errnrmy(4*e-2)+ errnrmz(4*e-2)+&
                                      errnrmx(4*e-1)+ errnrmy(4*e-1)+ errnrmz(4*e-1))
        enddo 

        end subroutine l2norm_dual_land
!      
!======================================================================! 
!========================================================= SMQHREGIONS !
!======================================================================! 
!    
        subroutine SmQhregions_land(egrad,nregion,eregion,esmgrad,area)  
    
!======================================================================!
!  Computes the smoothed and projected gradient in each region         !
!  using the Bank and Xu (2003b) approach using SmQh operators.        !
! 								       !		
!  Modified by: Yixin Ye                                               !
!               East china University of Technology                     ! 
!                                                                      !                       
!======================================================================!     
        use fem_params
        use mesh
        use datinfo
        implicit none
    
        integer                            :: ireg,j,n,kk,eleshft
        integer                            :: nregion,eregion(nele)
        real(double),    dimension(nele)   :: area
        complex(double), dimension(4*nele)   :: egrad 
        complex(double), dimension(4,4*nele) :: esmgrad 
        complex(double), allocatable, dimension(:)  :: Qhrhs, diag, Sm
    
        external    inprod_Qh_k2,inprod_sm_k1
     
        allocate ( Qhrhs(nnod),diag(nnod),Sm(nnod) )
!
!  Loop over regions and compute L2 projection of gradient, then Smooth 
!
! loop for exs then hxs:
        eleshft = 3

        do kk=1,4
           do ireg=1,nregion   
              Qhrhs     = (0d0,0d0)      
!
! Loop over elements, if element in region, loop over 4 nodes and adjust 
! boundary condition vector to include node in system
!
              bcnod = 2  ! default is to ignore all nodes from solution 
              do e=1,nele    
                 if (eregion(e).eq.ireg) then   ! element is in region ireg
                    do j = 1,4        
                       n = emap(j,e)
                       bcnod(n) = 0  ! node is included in computation
!  
! Integral of gradient value with basis functions over element
!
                       Qhrhs(n) = Qhrhs(n) + area(e)/4.d0*egrad(4*e-eleshft)         
                    enddo
                 endif ! eregrion(e).eq.i
              enddo   ! loop over elements
!
! Generate diagonal of Matrix as a preconditioner for L2 Projection system
!
              call gen_diagQh_k2(diag)
!
! Solve linear system for L2 projection using QMR
!
              Sm = (0d0,0d0) ! dummy array for guess
              call qmr(Qhrhs,SM,diag,nnod,inprod_Qh_k2,qhmaxiter,qhmaxerr)
!
! Do smoothing operations by applying Laplacian to QhGradU 
! using nsmiter iterations of QMR
!
              if (nsmiter.gt.0) then
      	         diag = (1.d0,0.d0)  ! no preconditioner
      	         Sm = (0d0,0d0)
      	         call qmr(Sm,Qhrhs,diag,nnod,inprod_Sm_k1,nsmiter,1d-25) 
              else
      	         Sm = Qhrhs ! copy over L2 projection array
              endif
! 
! Extract used nodes and put into esmgradx,z arrays:
! 
              do e=1,nele    
                 if (eregion(e).eq.ireg) then   ! element is in region ireg
                    do j = 1,4        
                       n = emap(j,e)
                       esmgrad(j,4*e-eleshft) = Sm(n)
                    enddo
                 endif ! eregrion(e).eq.i
              enddo   ! loop over elements
!
           enddo  ! loop over regions  
           eleshft = eleshft-1
        enddo !loop over kk
!
!  Deallocate arrays
!
        deallocate (Qhrhs,diag,Sm)
!
        end subroutine SmQhregions_land
!    
!======================================================================! 
!============================================================= GETGRAD !
!======================================================================! 
!
        subroutine getgradyz(rhs0,egradx,egrady,egradz,area)
!
!======================================================================!
!  Subroutine that computes the piecewise constant gradient of rhs1    !
!  within each element.						       !	
!   								       !                        
!======================================================================! 
    
        use fem_params
        use mesh   
        implicit none
    
        integer                              :: i,n(4)
        real(double),    dimension(4)        :: dndx,dndy,dndz,xe,ye,ze,a,b,c
        real(double) ,    dimension(nele)    :: area
        complex(double), dimension(4*nele)   :: egradx,egrady,egradz 
        complex(double), dimension(4*nnod)   :: rhs0    
!
!  Compute element gradient using linear basis
! 
        do e=1,nele   
           do i=1,4
              n(i)   = emap(i,e)
              xe(i)  = xn(n(i))
              ye(i)  = yn(n(i))
              ze(i)  = zn(n(i))
           enddo 
!      
!  Compute the coefficients of linear interpolation in an element.
!
           a(1) = ye(4)*(ze(3)-ze(2)) + ye(2)*(ze(4)-ze(3)) + ye(3)*(ze(2)-ze(4))
           a(2) = ye(4)*(ze(1)-ze(3)) + ye(1)*(ze(3)-ze(4)) + ye(3)*(ze(4)-ze(1))
           a(3) = ye(4)*(ze(2)-ze(1)) + ye(1)*(ze(4)-ze(2)) + ye(2)*(ze(1)-ze(4))
           a(4) = ye(3)*(ze(1)-ze(2)) + ye(1)*(ze(2)-ze(3)) + ye(2)*(ze(3)-ze(1))
           b(1) = xe(3)*(ze(4)-ze(2)) + xe(4)*(ze(2)-ze(3)) + xe(2)*(ze(3)-ze(4))
           b(2) = xe(3)*(ze(1)-ze(4)) + xe(4)*(ze(3)-ze(1)) + xe(1)*(ze(4)-ze(3))
           b(3) = xe(2)*(ze(4)-ze(1)) + xe(4)*(ze(1)-ze(2)) + xe(1)*(ze(2)-ze(4))
           b(4) = xe(2)*(ze(1)-ze(3)) + xe(3)*(ze(2)-ze(1)) + xe(1)*(ze(3)-ze(2))
           c(1) = xe(4)*(ye(3)-ye(2)) + xe(2)*(ye(4)-ye(3)) + xe(3)*(ye(2)-ye(4))
           c(2) = xe(4)*(ye(1)-ye(3)) + xe(1)*(ye(3)-ye(4)) + xe(3)*(ye(4)-ye(1))
           c(3) = xe(4)*(ye(2)-ye(1)) + xe(1)*(ye(4)-ye(2)) + xe(2)*(ye(1)-ye(4))
           c(4) = xe(3)*(ye(1)-ye(2)) + xe(1)*(ye(2)-ye(3)) + xe(2)*(ye(3)-ye(1))  
!
!  Compute the volume of the current element
!
           area(e)= ((xe(2)-xe(1))*(ye(3)*ze(4)-ze(3)*ye(4)) + (xe(1)- xe(3))*(ye(2)*ze(4)-ze(2)*ye(4)) + (xe(4)-xe(1))*(ye(2)*ze(3)-ze(2)*ye(3)) &
             &  +(xe(3)- xe(2))*(ye(1)*ze(4)-ze(1)*ye(4)) + (xe(2)- xe(4))*(ye(1)*ze(3)-ze(1)*ye(3)) + (xe(4)- xe(3))*(ye(1)*ze(2)-ze(1)*ye(2)))/6.d0       
!  
!  Compute the shape functions and their derivatives with respecting to
!  x- y- and z-directions
!      
           do i=1,4
              dndx(i) = a(i)/area(e)/6.d0
              dndy(i) = b(i)/area(e)/6.d0
              dndz(i) = c(i)/area(e)/6.d0
           enddo
! for Axs: 
           egradx(4*e-3) = sum(dndx(1:4)*rhs0(4*n(1:4)-3))  
           egrady(4*e-3) = sum(dndy(1:4)*rhs0(4*n(1:4)-3)) 
           egradz(4*e-3) = sum(dndz(1:4)*rhs0(4*n(1:4)-3)) 
! for Ays:
           egradx(4*e-2) = sum(dndx(1:4)*rhs0(4*n(1:4)-2))  
           egrady(4*e-2) = sum(dndy(1:4)*rhs0(4*n(1:4)-2)) 
           egradz(4*e-2) = sum(dndz(1:4)*rhs0(4*n(1:4)-2))
! for Azs:
           egradx(4*e-1) = sum(dndx(1:4)*rhs0(4*n(1:4)-1))  
           egrady(4*e-1) = sum(dndy(1:4)*rhs0(4*n(1:4)-1)) 
           egradz(4*e-1) = sum(dndz(1:4)*rhs0(4*n(1:4)-1))
! for phis:
           egradx(4*e  ) = sum(dndx(1:4)*rhs0(4*n(1:4)  ))  
           egrady(4*e  ) = sum(dndy(1:4)*rhs0(4*n(1:4)  )) 
           egradz(4*e  ) = sum(dndz(1:4)*rhs0(4*n(1:4)  ))
        enddo      
!    
        end subroutine getgradyz
!
!======================================================================! 
!========================================================== GETREGIONS !
!======================================================================!    
!    
        subroutine getregions(nregion,eregion)
    
!======================================================================!   
! Define regions (possibly unconnected) of unique conductivity.        !
! Discontinuities in grad(Ex,Hx) need to be allowed so L2 projection   !
! Qh and smoothing Sm need to be done over subregions of piecewise     !
! continuous sigma.            !
!                                                                      !
! Written by: Kerry Key                                                !
!              Scripps Institution of Oceanography                     !
!              2005-2006                                               !
! Modified By Yuguo Li, Summer 2008
!                                                                      !                       
!======================================================================!    
        use fem_params
        use mesh
!       use datinfo    
        implicit none
     
        integer                                :: i, nregion
        integer, dimension(nele)               :: eregion
        real(double),allocatable, dimension(:) :: irgregion
    
        allocate( irgregion(nele) )

        nregion = 0
        eregion = 0  ! set default to 0
        irgregion=0
    
        do e=1,nele
           i = 0
           do while (i.lt.nregion)  ! loop through existing regions
              i = i+1
              if (irg(e).eq.irgregion(i)) then
                 eregion(e)  = i  ! assign region number to element
                 i = nregion + 1  ! exit loop if region found
              endif
           enddo

           if (eregion(e).eq.0) then 
              nregion            = nregion+1
              eregion(e)         = nregion
              irgregion(nregion) = irg(e)
           endif
        enddo  ! loop over elements

        deallocate (irgregion)

        end subroutine getregions
!    
!======================================================================! 
!======================================================= GEN_DIAGQH_K2 !
!======================================================================! 
!
        subroutine gen_diagQh_k2(diag)

!======================================================================! 
!                                                                      !
!  Routine to generate the diagonal of the glogal FEM matrix for       !
!   L2 projection operator Qh                                          !
!                                                                      !
!  Written by: Yuguo Li, DEC 10, 2008                                  !
!              Scripps Institution of Oceanography                     !
!                                                                      ! 
!======================================================================! 
        use             fem_params
        use             mesh      
        implicit none
      
        real(double),   dimension(10)     :: mmtx
        integer,        dimension(4)      :: n
        integer                           :: i
        real(double),   dimension(4)      :: xe,ye,ze
        complex(double),dimension(nnod)   :: diag
!
!  Initialize the value on the diagonal of the FEM matrix
!
        diag = (0.d0,0.d0)
!
!  Loop over the elements in the FE mesh
!
        do e=1,nele
           iskip = 0
           do i=1,4
              n(i) = emap(i,e)
              xe(i) = xn(n(i))
              ye(i) = yn(n(i))
              ze(i) = zn(n(i))
              if (bcnod(n(i)).ge.1) iskip=1
           enddo

           if (iskip.eq.0) then
              call  stiff_k2(xe,ye,ze,mmtx) 
!
!  Extract the diagonals from the element matrices and append
!  the current values on the global matrix diagonal with 
!  their values.
!
              diag(n(1)) = diag(n(1)) + mmtx(1)
              diag(n(2)) = diag(n(2)) + mmtx(5)
              diag(n(3)) = diag(n(3)) + mmtx(8)
              diag(n(4)) = diag(n(4)) + mmtx(10)
           else
              do i=1,4
                 if (bcnod(n(i)).eq.2) diag(n(i))=(1.d0,0.d0)
              enddo
           endif             
        enddo 
!
        end subroutine gen_diagQh_k2     
!            
!======================================================================! 
!=============================================================STIFF_k1 !
!======================================================================! 
!
        subroutine stiff_k1(xe,ye,ze,mtx)

!======================================================================!
!                                                                      ! 
!  Compute the 3x3 element matrix K_{1e}                               !
!                                                                      ! 
!  Written by: Yuguo Li,  Dec 12 2006.                                 !
!                                                                      ! 
!======================================================================!    
        use             fem_params
        implicit none

        real(double)                :: del,g1
        real(double), dimension(10) :: mtx
        real(double), dimension(4)  :: xe,ye,ze,a,b,c
        integer                     :: i,j,k
      
        a(1) = ye(4)*(ze(3)-ze(2)) + ye(2)*(ze(4)-ze(3)) + ye(3)*(ze(2)-ze(4))
        a(2) = ye(4)*(ze(1)-ze(3)) + ye(1)*(ze(3)-ze(4)) + ye(3)*(ze(4)-ze(1))
        a(3) = ye(4)*(ze(2)-ze(1)) + ye(1)*(ze(4)-ze(2)) + ye(2)*(ze(1)-ze(4))
        a(4) = ye(3)*(ze(1)-ze(2)) + ye(1)*(ze(2)-ze(3)) + ye(2)*(ze(3)-ze(1))
        b(1) = xe(3)*(ze(4)-ze(2)) + xe(4)*(ze(2)-ze(3)) + xe(2)*(ze(3)-ze(4))
        b(2) = xe(3)*(ze(1)-ze(4)) + xe(4)*(ze(3)-ze(1)) + xe(1)*(ze(4)-ze(3))
        b(3) = xe(2)*(ze(4)-ze(1)) + xe(4)*(ze(1)-ze(2)) + xe(1)*(ze(2)-ze(4))
        b(4) = xe(2)*(ze(1)-ze(3)) + xe(3)*(ze(2)-ze(1)) + xe(1)*(ze(3)-ze(2))
        c(1) = xe(4)*(ye(3)-ye(2)) + xe(2)*(ye(4)-ye(3)) + xe(3)*(ye(2)-ye(4))
        c(2) = xe(4)*(ye(1)-ye(3)) + xe(1)*(ye(3)-ye(4)) + xe(3)*(ye(4)-ye(1))
        c(3) = xe(4)*(ye(2)-ye(1)) + xe(1)*(ye(4)-ye(2)) + xe(2)*(ye(1)-ye(4))
        c(4) = xe(3)*(ye(1)-ye(2)) + xe(1)*(ye(2)-ye(3)) + xe(2)*(ye(3)-ye(1))
!
!  Compute the area of the current element,multiplied by six
!
        del = dabs((xe(2)-xe(1))*(ye(3)*ze(4)-ze(3)*ye(4)) + (xe(1)- xe(3))*(ye(2)*ze(4)-ze(2)*ye(4)) + (xe(4)-xe(1))*(ye(2)*ze(3)-ze(2)*ye(3)) &
             &  +(xe(3)- xe(2))*(ye(1)*ze(4)-ze(1)*ye(4)) + (xe(2)- xe(4))*(ye(1)*ze(3)-ze(1)*ye(3)) + (xe(4)- xe(3))*(ye(1)*ze(2)-ze(1)*ye(2)))
      
        g1=1.d0/(6.d0*del)
!
!  Compute the elements of the element  matrix.
!
        mtx(1) = g1*(a(1)*a(1)+b(1)*b(1)+c(1)*c(1))  
        mtx(2) = g1*(a(1)*a(2)+b(1)*b(2)+c(1)*c(2)) 
        mtx(3) = g1*(a(1)*a(3)+b(1)*b(3)+c(1)*c(3)) 
        mtx(4) = g1*(a(1)*a(4)+b(1)*b(4)+c(1)*c(4)) 
        mtx(5) = g1*(a(2)*a(2)+b(2)*b(2)+c(2)*c(2))   
        mtx(6) = g1*(a(2)*a(3)+b(2)*b(3)+c(2)*c(3))
        mtx(7) = g1*(a(2)*a(4)+b(2)*b(4)+c(2)*c(4))
        mtx(8) = g1*(a(3)*a(3)+b(3)*b(3)+c(3)*c(3))
        mtx(9) = g1*(a(3)*a(4)+b(3)*b(4)+c(3)*c(4))
        mtx(10) = g1*(a(4)*a(4)+b(4)*b(4)+c(4)*c(4))

        end subroutine stiff_k1
!
!======================================================================! 
!============================================================ Stiff_k2 !
!======================================================================!
! 
        subroutine stiff_k2(xe,ye,ze,mtx)
!
!======================================================================!
!                                                                      ! 
!  Compute the 3x3 element matrix K_{2e}                               !
!                                                                      ! 
!  Written by: Yuguo Li, DEC 12, 2006.                                 !
!                                                                      ! 
!======================================================================!    
        use             fem_params
        implicit none
      
        real(double), dimension(10)  :: mtx
        real(double), dimension(4)  :: xe,ye,ze
        real(double), dimension(4)  :: a,b,c
        real(double)                :: del
!
!  Compute the area of the current element,multiplied by six
!
        del = dabs((xe(2)-xe(1))*(ye(3)*ze(4)-ze(3)*ye(4)) + (xe(1)- xe(3))*(ye(2)*ze(4)-ze(2)*ye(4)) + (xe(4)-xe(1))*(ye(2)*ze(3)-ze(2)*ye(3)) &
             &  +(xe(3)- xe(2))*(ye(1)*ze(4)-ze(1)*ye(4)) + (xe(2)- xe(4))*(ye(1)*ze(3)-ze(1)*ye(3)) + (xe(4)- xe(3))*(ye(1)*ze(2)-ze(1)*ye(2)))
      

        mtx(1) = del/60.d0
        mtx(2) = del/120.d0
        mtx(3) = mtx(2)
        mtx(4) = mtx(2)
        mtx(5) = mtx(1)
        mtx(6) = mtx(2)
        mtx(7) = mtx(2)
        mtx(8) = mtx(1)
        mtx(9) = mtx(2)
        mtx(10) = mtx(1)

        end subroutine stiff_k2
!
!======================================================================! 
!======================================================== INPROD_QH_K2 !
!======================================================================! 
!
        subroutine inprod_Qh_k2(ndf,invec,outvec,matflops)
!
!======================================================================!
! Computes inner product for L2 projection                             !
! Product of mass matrix and vector                                    !	                                    ! 
!  Written by: Yuguo Li                                                !                                                                     ! 
!======================================================================!
        use           fem_params
        use           mesh      
        implicit none
      
        integer                          :: ndf,matflops,i
        integer, dimension(4)            :: n
        real(double),    dimension(10)   :: mtx
        real(double),    dimension(4)    :: xe,ye,ze
        complex(double), dimension(ndf)  :: invec,outvec
        real(double)                     :: del
!
!  Set the number of FLOP's for each call to INPROD_MT
!
        matflops = 87*nele
!
!  Initialize the value of the output vector
!
        outvec = (0.d0,0.d0)
!
!  Loop over elements in the FE mesh
!
        do e=1,nele
           iskip = 0
           do i=1,4
              n(i) = emap(i,e)
              xe(i) = xn(n(i))
              ye(i) = yn(n(i))
              ze(i) = zn(n(i))
              if (bcnod(n(i)).eq.2) iskip=1
           enddo
            
           if (iskip.eq.0) then
!
!  Compute the element mass matrix.  Note that 
!  the indices of mmtx correspond to matrix elements
!  (1,1),(1,2),(1,3),(2,2),(2,3) and (3,3) where the mtx is
!  symmetric.
     
              call  stiff_k2(xe,ye,ze,mtx) 
!  
!  Contract the complete element matrix with the input vector and 
!  then add the result to the output vector
!

 
              outvec(n(1)) = outvec(n(1))                                     &
                               + (mtx(1)) * invec(n(1)) &
                               + (mtx(2)) * invec(n(2)) &
                               + (mtx(3)) * invec(n(3)) &
                               + (mtx(4)) * invec(n(4)) 
                
              outvec(n(2)) = outvec(n(2))                                     &
                               + (mtx(2)) * invec(n(1)) &
                               + (mtx(5)) * invec(n(2)) &
                               + (mtx(6)) * invec(n(3)) &
                               + (mtx(7)) * invec(n(4))
                
              outvec(n(3)) = outvec(n(3))                                     &
                               + (mtx(3)) * invec(n(1)) &
                               + (mtx(6)) * invec(n(2)) &
                               + (mtx(8)) * invec(n(3)) &
                               + (mtx(9)) * invec(n(4))

              outvec(n(4)) = outvec(n(4))                                     &
                               + (mtx(4)) * invec(n(1)) &
                               + (mtx(7)) * invec(n(2)) &
                               + (mtx(9)) * invec(n(3)) &
                               + (mtx(10)) * invec(n(4))
           else
              do i=1,4
                 if (bcnod(n(i)).eq.2) outvec(n(i))=invec(n(i))
              enddo
           endif

        enddo
   
        end subroutine inprod_Qh_k2
!     
!======================================================================! 
!=========================================================== INPROD_K1 !
!======================================================================! 
!
        subroutine inprod_k1(ndf,invec,outvec,matflops)
!
!======================================================================!
! Computes inner product for L2 projection                             !
! Product of mass matrix and vector                                    !	                                    ! 
!  Written by: Yuguo Li                                                !                                                                     ! 
!======================================================================!
        use           fem_params
        use           mesh 
        use           datinfo     
        implicit none
      
        complex(8) :: k0,ommu
        complex(8) :: ke1(16,16)    
        integer                          :: ndf,matflops,i,jj
        integer, dimension(4)            :: n
        real(double),    dimension(10)   :: mtx
        real(double),    dimension(4)    :: xe,ye,ze
        complex(double), dimension(ndf)  :: invec,outvec
        real(double)                     :: del
        real(double)  ::omega,sige(6)
        common /kxw/omega
!
!  Set the number of FLOP's for each call to INPROD_MT
!
        matflops = 87*nele

        ommu=ic*omega*mu0
!
!  Initialize the value of the output vector
!
        outvec = (0.d0,0.d0)
!
!  Loop over elements in the FE mesh
!
        do e=1,nele
!          iskip = 0
           jj=irg(e)
           do i=1,4
              n(i) = emap(i,e)
              xe(i) = xn(n(i))
              ye(i) = yn(n(i))
              ze(i) = zn(n(i))
!             if (bcnod(n(i)).eq.2) iskip=1
           enddo
            
!          if (iskip.eq.0) then
!
!  Compute the element mass matrix.  Note that 
!  the indices of mmtx correspond to matrix elements
!  (1,1),(1,2),(1,3),(2,2),(2,3) and (3,3) where the mtx is
!  symmetric.
!
           sige=sg(jj,:)

           call gen_K1(xe,ye,ze,del,ommu,sige,Ke1) 
     
!          call  stiff_k2(xe,ye,ze,mtx) 
!  
!  Contract the complete element matrix with the input vector and 
!  then add the result to the output vector
!
           do i=1,4
              if (bcnod(n(i)).eq.0) then
                 outvec(4*n(i)-3) = outvec(4*n(i)-3)                                     &
                               + (ke1(i,1)) * invec(4*n(1)-3) &
                               + (ke1(i,5)) * invec(4*n(1)-2) &
                               + (ke1(i,9)) * invec(4*n(1)-1) &
                               + (ke1(i,13)) * invec(4*n(1)  ) &
                               + (ke1(i,2)) * invec(4*n(2)-3) &
                               + (ke1(i,6)) * invec(4*n(2)-2) &
                               + (ke1(i,10)) * invec(4*n(2)-1) &
                               + (ke1(i,14)) * invec(4*n(2)  ) &
                               + (ke1(i,3)) * invec(4*n(3)-3) &
                               + (ke1(i,7)) * invec(4*n(3)-2) &
                               + (ke1(i,11)) * invec(4*n(3)-1) &
                               + (ke1(i,15)) * invec(4*n(3)  ) &
                               + (ke1(i,4)) * invec(4*n(4)-3) &
                               + (ke1(i,8)) * invec(4*n(4)-2) &
                               + (ke1(i,12)) * invec(4*n(4)-1) &
                               + (ke1(i,16)) * invec(4*n(4)  ) 


                 outvec(4*n(i)-2) = outvec(4*n(i)-2)                                     &
                               + (ke1(i+4,1)) * invec(4*n(1)-3) &
                               + (ke1(i+4,5)) * invec(4*n(1)-2) &
                               + (ke1(i+4,9)) * invec(4*n(1)-1) &
                               + (ke1(i+4,13)) * invec(4*n(1)  ) &
                               + (ke1(i+4,2)) * invec(4*n(2)-3) &
                               + (ke1(i+4,6)) * invec(4*n(2)-2) &
                               + (ke1(i+4,10)) * invec(4*n(2)-1) &
                               + (ke1(i+4,14)) * invec(4*n(2)  ) &
                               + (ke1(i+4,3)) * invec(4*n(3)-3) &
                               + (ke1(i+4,7)) * invec(4*n(3)-2) &
                               + (ke1(i+4,11)) * invec(4*n(3)-1) &
                               + (ke1(i+4,15)) * invec(4*n(3)  ) &
                               + (ke1(i+4,4)) * invec(4*n(4)-3) &
                               + (ke1(i+4,8)) * invec(4*n(4)-2) &
                               + (ke1(i+4,12)) * invec(4*n(4)-1) &
                               + (ke1(i+4,16)) * invec(4*n(4)  )
                
                 outvec(4*n(i)-1) = outvec(4*n(i)-1)                                     &
                               + (ke1(i+8,1)) * invec(4*n(1)-3) &
                               + (ke1(i+8,5)) * invec(4*n(1)-2) &
                               + (ke1(i+8,9)) * invec(4*n(1)-1) &
                               + (ke1(i+8,13)) * invec(4*n(1)  ) &
                               + (ke1(i+8,2)) * invec(4*n(2)-3) &
                               + (ke1(i+8,6)) * invec(4*n(2)-2) &
                               + (ke1(i+8,10)) * invec(4*n(2)-1) &
                               + (ke1(i+8,14)) * invec(4*n(2)  ) &
                               + (ke1(i+8,3)) * invec(4*n(3)-3) &
                               + (ke1(i+8,7)) * invec(4*n(3)-2) &
                               + (ke1(i+8,11)) * invec(4*n(3)-1) &
                               + (ke1(i+8,15)) * invec(4*n(3)  ) &
                               + (ke1(i+8,4)) * invec(4*n(4)-3) &
                               + (ke1(i+8,8)) * invec(4*n(4)-2) &
                               + (ke1(i+8,12)) * invec(4*n(4)-1) &
                               + (ke1(i+8,16)) * invec(4*n(4)  )

                 outvec(4*n(i)) = outvec(4*n(i))                                     &
                               + (ke1(i+12,1)) * invec(4*n(1)-3) &
                               + (ke1(i+12,5)) * invec(4*n(1)-2) &
                               + (ke1(i+12,9)) * invec(4*n(1)-1) &
                               + (ke1(i+12,13)) * invec(4*n(1)  ) &
                               + (ke1(i+12,2)) * invec(4*n(2)-3) &
                               + (ke1(i+12,6)) * invec(4*n(2)-2) &
                               + (ke1(i+12,10)) * invec(4*n(2)-1) &
                               + (ke1(i+12,14)) * invec(4*n(2)  ) &
                               + (ke1(i+12,3)) * invec(4*n(3)-3) &
                               + (ke1(i+12,7)) * invec(4*n(3)-2) &
                               + (ke1(i+12,11)) * invec(4*n(3)-1) &
                               + (ke1(i+12,15)) * invec(4*n(3)  ) &
                               + (ke1(i+12,4)) * invec(4*n(4)-3) &
                               + (ke1(i+12,8)) * invec(4*n(4)-2) &
                               + (ke1(i+12,12)) * invec(4*n(4)-1) &
                               + (ke1(i+12,16)) * invec(4*n(4)  )

              else
                 outvec(4*n(i)-3)=invec(4*n(i)-3)
                 outvec(4*n(i)-2)=invec(4*n(i)-2)
                 outvec(4*n(i)-1)=invec(4*n(i)-1)
                 outvec(4*n(i)  )=invec(4*n(i)  )
              endif 
           enddo

        enddo
   
        end subroutine inprod_k1            
!            
!======================================================================! 
!======================================================== INPROD_SM_K1 !
!======================================================================! 
!
        subroutine inprod_Sm_k1(ndf,invec,outvec,matflops)
!
!======================================================================!
! Computes inner product for smoothness operator		       !
! Product of stiffness matrix and vector                               !
!  Written by: Yuguo Li, Dec 10 2008                                   !
!								       ! 
!======================================================================!
        use           fem_params
        use           mesh
        implicit none
      
        integer                          :: ndf,matflops,i
        integer, dimension(4)            :: n
        real(double),    dimension(10)   :: mtx
        complex(double), dimension(ndf)  :: invec,outvec
        real(double), dimension(4)       :: xe,ye,ze 
!
!  Set the number of FLOP's for each call to INPROD_MT
!
        matflops = 87*nele
!
!  Initialize the value of the output vector
!
        outvec = (0.d0,0.d0)
!
!  Loop over elements in the FE mesh
!
        do e=1,nele
           iskip = 0
           do i=1,4
              n(i)  = emap(i,e)
              xe(i) = xn(n(i))
              ye(i) = yn(n(i))
              ze(i) = zn(n(i))
              if (bcnod(n(i)).eq.2) iskip=1
           enddo
            
           if (iskip.eq.0) then
!
!  Compute the element matrix k1e. 
!
              call stiff_k1(xe,ye,ze,mtx) 

!
!  Contract the complete element matrix with the input vector and 
!  then add the result to the output vector
! 
              outvec(n(1)) = outvec(n(1))                                     &
                               + (mtx(1)) * invec(n(1)) &
                               + (mtx(2)) * invec(n(2)) &
                               + (mtx(3)) * invec(n(3)) &
                               + (mtx(4)) * invec(n(4))
                
              outvec(n(2)) = outvec(n(2))                                     &
                               + (mtx(2)) * invec(n(1)) &
                               + (mtx(5)) * invec(n(2)) &
                               + (mtx(6)) * invec(n(3)) &
                               + (mtx(7)) * invec(n(4))
                
              outvec(n(3)) = outvec(n(3))                                     &
                               + (mtx(3)) * invec(n(1)) &
                               + (mtx(6)) * invec(n(2)) &
                               + (mtx(8)) * invec(n(3)) &
                               + (mtx(9)) * invec(n(4))

              outvec(n(4)) = outvec(n(4))                                     &
                               + (mtx(4)) * invec(n(1)) &
                               + (mtx(7)) * invec(n(2)) &
                               + (mtx(9)) * invec(n(3)) &
                               + (mtx(10)) * invec(n(4))

           else
              do i=1,4
                 if (bcnod(n(i)).eq.2) outvec(n(i))=invec(n(i))
              enddo
           endif 

        enddo

        end subroutine inprod_Sm_k1
!
!======================================================================!
!================================================================= QMR !
!======================================================================!
!
        subroutine qmr(rhs,guess,diag,ndf,inprod,MAXIT,MAXERR)
!----------------------------------------------------------------------!
!                                                                      !
!  Routine to solve a complex symmetric system of linear equations     !
!  using the Quasi-Residual-Method (QMR).  The specific algorithm used !
!  here is the coupled two-term recurrence scheme without 'look-ahead' !
!  outlined in Freund  and Nachtigal, SIAM J Sci. Comput., v15,        !
!  313-337 (1994), algorithm 8.1.   Routine implements a diagonal      !
!  preconditioning matrix into the Lanczos iterates.  Multiplication   !
!  by the preconditioner of the right-hand-side vector upon entering   !
!  or multiplication of the solution upon exiting the routine are NOT  !
!  necessary.  Nor is it necessary to modify the routine 'inprod'      !
!  which computes the action of the FD matrix upon an arbitrary vector.!
!                                                                      !
!  Written by: C J Weiss, Sandia National Labs., April 25, 2000.       !
!                                                                      !
!                           Modifications                              !
!                           =============                              !
!                                                                      !
!  Oct 3 2000: Inverse of the diagonal computed and stored to reduce   !
!               cost of Jacobi scaling.                                !
!                                                                      !
!  Oct 3 2000: Reduced storage cost by eliminating a specific vector   !
!               to store the current solution and instead storing this !
!               solution in the vector containing the RHS since the    !
!               RHS is only required once, before the main iterative   !
!               loop.                                                  !
!                                                                      !
!  Sep 17 2002: Removed all the _x,_y,_z rubbish for use in 2DMT FEM   !
!                modeling.                                             !
!                                                                      !
!  2005-2006:   KWK.  Minor changes for adding a starting "guess" and  !
!               output a log file. Tolerance stopping                  !
!               Added functionality for performing smoothing           !
!               operations
!
!
!
!  KWK Notes:
!  ========
!
!  rhs:       On input rhs is b in eqn Ax =b.
!             On output rhs is x in eqn Ax = b.
!
!  guess:     A starting guess for x in Ax = b.  If doing a few
!             smoothing iterations then guess is the vector to
!             be smoothed by partially solving nabla^2 (guess) = 0
!
!----------------------------------------------------------------------!

        real(8)                       :: MAXERR
        integer                       :: MAXIT
        complex(8), dimension(ndf)    :: guess

        integer                       :: ndf
        complex(8), dimension(ndf)    :: rhs,diag

        complex(8), dimension(:)      :: inv_diag
        complex(8), dimension(:)      :: p,d,vt,v,tmp,r,s

        allocatable      :: inv_diag
        allocatable      :: p,d,vt,v,tmp,r,s

        real(4),    dimension(2)    :: time

        complex(8)  :: eps,del,beta,eta,k1
        real(8)     :: k2,rho_n,rho_np, rho_np0
        real(8)     :: c_nm,c_n,theta_n,theta_nm
        integer     :: iter
        real(8)     :: flops,matflops
        real(4)     :: t1,t2,t3,etime_
        real(4)     :: t11,t22,t(6)
        real(8)     :: rnorm,rhsnorm, guessnorm
!
! KWK 05/01/2006 tolerance: if no longer converging, then stop
!
        integer, parameter :: ntol = 50 ! 20 iterations
        integer            :: itol ! counter
        real(8), parameter :: tol = 1.d-6 ! relative tolerance for quitting
        real(8)            :: rlast

! Variables !----------------------------------------------------------!
!                                                                      !
!  rhs:            Upon entering, the right hand side for the linear   !
!                    system;  once inside the main loop, the current   !
!                    solution estimate.                                !
!  diag:           Diagonal of coefficient matrix.                     !
!  inv_diag:       Inverse of the diagonal of coefficient matrix, used !
!                    here as a Jacobi-style preconditioner.            !
!  ndf:            Number of degrees of freedom in the linear system.  !
!  d:              Direction vector for solution updates.              !
!  vt:             Unnormalized Lanczos vector.                        !
!  v,p:            The two 'Lanczos vectors'.                          !
!  tmp:            Action of the coeficient matrix on a given vector.  !
!  r:              Residual vector.                                    !
!  s:              Supplementary vector used to compute residual via   !
!                    vector updates instead of matrix-vector products  !
!                    at each interation.                               !
!  k1,k2:          Convenient constants for complex double and real    !
!                  double precision, respectively.                     !
!  eps,del,beta,eta,theta_nm,theta_n:  ???                             !
!  rho_n,rho_np:   L2 norms of the 'v' Lanczos vector for interations  !
!                    n and n+1.                                        !
!  c_nm,cn:        Elements of the Given's rotation matrix used in the !
!                    QR decomposition of the coefficient matrix for    !
!                    iterations n and n+1.                             !
!  iter:           Iteration counter.                                  !
!  dtime,time:     System call for elapsed time in QMR routine.        !
!  t1,t2:          Starting and ending times for QMR routine.          !
!  rnorm:          L2 norm of the residual vector.                     !
!  rhsnorm:        L2 norm of the RHS vector.                          !
!                                                                      !
!----------------------------------------------------------------------!
!
!  Allocate storage for the local arrays used in the QMR routine.
!
        allocate (inv_diag(ndf), p(ndf),   d(ndf) )
        allocate (      vt(ndf), v(ndf), tmp(ndf) )
        allocate (       r(ndf), s(ndf)           )
!
!  Initialize time counter...
!
!       t1 = etime_(time)
!
!  Open the output file.
!
        open(unit=11,file='qmr.iter',position='APPEND')
!
!  Compute the norm of the RHS vector
!
        rhsnorm = sum(rhs*dconjg(rhs))
        rhsnorm = dsqrt(rhsnorm)
!
! Compute norm of guess vector:
!
        guessnorm = sum(guess*dconjg(guess))
        guessnorm = dsqrt(guessnorm)
!
! Stop if both vectors are 0 since residual is then 0
!
        if ((rhsnorm.eq.0).and.(guessnorm.eq.0)) then
           write(11,*) ' QMR: trivial solution,returning to main routine.'
           rhs = (0d0,0d0)
           close(11)
           return
        endif
!
!  Compute the inverse of the diagonal
!
        inv_diag = 1.d0/diag
!
!  Politely introduce yourself.
!
!       write(11,*) '---------------------------------'
!       write(11,*) '  2-term coupled recurrence QMR  '
!       write(11,*) '    with Jacobi preconditioner   '
!       write(11,*) '================================='
        write(11,*) ' '
        write(11,*) ' Target residual norm: ',MAXERR*rhsnorm
        write(11,*) ' RHS norm:             ',rhsnorm
        write(11,*) ' Maximum number of iterations: ',MAXIT
        write(11,*) ' '
!       write(11,*) 'A list of iterate numbers and corresponding residual'
!       write(11,*) 'norms is located in the file: ',ITERFILE
!       write(11,*) ' '

        write(11,99) ' Number of unknowns:  ', ndf
        write(11,*) '  '
 99     format(1x,a21,1x,i8)
!
!  Initialize some vectors
!
        p = (0.d0,0.d0)
        d = (0.d0,0.d0)
        rlast = 0d0
!
!  Squirrel away the RHS vector into the vector 's'
!  since 's' isn't used until the main loop.  Store in
!  the vector 'rhs' the initial guess at the solution.
!
        s = rhs
        rhs = (0.d0,0.d0)
!----------------------------! ATTENTION !-----------------------------!
!                                                                      !
!  From this point forward in the algorithm, the vector 'rhs' will     !
!  contain the current solution.                                       !
!                                                                      !
!----------------------------------------------------------------------!
!
!  Compute the initial residual vector.
!
!   这里直接赋值等于不行么？  何必浪费计算资

        call inprod(ndf,rhs,tmp,matflops)

        r = s - tmp              ! 这里直接等于rhs吧
        rho_np = sum(r*dconjg(r))
        rho_np = dsqrt(rho_np)

        write(11,*) ' Initial residual using zero: ',rho_np
!
!  Now compute residual using guess solution
!
        rho_np0 = rho_np
        rhs = guess

        call inprod(ndf,rhs,tmp,matflops)

        r = s - tmp
        rho_np = sum(r*dconjg(r))
        rho_np = dsqrt(rho_np)

        write(11,*) ' Initial residual using guess: ',rho_np
        write(11,*) ' '
!
!  Check for both residuals =0, meaning 0 vector solution:
!
        if ((rho_np.eq.0).and.(rho_np0.eq.0)) then
           write(11,*) ' QMR: trivial 0 vector solution, returning'
           close(11)
           return
        endif
!
! If 0 vector residual better than guess use it.  Make
! sure this is not a smoothing iteration where rhs=0
!
        if ( (rho_np.gt.rho_np0).and.(rhsnorm.ne.0) ) then
           ! zero vector better than guess, reinitialize
           write(11,*) ' Guess not as good, using zero vector instead...'
           rhs = (0.d0,0.d0)
           call inprod(ndf,rhs,tmp,matflops)
           r = s - tmp
           rho_np = sum(r*dconjg(r))
           rho_np = dsqrt(rho_np)
        endif
!
!  Initialize these chaps...
!
        vt       = s - tmp
        v        = vt/rho_np
        c_n      =  1.d0
        eps      = ( 1.d0, 0.d0)
        theta_n  =  0.d0
        eta      = (-1.d0, 0.d0)
!
!  ...and set the vector 's' to its proper starting value.
!
        s = (0.d0,0.d0)

!----------------------------------------------------------------------!
!                       Start main iterative loop                      !
!----------------------------------------------------------------------!

        write(11,*) ' Now computing the QMR iterates...'
        write(11,*) ' '

        do 100 iter=1,MAXIT
           c_nm     = c_n
           rho_n    = rho_np
           theta_nm = theta_n
!
! step 1, flops: 3*ndf
!
!          t11  = etime_(time) !<------------------------------------------! T
                                                                       ! I
           del  = sum(v*v*inv_diag)                                       ! M
                                                                       ! E
!          t22  = etime_(time) !<------------------------------------------! R
!          t(1) = t22-t11

           if (del.eq.(0.d0,0.d0)) then
              if (iter.eq.1) then
                 write(11,*) ' QMR: trivial solution... returning to main routine.'
                 close(11)
                 return
              else
                 write(11,*) ' fatal error: del = 0'
                 stop
              endif
           endif
!
! step 2, flops: 2 + 3*ndf
!
!          t11  = etime_(time) !<------------------------------------------! T
                                                                       ! I
           k1 = rho_n * del / eps                                         ! M
           p  = v*inv_diag - p*k1                                         ! E
                                                                       ! R
!          t22  = etime_(time) !<------------------------------------------!
!          t(2) = t22-t11
!
! step 3, flops: matflop + 2 * 8*ndf
!
!          t11  = etime_(time) !<------------------------------------------! T
                                                                       ! I
           call inprod(ndf,p,tmp,matflops)                                ! M
                                                                       ! E
!          t22  = etime_(time) !<------------------------------------------! R
!          t(3) = t22-t11

!          t11  = etime_(time) !<------------------------------------------! T
                                                                       ! I
           eps  = sum(p*tmp)                                              ! M
           beta = eps / del                                               ! E
           vt   = tmp - v*beta                                            ! R
           rho_np = sum(vt*dconjg(vt))                                    !
           rho_np = dsqrt(rho_np)                                         !
                                                                       !
!          t22  = etime_(time) !<------------------------------------------!
!          t(4) = t22-t11
!
! step 4, flops: 13 + 8*ndf
!
!          t11  = etime_(time) !<------------------------------------------! T
                                                                       ! I
           theta_n = rho_np*rho_np/c_nm/(dble(beta)**2+dimag(beta)**2)    ! M
           c_n = 1.d0/(1.d0 + theta_n)                                    ! E
           eta = -eta*rho_n*c_n/(beta*c_nm)                               ! R
           k2  = theta_nm*c_n                                             !
                                                                       !
           d   =   p*eta + d*k2                                           !
           rhs = rhs     + d                                              !
           s   = tmp*eta + s*k2                                           !
           r   = r       - s                                              !
                                                                       !
!          t22  = etime_(time) !<------------------------------------------!
!          t(5) = t22-t11
!
! step 5, flops: 1 + 5*ndf
!
!          t11  = etime_(time) !<------------------------------------------! T
                                                                       ! I
           v   = vt/rho_np                                                ! M
           rnorm = sum(r*conjg(r))                                        ! E
           rnorm = dsqrt(rnorm)                                            ! R
                                                                       !
!          t22  = etime_(time) !<------------------------------------------!
!          t(6) = t22-t11
!
!  Every 50 iterations, write out progress to the screen
!
           if (mod(iter,50).eq.0) then
!             t2 = etime_(time)
!             write(11,98) iter, MAXIT,rnorm,t2-t1
	      write(11,981) iter, MAXIT,rnorm
           endif
!
!  Every iteration, write the iteration count and the residual norm to
!  the output file.
!
!          write(11,*) iter,rnorm
!
!  Check to see if the solution has converged...
!
           if (rnorm.lt.MAXERR*rhsnorm) then
!             t2 = etime_(time)
              k2 = flops(iter,matflops,ndf)

              write(11,*) ' Solution converged in ',iter,' iterations'
              write(11,*) '  Residual:',rnorm
              write(11,*) ' '
!             write(11,*) ' Time spent in QMR [s]: ',t2 - t1
!             write(11,*) ' Seconds per iteration:    ',(t2-t1)/real(iter)
!             write(11,*) '  '
!             write(11,*) ' Flop Rate [Mflop/s]:',k2/(t2-t1)/1.e6

              write(11,99) ' Number of unknowns:  ', ndf
              write(11,*) '  '
!
!  Deallocate storage for the local arrays used in the QMR routine.
!
              deallocate (inv_diag,p,d,vt,v,tmp,r,s)
!
!  Au revoir, QMR.  It's back to the main routine...
!
              close(11)
!             write(6,*) 'close'
              return
!
! KWK 05/02/2006:
! Catch to stop if residual not reducing for ntol iterations
!
           else
              if (abs(rnorm-rlast)/rnorm.le.tol) then  ! if changed only a little
                 itol = itol+1
                 if (itol.gt.ntol) then
                    write(11,*) ' '
                    write(11,*) ' QMR residual not changing for last ',itol,' iterations'
                    write(11,*) ' Returning ....'
                    write(11,*) ' Residual:',rnorm
                    write(11,*) ' '
!                    write(11,*) ' Time spent in QMR [s]: ',t2 - t1
!                    write(11,*) ' Seconds per iteration:    ',(t2-t1)/real(iter)
!                    write(11,*) ' '
!                    write(11,*) ' Flop Rate [Mflop/s]:',k2/(t2-t1)/1.e6
                    write(11,*) ' Number of unknowns:  ', ndf
                    close(11)
!                   write(6,*) 'close'
                    deallocate (inv_diag,p,d,vt,v,tmp,r,s)
                    return
                 endif
              else
                 itol = 0
              endif
              rlast = rnorm
           endif
 100    continue
!
!  If the residual is still unacceptibly large after the maximum
!  number of iterations, then just exit gracefully.
!
!       t2 = etime_(time)
        k2 = flops(iter,matflops,ndf)
!       write(11,*) ' '
!       write(11,*) '********** WARNING ************'
        write(11,*) ' Used Maximum # iterations    '
!       write(11,*) '*******************************'
!       write(11,*) ' '
        write(11,*) ' Returning best solution so far...'
        write(11,*) ' Iteration:',iter
        write(11,*) ' Residual:',rnorm
        write(11,*) ' Target residual norm: ',MAXERR*rhsnorm
        if (iter.gt.100) then
           write(11,*) '  ... better luck next time, bucko.'
        endif
!       write(11,*) ' '
!       write(11,*) ' Time spent in QMR [s]: ',t2 - t1
!       write(11,*) ' Seconds per iteration:    ',(t2-t1)/real(iter)
!       write(11,*) ' '
!       write(11,*) ' Flop Rate [Mflop/s]:',k2/(t2-t1)/1.e6
        write(11,99) ' Number of unknowns:  ', ndf
        write(11,*) '  '

        close(11)
!       write(6,*) 'close'
!
!  Deallocate storage for the local arrays used in the QMR routine.
!
        deallocate (inv_diag,p,d,vt,v,tmp,r,s)

! 98    format('  it: ',i5,'/',i5,'    r: ',e12.6,'    t:',f8.4,1x,f8.2)
 981    format('  it: ',i5,'/',i5,'    r: ',e13.6,'   t:',f8.4 )

 101    format(10(e11.4,1x))

        end subroutine qmr
!      
!======================================================================!
!=============================================================== FLOPS !
!======================================================================!
!
        function flops(i,matflops,ndf)
        use       fem_params

        real(double) :: flops,matflops
        integer      :: i,ndf
!
!  Define a function that returns the number of floating point
!  operations as a function of QMR iterates.
!
        flops = 3.d0 + matflops + 11.d0*ndf +    &
                 i*(18.d0 + matflops + 27.d0*ndf)
        end function flops
!
!======================================================================!
!============================================================= indexx !
!======================================================================!
!
        SUBROUTINE indexx(n,arr,indx)
!
! Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that 
! arr(indx(j)) is in ascending order for j=1,2,...,N. 
! The input quantities n and arr are not changed.
!
        implicit none

        INTEGER n,indx(n),M,NSTACK
        REAL(8) arr(n)
        PARAMETER (M=7,NSTACK=1000)
        INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
        REAL(8) a

        do j=1,n
           indx(j)=j
        enddo

        jstack=0
        l=1
        ir=n

1       if(ir-l.lt.M)then
           do 13 j=l+1,ir
              indxt=indx(j)
              a=arr(indxt)
              do 12 i=j-1,l,-1
                 if(arr(indx(i)).le.a)goto 2
                 indx(i+1)=indx(i)
12            continue
              i=l-1
2             indx(i+1)=indxt
13         continue

           if(jstack.eq.0)return
              ir=istack(jstack)
              l=istack(jstack-1)
              jstack=jstack-2
        else
           k=(l+ir)/2
           itemp=indx(k)
           indx(k)=indx(l+1)
           indx(l+1)=itemp

           if(arr(indx(l)).gt.arr(indx(ir)))then
              itemp=indx(l)
              indx(l)=indx(ir)
              indx(ir)=itemp
           endif

        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif

        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif

        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)

3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3

4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4

        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)

        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'

        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif

        endif
        goto 1
        END
!
!
!==================================================================================
!==================================================================================
subroutine ICBICG_STAB_CSR(ND,NNZ,A,icol,irow,B)
!***********************************************************************************
!***********************************************************************************
!**          ±Ÿ³ÌÐò²Î¿Œ¡¶An iterative solution method for linear systems
!**    of which the coefficient matrix is a symmetric M-matrix¡·By J. A. Meijerink 
!**              ²ÉÓÃ²»ÍêÈ«CholeskyÔ€ŽŠÀí  ÁõÓ±(2010-7-28)    
!**                    °æÈš£º hi.baidu.com/aliouying
!**  
!***********************************************************************************
!***********************************************************************************
!**    modified to ICBICG_ATAB(Biconjugate gradient stabilized method using an 
!**           imcomplete Cholesky preconditioned technique)
!**                   by Yanbo(2015.1.6)  OUC
!**
!***********************************************************************************
!***********************************************************************************
!** AŸØÕóŽæŽ¢žñÊœ£º°ŽÐÐË³ÐòÑ¹ËõŽ¢ŽæÉÏÈýœÇŸØÕó
!** A(NNZ)save nonzero elements  of upper triangular matrix A by row;
!** icol(nnz)save the column number of nonzero elements in A(nnz) by row;
!** irow(nd+1)save the index values of diagonal elements(nonzero elements)
!** of input upper triangular matrix A by row.
!** The storage method described above is calle the row compression storage technology (CSR)
!** B(nd)save the right-hand side,and nrhs is the colunm number of the right-hand side.
implicit none
!***********************************************************************************
!**   ÊäÈëºÍÊä³ö
!***********************************************************************************
integer,intent(in)       :: ND,NNZ
integer,intent(in)       :: icol(NNZ),irow(ND+1) !CSRÑ¹ËõŽæŽ¢Ë÷ÒýŸØÕópointB,pointE
complex(8),intent(in)    :: A(NNZ)               !ŸØÕó±äÁ¿A*x=B
complex(8),intent(inout) :: B(ND)           !ŸØÕó±äÁ¿A*x=B
!***********************************************************************************
!**   ³ÌÐòÖÐŒä±äÁ¿
!***********************************************************************************
complex(8) ,dimension(:),allocatable  :: R,RR
complex(8) ,dimension(:),allocatable  :: P,PP
complex(8) ,dimension(:),allocatable  :: S,SS
complex(8) ,dimension(:),allocatable  :: v,t
complex(8) ,dimension(:),allocatable  :: X   
complex(8) ,dimension(:),allocatable  :: D 
real(8),dimension(:),allocatable      :: res_err
complex(8), external                  :: m_dot_product
integer    :: i,inrhs,j
integer    :: IT                                !ÊÕÁ²×ÜŒÆµüŽúŽÎÊý
real(8)    :: EPS                               !ÊÕÁ²Îó²î
integer    :: ITMAX                             !the maximal iterative number
real(8)    :: res_err0                          !the initial residual error before solving equation.
complex(8) :: ro,ro1
complex(8) :: alpha,womega,belta
real(8)    :: time0,time1
!***********************************************************************************
!***********************************************************************************
!**   ÉèÖÃ³õÊŒÖµ(³õÊŒ»¯±äÁ¿)
 EPS = 1.0D-20 
 ITMAX=2000
 allocate(res_err(ITMAX))
 allocate(R(nd),RR(nd),P(ND),PP(nd),S(ND),SS(nd),v(ND),t(nd),D(ND),X(ND))
 X=(0.d0,0.d0)                                       !ÊäÈë³õÊŒ¹ÀŒÆÖµºÍÊä³ö×îÖÕœá¹û
!***********************************************************************************
!***********************************************************************************
call cpu_time(time0)
call Creat_diagonal_C(ND,NNZ,A,D,icol,irow)
call cpu_time(time1)
write(*,fmt='(5x,a40,f15.7)') 'The time of forming the diagonal matrix:', (time1-time0)
!*************************************************************************************
!**   ¿ªÊŒŒÆËã
!*************************************************************************************
!write(*,fmt='(5x,a40,i10)') 'The sum number of the right-hand side', nrhs 
!do inrhs=1,nrhs                                        !loop the right-hand side 
!    write(*,fmt='(5x,a40,i10)') 'The loop number of the right-hand side',inrhs 
    call MatMulVec_CSR(ND,NNZ,A,X,R,icol,irow)     
    R = B- R
    !r(0)=b-A*x(0)

    IT=0
    res_err0=sqrt(sum(abs(R)**2))
    write(*,fmt='(10x,a20)')  'the initial residual error'
    write(*,fmt='(10x,e20.10)') res_err0
    !if( res_err0< EPS ) then
    !    write(*,*) 'The initial input solution is the final solution'
    !    b=x
    !    exit
    !endif
    RR=R
    ro=m_dot_product(RR,R,ND)
    write(*,*) abs(ro)

    if(abs(ro) <1.0d-25) then
    write(*,*) 'This algorithm is failed,please check input parameters or change other iterative method.'
    write(*,*) 'The ro is equal zero',ro
    deallocate(res_err)
    deallocate(R,RR,S,SS,P,PP,v,T,X,D)
    stop
    end if

    P=R
    do i = 1,ITMAX
        IT=i  !record the iterative numbers
        call forward_solve(ND,NNZ,A,D,P,PP,icol,irow)
        !C(T)D(-1)C(T)PP=P
        forall(j=1:ND)              
        PP(j) = PP(j)*D(j)       
        end forall
     
        !ÕâÀï¿ÉÒÔÖ±œÓÓÃZ=V*YŽúÌæ£¬Ö÷ÒªÄ¿µÄÊÇÎªÁË·œ±ãÌæ»»V,¶øÇÒforallµÄÐ§ÂÊ²»Ží^_^
         call back_solve(ND,NNZ,A,D,PP,PP,icol,irow)
       
        !C(T)PP=PP
         call MatMulVec_CSR(ND,NNZ,A,PP,v,icol,irow)
         alpha = ro/m_dot_product(RR,v,ND)
         S=R - alpha*v
     
         res_err(i) = sqrt(sum(abs(S)**2))
      
        if( res_err(i) < EPS ) then
             X = X + alpha*PP
             b=X
            write(*,fmt='(10x,a20,a20)') 'the iterative number','the residual error'
            write(*,fmt='(10x,I15,5x,e20.10)') IT,res_err(i)
            exit
        endif

          call forward_solve(ND,NNZ,A,D,S,SS,icol,irow)

          forall(j=1:ND)              
           SS(j) = SS(j)*D(j)        
          end forall

         !ÕâÀï¿ÉÒÔÖ±œÓÓÃZ=V*YŽúÌæ£¬Ö÷ÒªÄ¿µÄÊÇÎªÁË·œ±ãÌæ»»V,¶øÇÒforallµÄÐ§ÂÊ²»Ží^_^
        call back_solve(ND,NNZ,A,D,SS,SS,icol,irow)
        call MatMulVec_CSR(ND,NNZ,A,SS,T,icol,irow)
        womega= m_dot_product(T,S,ND)/m_dot_product(T,T,ND)
       
        if(abs(womega) < 1.0d-25) then
        write(*,*) 'This algorithm is failed,please check input parameters or change other iterative method.'
        write(*,*) 'The womega is equal zero',womega
        deallocate(res_err)
        deallocate(R,RR,S,SS,P,PP,v,T,X,D)
        stop
        end if

        X = X+ alpha*PP+womega*SS
        R=S-womega*T
        res_err(i) = sqrt(sum(abs(R)**2))

        if( res_err(i) < EPS ) then
            b=x
            write(*,fmt='(10x,a20,a20)') 'the iterative number','the residual error'
            write(*,fmt='(10x,I15,5x,e20.10)') IT,res_err(i)
            exit
        endif

        ro1=m_dot_product(RR,R,ND)

        if(abs(ro1)<1.0d-100) then
            write(*,*) 'This algorithm is failed,please check input parameters or change other iterative method.'
            write(*,*) 'The ro1 is equal zero',ro1
            deallocate(res_err)
            deallocate(R,RR,S,SS,P,PP,v,T,X,D)
            stop
        end if

        belta=(ro1/ro)*(alpha/womega)
        P=R+belta*(P-womega*V)
        ro= ro1

        write(*,fmt='(10x,a20,a20)') 'the iterative number','the residual error'
        write(*,fmt='(10x,I15,5x,e20.10)') IT,res_err(i)

        if(i==ITMAX) then
            b=x
            write(*,*) 'The solution is not convergent'
            write(*,fmt='(10x,a20,a20)') 'the iterative number','the residual error'
            write(*,fmt='(10x,I15,5x,e20.10)') IT,res_err(i)
            deallocate(res_err)
            deallocate(R,RR,S,SS,P,PP,v,T,X,D)
            stop
        end if 
        
    enddo
!end do   
                                !end right-hand side
deallocate(res_err)
deallocate(R,RR,S,SS,P,PP,v,T,X,D)
return
end subroutine ICBICG_STAB_CSR
!================this is split line================================================!
subroutine Creat_diagonal_C(ND,NNZ,MAT,VEC_out,icol,irow)
implicit none
!***********************************************************************************
!**   ÊäÈëºÍÊä³ö
!***********************************************************************************
integer,intent(in)    :: ND,NNZ
integer,intent(in)    :: icol(NNZ),irow(ND+1)       !CSRÑ¹ËõŽæŽ¢Ë÷ÒýŸØÕópointB,pointE
complex(8),intent(in) :: MAT(NNZ)                   !ŸØÕó±äÁ¿MAT
complex(8)            :: VEC_out(ND)                !Êä³ö×îÖÕœá¹û
!***********************************************************************************
!**   ³ÌÐòÖÐŒä±äÁ¿
!***********************************************************************************
integer :: i,k
VEC_out=0.d0
VEC_out(1) = MAT(irow(1))
do i = 2 , ND
        do k = irow(i-1)+1,irow(i)-1
             VEC_out(icol(k)) = VEC_out(icol(k)) + MAT(k)*MAT(k)/VEC_out(i-1)
        enddo
    VEC_out(i) = MAT(irow(i)) - VEC_out(i)
enddo

return
end subroutine Creat_diagonal_C
!================this is split line================================================!
subroutine MatMulVec_CSR(ND,NNZ,MAT,VEC,VEC_out,icol,irow)
implicit none
!***********************************************************************************
!**   ÊäÈëºÍÊä³ö
!***********************************************************************************
integer,intent(in)    :: ND,NNZ
integer,intent(in)    :: icol(NNZ),irow(ND+1)     !CSRÑ¹ËõŽæŽ¢Ë÷ÒýŸØÕópointB,pointE
complex(8),intent(in) :: MAT(NNZ),VEC(ND)         !ŸØÕó±äÁ¿MATºÍÏòÁ¿VEC
complex(8)            :: VEC_out(ND)              !Êä³ö×îÖÕœá¹û
!***********************************************************************************
!**   ³ÌÐòÖÐŒä±äÁ¿
!***********************************************************************************
integer :: i,j
VEC_out = (0.0d0,0.0d0)
do i = 1 , ND
     VEC_out(i) = VEC_out(i) + MAT(irow(i))*VEC(i)
    do j = irow(i)+1, irow(i+1)-1
        VEC_out(i) = VEC_out(i) + MAT(j)*VEC(icol(j))
        VEC_out(icol(j)) = VEC_out(icol(j)) + MAT(j)*VEC(i)
    enddo
enddo
return
end subroutine MatMulVec_CSR
!==================================================================================
!================this is split line================================================!
function m_dot_product(A_in,B_in,N)
implicit none
integer    :: N
complex(8) :: A_in(N), B_in(N)
complex(8) :: m_dot_product
integer    :: i
m_dot_product = 0.0
do i = 1 , N
    m_dot_product = m_dot_product + A_in(i)*B_in(i)
enddo
return
end function m_dot_product
!================this is split line================================================!
!==================================================================================
!================this is split line================================================!
subroutine forward_solve(ND,NNZ,MAT,D,VEC,VEC_out,icol,irow)
implicit none
!***********************************************************************************
!**   ÊäÈëºÍÊä³ö
!***********************************************************************************
integer,intent(in)      :: ND,NNZ
integer,intent(in)      :: icol(NNZ),irow(ND+1)         !CSRÑ¹ËõŽæŽ¢Ë÷ÒýŸØÕópointB,pointE
complex(8),intent(in)   :: MAT(NNZ),VEC(ND)             !ŸØÕó±äÁ¿MATºÍÏòÁ¿VEC
complex(8),intent(in)   :: D(nd)                        !diagonal matrix
complex(8),intent(out)  :: VEC_out(ND)                  !Êä³ö×îÖÕœá¹û
!***********************************************************************************
!**   ³ÌÐòÖÐŒä±äÁ¿
!***********************************************************************************
integer    :: i,j
complex(8) :: sum
VEC_out = VEC
do i = 1,ND
   VEC_out(i) = VEC_out(i)/D(i)
   do j = irow(i)+1,irow(i+1)-1    
       VEC_out(icol(j)) = VEC_out(icol(j)) - MAT(j)*VEC_out(i)
   enddo
enddo
return
    end subroutine forward_solve
    
!==================================================================================
!================this is split line================================================!
 subroutine back_solve(ND,NNZ,MAT,D,VEC,VEC_out,icol,irow)
implicit none
!***********************************************************************************
!**   ÊäÈëºÍÊä³ö
!***********************************************************************************
integer,intent(in)     :: ND,NNZ
integer,intent(in)     :: icol(NNZ),irow(ND+1)     !CSRÑ¹ËõŽæŽ¢Ë÷ÒýŸØÕópointB,pointE
complex(8),intent(in)  :: MAT(NNZ),VEC(ND)         !ŸØÕó±äÁ¿MATºÍÏòÁ¿VEC
complex(8),intent(in)  :: D(nd)                    !diagonal matrix
complex(8),intent(out) :: VEC_out(ND)              !Êä³ö×îÖÕœá¹û
!***********************************************************************************
!**   ³ÌÐòÖÐŒä±äÁ¿
!***********************************************************************************
integer    :: i,j
complex(8) :: sum


do i = ND,1,-1
   sum = (0.0d0,0.0d0)
   
   do j = irow(i+1)-1,irow(i)+1,-1
       sum = sum + MAT(j)*VEC_out(icol(j))
   enddo
   VEC_out(i) = (VEC(i)-sum)/D(i)
enddo
return
    end subroutine back_solve

!=========================================================================================
!========================================================================================
subroutine SSOR_PCG_CSR(ND,NNZ,A,B,col,irw)
!***********************************************************************************
!**   Â±ÅžÂ³ÃÃÃ²Â²ÃÂ¿ÅÃÃÃÃÃÃÂ¡Â¶ÃÃÃâ¬ÅœÅ ÃÃ­Â¹Â²Ã©Ã®ÃÃÂ¶ÃÂ·Å¡ÃÃ³ÅÃ¢ÃÃÃÃÃÂªÂ·ÅÂ³ÃÃÃ©ÅÂ°Â³ÃÃÃ²ÃÃšÅÃÂ¡Â·ÅŸÃÅÃžÂµÃÂµÃŒÅœÃºÅŸÃ±ÃÅ
!**    Â²ÃÃÃSSORÃâ¬ÅœÅ ÃÃ­Â£Â¬Ãâ¬ÅœÅ ÃÃ­ÅžÃÃÃ³ÃÂª:M=(2-W)^{-1}*(D/W+L)*(D/W)^{-1}*(D/W+L)^{T}
!**                   ÃÃµÃÂ±     liouying@ouc.edu.cn
!**                          2013-7-31
!***********************************************************************************
 implicit none
!***********************************************************************************
!**   ÃÃ€ÃÃ«ÂºÃÃÃ€Â³Ã¶
!***********************************************************************************
 integer,intent(in) :: ND,NNZ
 integer,intent(in) :: col(NNZ),irw(ND+1) ! CSRÃÂ¹ÃÃµÅœÃŠÅœÂ¢ÃÃ·ÃÃœÅžÃÃÃ³pointB,pointE
 complex(8),intent(in) :: A(NNZ)    !ÅžÃÃÃ³Â±Ã€ÃÂ¿A*x=B
 complex(8),intent(inout) :: B(ND)         !ÊäÈë³õÊŒ¹ÀŒÆÖµºÍÊä³ö×îÖÕœá¹û
! real(8),intent(in) :: EPS                 !ÃÃÃÂ²ÃÃ³Â²Ã®
! integer, intent(in) :: ITMAX              !ÃÃ®ÅœÃ³ÂµÃŒÅœÃºÅœÃÃÃœ
! real(8), intent(in) :: W                  !ÃÃÂ³ÃÃÃ²ÃÃ
! integer, intent(inout) :: IT                !ÃÃÃÂ²ÃÃÅÃÂµÃŒÅœÃºÅœÃÃÃœ
! real(8), intent(out) :: reg_err           !ÂµÃŒÅœÃºÂºÃ³ÂµÃÂ²ÃÂ²Ã®

        real(8) :: W,EPS,reg_err
        INTEGER :: IT,ITMAX
        complex(8):: X(ND)        
!***********************************************************************************
!**   Â³ÃÃÃ²ÃÃÅÃ€Â±Ã€ÃÂ¿
!***********************************************************************************
 complex(8) :: Y(ND),Z(ND),D(ND),V(ND),TEMP(ND)
 integer :: i
 complex(8) :: delta, tao, beta,delta1,delta0
 complex(8), external :: m_dot_product

        ITMAX = 5000
        EPS = 1.0D-18                      
        W = 1.d0
        X =(0.d0,0.d0) 
!***********************************************************************************
!**   ÃÃšÃÃÂ³ÃµÃÅÃÂµ
!***********************************************************************************

 forall(i=1:ND)
    V(i) = (2.0-W)*A(irw(i))/W      !ÃÃ¢ÃÃ¯VÂµÂ±Â³ÃÃÃÃÂ»ÃÂ¬ÃÃ²ÃÂ¿Â£Â¬ÃÂµÅÃÃÃÃÃŒÃÃÃÂ»ÅŸÃ¶Â¶ÃÅÃÅžÃÃÃ³Â£Â¬ÅÂ°
 end forall 

! open(10,file='A.dat',status='REPLACE')
!  do i=1,ND
!    write(10,*)v(i)
!  enddo
! close(10)
!***********************************************************************************
!**   Â¿ÂªÃÅÅÃÃÃ£
!***********************************************************************************
 call MatMulVec_CSR(ND,NNZ,A,X,Y,col,irw)
 Y = Y - B

 call forward_solve_w(ND,NNZ,A,Y,Y,col,irw,w)

 forall(i=1:ND)              !ÃÃ¢ÃÃ¯Â¿ÃÃÃÃÂ±ÅÃÃÃZ=V*YÅœÃºÃÃŠÂ£Â¬ÃÃ·ÃÂªÃÂ¿ÂµÃÃÃÃÂªÃÃÂ·ÅÂ±Ã
    Z(i) = -V(i)*Y(i)      
 end forall

! open(12,file='Y.dat',status='REPLACE')
!  do i=1,ND
!    write(12,*)Y(i)
!  enddo
! close(12)

 call back_solve_w(ND,NNZ,A,Z,D,col,irw,w)

 open(unit=11,file='ssorpcg.iter',position='APPEND')

 IT = 0
 delta = m_dot_product(Y,V*Y,ND)

   delta0=delta
!   write(*,*) 'the initial residual error:',delta0

 beta = 1
 do i = 1,ITMAX
!    if( abs(delta) < EPS) then
!        return
!    endif
    tao = delta/m_dot_product(D,2.0*Z-V*D,ND)
    X = X + tao*D
! Â¶Ã ÃÃÃÃÃÂ²Â·ÅÃÅ
    !if( sqrt(sum(abs(tao*D/X)**2)) < EPS ) then
    !    return
    !endif
    call forward_solve_w(ND,NNZ,A,Z-V*D,TEMP,col,irw,w)
    Y = Y + tao*(TEMP+D)
    delta1 = m_dot_product(Y,V*Y,ND)
    beta = delta1/delta
    Z = -V*Y + beta*Z
    call back_solve_w(ND,NNZ,A,Z,D,col,irw,w)
    delta = delta1
    IT = IT + 1
    reg_err = abs(delta/delta0)
!    write(*,*) 'PCG iteration:',IT,'of res:',reg_err 
    write(11,*) IT,reg_err    
    if( reg_err < EPS.or.IT==ITMAX ) Then
        B=X
        return
    Endif
 enddo

 close(11)    
 return
 end subroutine SSOR_PCG_CSR    

!================this is split line================================================!
subroutine forward_solve_w(ND,NNZ,MAT,VEC,VEC_out,jcol,irw,W)
implicit none
!***********************************************************************************
!**   ÊäÈëºÍÊä³ö
!***********************************************************************************
integer,intent(in) :: ND,NNZ
integer,intent(in) :: jcol(NNZ),irw(ND+1)         !CSRÑ¹ËõŽæŽ¢Ë÷ÒýŸØÕópointB,pointE
complex(8),intent(in) :: MAT(NNZ),VEC(ND)         !ŸØÕó±äÁ¿MATºÍÏòÁ¿VEC
real(8), intent(in) :: W
complex(8),intent(out) :: VEC_out(ND)             !Êä³ö×îÖÕœá¹û
!***********************************************************************************
!**   ³ÌÐòÖÐŒä±äÁ¿
!***********************************************************************************
integer :: i,j

VEC_out = VEC

do i = 1,ND
    VEC_out(i) = VEC_out(i)/MAT(irw(i))*W
    do j = irw(i)+1,irw(i+1)-1    
        VEC_out(jcol(j)) = VEC_out(jcol(j)) - MAT(j)*VEC_out(i)
    enddo
enddo

return
end subroutine forward_solve_w

!================this is split line================================================!
subroutine back_solve_w(ND,NNZ,MAT,VEC,VEC_out,jcol,irw,W)
implicit none
!***********************************************************************************
!**   ÊäÈëºÍÊä³ö
!***********************************************************************************
integer,intent(in) :: ND,NNZ
integer,intent(in) :: jcol(NNZ),irw(ND+1)   !CSRÑ¹ËõŽæŽ¢Ë÷ÒýŸØÕópointB,pointE
complex(8),intent(in) :: MAT(NNZ),VEC(ND)   !ŸØÕó±äÁ¿MATºÍÏòÁ¿VEC
real(8), intent(in) :: W
complex(8),intent(out) :: VEC_out(ND)              !Êä³ö×îÖÕœá¹û
!***********************************************************************************
!**   ³ÌÐòÖÐŒä±äÁ¿
!***********************************************************************************
integer :: i,j
complex(8) :: sum

do i = ND,1,-1
    sum = (0.0d0,0.0d0)
    do j = irw(i+1)-1,irw(i)+1,-1
        sum = sum + MAT(j)*VEC_out(jcol(j))
    enddo
    VEC_out(i) = (VEC(i)-sum)/MAT(irw(i))*W
enddo
    
return
end subroutine back_solve_w

    
    
 
