      program cadd_nodisl_k_test
      
      use mod_types, only: dp
      use mod_io, only: writeOutput, initSimulation, writeRestart_ptr,
     &                  writeDump_ptr
      use mod_integrate, only: loopVerlet
      use mod_misc, only: updateMiscIncrementCurr
      use mod_poten_assign, only: getPotForcesAll_ptr
      use mod_neighbors, only: updateNeighborsCheck, 
     &   updateNeighborsNoCheck, updateNeighIncrementCurr
      use mod_math, only: isMultiple
      use mod_kdispfield, only: applyKDispIso
      use mod_materials, only: materials
      use mod_fe_elements, only: fematerials
      use mod_nodes, onlY: nodes
      use mod_damping, only: normaldamping, dampingdata
      use mod_fe_main_2d_assign, only: solveAll_ptr,
     &                                 updateFENodalPosnAll_ptr
      use mod_pad_atoms, only: updatePad
      
      implicit none
      
      integer :: i
      integer :: mnumfe, mnum
      real(dp) :: KIstart, KIincr, KIend, KIapply, KIcurr, KII
      real(dp) :: xc, yc
      real(dp) :: mu, nu
      integer :: nstepsK, natomisticsteps
      real(dp) :: dt, disptol
      
C     read, initialize
      call initSimulation('cadd_nodisl_k_test_small','cadd_nodisl')
      write(*,*) 'Blah'
      
C     material stuff
      mnumfe = 1
      mnum = fematerials%list(mnumfe)
      mu = materials(mnum)%mu
      nu = materials(mnum)%nu
      
C     K-field
      KIstart = 5.0_dp
      KIincr = 0.25_dp
      KIend = 9.0_dp
      KII = 0.0_dp
      nstepsK = nint((KIend - KIstart)/KIincr)
      
C     crack center (slightly offset from atom)
      xc = 0.5_dp
      yc = 0.3_dp
      
C     atomistic stuff
      natomisticsteps = 100
      dt = 0.02_dp
      disptol = 1.0e-5_dp
      
CC     apply K, equilibrate, dump
C     do i = 1, nstepsK
C         if (i == 0) then
C             KIapply = KIstart
C         else
C             KIapply = KIincr
C         end if
C         
C         KIcurr = KIstart + i*KIincr
C         write(*,*) 'Current KI', KIcurr
C         
C         call applyKDispIso(KIapply,KII,mu,nu,xc,yc,'all')
C         call equilibrateCADDNoDisl(natomisticsteps,dt,normaldamping,
C    &                               disptol)
C         call updateMiscIncrementCurr(1)
C         call writeDump_ptr()
C     end do

      contains
************************************************************************
      subroutine equilibrateCADDNoDisl(natomisticsteps,dt,damp,
     &                                 disptol)

C     input variables
      integer :: natomisticsteps
      real(dp) :: dt
      type(dampingdata) :: damp
      real(dp) :: disptol      

C     local variables
      real(dp), allocatable :: allposnold(:,:)
      real(dp), allocatable :: allposnnew(:,:)
      real(dp) :: dispnorm

C     we will check nodal positions to see if we're done
      allposnold = nodes%posn(1:2,:)    
      dispnorm = huge(0.0_dp)
      
      do while (dispnorm < disptol)
          write(*,*) 'Starting iteration'
          
          call loopVerlet(natomisticsteps,dt,'all',damp) ! step 1 (see Algorithm.txt)
          call solveAll_ptr() ! step 3
          call updatePad() ! step 4
          
C         need to update positions to compare
          call updateFENodalPosnAll_ptr()
          allposnnew = nodes%posn(1:2,:)
          dispnorm = maxval(sum(abs(allposnnew-allposnold),2)) ! infinity-norm
          
          write(*,*) 'Displacement inf. norm', dispnorm
          allposnold = allposnnew
      end do
      
      end subroutine
************************************************************************
      end program