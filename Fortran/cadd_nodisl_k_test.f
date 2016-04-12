      program cadd_nodisl_k_test
      
      use mod_types, only: dp
      use mod_io, only: writeOutput, initSimulation, writeRestart_ptr,
     &                  writeDump_ptr
      use mod_integrate, only: loopVerlet
      use mod_misc, only: misc, updateMiscIncrementCurr
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
      use mod_utils, only: prettyPrintMat, prettyPrintVec
      use mod_find_crack_atomistic, only: findCrack
      
      implicit none
      
      integer :: iunit
      integer :: i
      integer :: mnumfe, mnum
      real(dp) :: KIstart, KIincr, KIend, KIapply, KIcurr, KII
      real(dp) :: xc, yc
      real(dp) :: mu, nu
      integer :: nstepsK, natomisticsteps, natomisticstepstot
      real(dp) :: dt, forcetol
      character(len=15) :: gammasuffix, dtsuffix, stepssuffix
      character(len=:), allocatable :: filename
      real(dp) :: crackpos(2)
      
C     read, initialize
      call initSimulation('cadd_nodisl_k_test_medium','cadd_nodisl')
      call writeDump_ptr()
      write(*,*) 'Wrote dump'
      
C     material stuff
      mnumfe = 1
      mnum = fematerials%list(mnumfe)
      mu = materials(mnum)%mu
      nu = materials(mnum)%nu
      
C     K-field
      KIstart = 8.55_dp
      KIincr = 0.1_dp
      KIend = 9.55_dp
      KII = 0.0_dp
      nstepsK = nint((KIend - KIstart)/KIincr)
      
C     crack center (slightly offset from atom)
      xc = 0.5_dp
      yc = 0.3_dp
      
C     atomistic stuff
      natomisticsteps = 20
      dt = 0.02_dp
      normaldamping%gamma = 0.1
      forcetol = 1.0e-4_dp
      
C     file
      write (gammasuffix,'(I0)') nint(normaldamping%gamma*100.0_dp)
      write (dtsuffix,'(I0)') nint(dt*1000)
      write (stepssuffix,'(I0)') natomisticsteps
      filename = trim(misc%simname)//'_gamma_'//trim(gammasuffix)
     &                             //'_dt_'//trim(dtsuffix)
     &                             //'_steps_'//trim(stepssuffix)
      open(newunit=iunit,file=filename)

      crackpos = [0.0_dp,0.0_dp]      
      
C     apply K, equilibrate, dump
      do i = 0, nstepsK          
          if (i == 0) then
              KIapply = KIstart
          else
              KIapply = KIincr
          end if    
          KIcurr = KIstart + i*KIincr
          write(*,*) 'Current KI', KIcurr
          
          call applyKDispIso(KIapply,KII,mu,nu,crackpos(1),
     &                                         crackpos(2),'all')
          call equilibrateCADDNoDisl(natomisticsteps,dt,normaldamping,
     &                               forcetol,natomisticstepstot)
          write(iunit,*) KIcurr, natomisticstepstot
          call updateMiscIncrementCurr(1)
          call writeDump_ptr()
          crackpos = findCrack()
          write(*,*) 'Crack position', crackpos
      end do
      
      close(iunit)

      contains
************************************************************************
      subroutine equilibrateCADDNoDisl(natomisticsteps,dt,damp,
     &                                 forcetol,natomisticstepstot)

C     input variables
      integer :: natomisticsteps
      real(dp) :: dt
      type(dampingdata) :: damp
      real(dp) :: forcetol
      
C     output variables
      integer :: natomisticstepstot

C     local variables
      real(dp) :: forcenorm
      integer :: counter

C     we will check forces on atoms
      forcenorm = huge(0.0_dp)
      counter = 0
      
      do while (forcenorm > forcetol)
          ! write(*,*) 'Starting iteration'
          call loopVerlet(natomisticsteps,dt,'all',damp) ! step 1 (see Algorithm.txt)
          counter = counter + 1
          
          call solveAll_ptr() ! step 3
          call updatePad() ! step 4
          
          forcenorm = maxval(sum(abs(nodes%potforces),1)) ! infinity norm
          write(*,*) 'forcenorm', forcenorm
      end do
      natomisticstepstot = counter*natomisticsteps
      
      end subroutine equilibrateCADDNoDisl
************************************************************************
      end program