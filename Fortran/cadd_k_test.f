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
      use mod_disl_detect_pass, only: detectAndPassDislocations
      use mod_dd_main, only: runDDStep
      use mod_dd_integrate, only: getResolvedStressOnDisl
      use mod_disl_try, only: disl
      use mod_disl_escaped, only: escapeddisl
      use mod_disl_ghost, only: ghostdisl
      
      implicit none
      
      integer :: iunit
      integer :: i
      integer :: mnumfe, mnum
      real(dp) :: KIstart, KIincr, KIend, KIapply, KII
      real(dp) :: xc, yc
      real(dp) :: mu, nu
      integer :: nstepsK, natomisticsteps, natomisticstepstot
      real(dp) :: dt, dtdd, forcetol
      character(len=15) :: gammasuffix, dtsuffix, stepssuffix
      character(len=:), allocatable :: filename
      integer :: counter
      
C     read, initialize
      call initSimulation('cadd_k_test_large','cadd')
      call writeDump_ptr()
      write(*,*) 'Wrote dump'
      
C     material stuff
      mnumfe = 1
      mnum = fematerials%list(mnumfe)
      mu = materials(mnum)%mu
      nu = materials(mnum)%nu
      
C     K-field
      KIstart = 7.5_dp
      KIincr = 0.25_dp
      KIend = 9.5_dp
      KII = 0.0_dp
      nstepsK = nint((KIend - KIstart)/KIincr)
      
C     crack center (slightly offset from atom)
      xc = 0.5_dp
      yc = 0.3_dp
      
C     atomistic stuff
      natomisticsteps = 20
      dt = 0.02_dp
      dtdd = natomisticsteps*dt
      forcetol = 5.0e-3_dp
      
C     file
      write (gammasuffix,'(I0)') nint(normaldamping%gamma*100.0_dp)
      write (dtsuffix,'(I0)') nint(dt*1000)
      write (stepssuffix,'(I0)') natomisticsteps
      filename = trim(misc%simname)//'_gamma_'//trim(gammasuffix)
     &                             //'_dt_'//trim(dtsuffix)
     &                             //'_steps_'//trim(stepssuffix)
      open(newunit=iunit,file=filename)

      
C     apply K, equilibrate, dump
      do i = 0, nstepsK
          KIapply = KIstart + i*KIincr
          write(*,*) 'Current KI', KIapply
          
          call applyKDispIso(KIapply,KII,mu,nu,xc,yc,'all')
          call equilibrateCADD(natomisticsteps,dt,dtdd,normaldamping,
     &                    counter,forcetol,natomisticstepstot)
          write(iunit,*) KIapply, counter, natomisticstepstot
          call updateMiscIncrementCurr(1)
          call writeDump_ptr()
      end do
      
      call writeRestart_ptr()
      
      close(iunit)

      contains
************************************************************************
      subroutine equilibrateCADD(natomisticsteps,dt,dtdd,damp,counter,
     &                           forcetol,natomisticstepstot)

C     input variables
      integer :: natomisticsteps
      real(dp) :: dt, dtdd
      type(dampingdata) :: damp
      real(dp) :: forcetol
      
C     output variables
      integer :: natomisticstepstot

C     local variables
      real(dp) :: forcenorm
      integer :: counter
      logical :: detected
      integer :: mnumfe

C     we will check forces on atoms to see if we're done 
      forcenorm = huge(0.0_dp)
      counter = 0
      mnumfe = 1
      
      do while (forcenorm > forcetol)
C         write(*,*) 'Starting iteration'
          counter = counter + 1
          
          call loopVerlet(natomisticsteps,dt,'all',damp) ! step 1 (see Algorithm.txt)
          
          call detectAndPassDislocations(detected) ! step 2
          if (detected) then
              call updateMiscIncrementCurr(10)
              write(*,*) 'After passing'
              write(*,*) 'increment', misc%incrementcurr
              call writeDump_ptr()
          end if
          
          if (escapeddisl(mnumfe)%nescapeddisl > 0) then
              write(*,*) 'ESCAPE'
          end if
          
          call solveAll_ptr() ! step 3
          call updatePad() ! step 4          
          call runDDStep(dtdd) ! step 5 
          
          forcenorm = maxval(sum(abs(nodes%potforces),1)) ! 1-norm
          write(*,*) 'forcenorm', forcenorm
          
      end do
      natomisticstepstot = counter*natomisticsteps
      
      end subroutine equilibrateCADD
************************************************************************
      end program