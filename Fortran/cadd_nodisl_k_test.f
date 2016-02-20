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
      
      implicit none
      
      integer :: iunit
      integer :: i
      integer :: mnumfe, mnum
      real(dp) :: KIstart, KIincr, KIend, KIapply, KIcurr, KII
      real(dp) :: xc, yc
      real(dp) :: mu, nu
      integer :: nstepsK, natomisticsteps, natomisticstepstot
      real(dp) :: dt, disptol
      character(len=15) :: gammasuffix, dtsuffix, stepssuffix
      character(len=:), allocatable :: filename
      
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
      disptol = 1.0e-5_dp
      
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
          if (i == 0) then
              KIapply = KIstart
          else
              KIapply = KIincr
          end if
          
          KIcurr = KIstart + i*KIincr
          write(*,*) 'Current KI', KIcurr
          
          call applyKDispIso(KIapply,KII,mu,nu,xc,yc,'all')
          call equilibrateCADDNoDisl(natomisticsteps,dt,normaldamping,
     &                               disptol,natomisticstepstot)
          write(iunit,*) KIcurr, natomisticstepstot
          call updateMiscIncrementCurr(1)
          call writeDump_ptr()
      end do
      
      close(iunit)

      contains
************************************************************************
      subroutine equilibrateCADDNoDisl(natomisticsteps,dt,damp,
     &                                 disptol,natomisticstepstot)

C     input variables
      integer :: natomisticsteps
      real(dp) :: dt
      type(dampingdata) :: damp
      real(dp) :: disptol 
      
C     output variables
      integer :: natomisticstepstot

C     local variables
      real(dp), allocatable :: allposnold(:,:)
      real(dp), allocatable :: allposnnew(:,:)
      real(dp) :: dispnorm
      integer :: counter

C     we will check nodal positions to see if we're done
      allposnold = nodes%posn(1:2,:)    
      dispnorm = huge(0.0_dp)
      counter = 0
      
      do while (dispnorm > disptol)
          write(*,*) 'Starting iteration'
          call loopVerlet(natomisticsteps,dt,'all',damp) ! step 1 (see Algorithm.txt)
          counter = counter + 1
          
          call solveAll_ptr() ! step 3
          call updatePad() ! step 4
          
C         need to update positions to compare
          call updateFENodalPosnAll_ptr()
          allposnnew = nodes%posn(1:2,:)
          dispnorm = getDispNorm(allposnold,allposnnew,nodes%types)
          allposnold = allposnnew
      end do
      natomisticstepstot = counter*natomisticsteps
      
      end subroutine equilibrateCADDNoDisl
************************************************************************
      function getDispNorm(posnold,posnnew,types) result(dispnormmax)
      
      implicit none
      
C     input variables
      real(dp) :: posnold(:,:)
      real(dp) :: posnnew(:,:)
      integer :: types(:,:)
      
C     output variables
      real(dp) :: dispnormmax
      
C     local variables
      integer :: i
      integer :: nodetype
      real(dp) :: dispnorm, dispnormatoms, dispnormfenodes
      
      dispnormatoms = 0.0_dp
      dispnormfenodes = 0.0_dp
      
      do i = 1, size(types,2)
          nodetype = types(2,i)
          dispnorm = sum(abs(posnnew(:,i) - posnold(:,i))) ! infinity norm
          if (nodetype == 0) then ! fe node
              if (dispnorm > dispnormfenodes) then
                  dispnormfenodes = dispnorm
              end if
          else if ((nodetype == 1).or.(nodetype == 2)) then
              if (dispnorm > dispnormatoms) then
                  dispnormatoms = dispnorm
              end if
          end if
      end do
      
      write(*,*) 'Displacement inf. norm atoms', dispnormatoms
      write(*,*) 'Displacement inf. norm fe nodes', dispnormfenodes
      
      dispnormmax = max(dispnormatoms,dispnormfenodes)
      
      end function getDispNorm
************************************************************************
      end program