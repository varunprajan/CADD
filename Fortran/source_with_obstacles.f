      program source_with_obstacle
      
      use mod_types, only: dp
      use mod_io, only: writeOutput, initSimulation, writeRestart_ptr,
     &                  writeDump_ptr
      use mod_misc, only: updateMiscIncrementCurr, misc
      use mod_math, only: isMultiple, linspace
      use mod_materials, only: materials
      use mod_nodes, only: nodes
      use mod_groups, only: groups, allgroupname, getGroupNum
      use mod_fe_elements, only: fematerials
      use mod_fe_main_2d, only: getFEStressAtPoint
      use mod_fe_main_2d_assign, only: solveAll_ptr
      use mod_disl_try, only: disl, dislt, addDislocation,
     &  sortedplanedata
      use mod_disl_fields2, only: getPKTildeStressAll,
     &  getTildeStressAtPointAll, getGhostStressAtPointAll,
     &  getLatentStressAtPointAll, getRealStressAtPointAll
      use mod_dd_main, only: runDDStep, getResolvedStressOnSource
      use mod_dd_integrate, only: getResolvedStressOnDisl
      use mod_slip_sys ,only: resolveStress
      use mod_disl_misc, only: dislmisc
      
      implicit none
      
      real(dp) :: tauinit, taufinal, tauincr, taudisl
      real(dp), allocatable :: tauvec(:)
      real(dp) :: mu, gammaincr
      integer :: i
      integer :: ntau, nsteps
      real(dp) :: dt, dt0
      integer :: dislnum
      integer :: mnum, mnumfe
      real(dp) :: vmax
      character(len=15) :: facsuffix, vmaxsuffix
      character(len=:), allocatable :: filename
      integer :: iunit
      integer :: fac
      
      call initSimulation('source_with_obstacles','dd')
      call writeDump_ptr()
 
C     material properties
      mnumfe = 1
      mnum = fematerials%list(mnumfe)
      mu = materials(mnum)%mu
      ! materials(mnumfe)%dislvmax = 20.0_dp ! ad hoc cutoff
      vmax = materials(mnum)%dislvmax
      
C     timestep (nanoseconds)
      dt0 = 0.5_dp
      fac = 50
      dt = dt0/real(fac,dp)
      nsteps = 50*fac
      
C     file
      write (facsuffix,'(I0)') fac
      write (vmaxsuffix,'(I0)') nint(vmax)
      if (dislmisc%gradientcorrection) then
          filename = trim(misc%simname)//'_'//trim(facsuffix)//
     &               '_vmax_'//trim(vmaxsuffix)//'_corr'
      else
          filename = trim(misc%simname)//'_'//trim(facsuffix)//
     &               '_vmax_'//trim(vmaxsuffix)//'_uncorr'
      end if
      open(newunit=iunit,file=filename)
      
C     loop over shear stresses
      tauinit = 0.0_dp
      taufinal = 400.0_dp
      ntau = 150
      allocate(tauvec(ntau))
      tauvec = linspace(tauinit,taufinal,ntau)
      do i = 1, ntau - 1
          write(*,*) 'i', i
C         stress, strain
          tauincr = tauvec(i+1) - tauvec(i)
          gammaincr = tauincr/mu
          call applyShearDisplacements(allgroupname,gammaincr)
          
C         dd routine
          call runDDCycle(nsteps,dt)
          
C         report stress at pinned dislocation
          dislnum = getMaxDislX(mnumfe)
          if (dislnum /= 0) then
              taudisl = getResolvedStressOnDisl(mnumfe,dislnum)
              write(iunit,*) tauvec(i+1), taudisl, disl(mnumfe)%ndisl ! applied stress vs. stress on leading dislocation vs. number of dislocations
          end if
          
C         dump
          call updateMiscIncrementCurr(1)
          call writeDump_ptr()
      end do
      
      close(iunit)
          
      contains
************************************************************************
      subroutine runDDCycle(nsteps,dt)
      
C     input variables
      integer :: nsteps
      real(dp) :: dt
      
C     local variables
      integer :: i
      
      mnumfe = 1
      do i = 1, nsteps
          ! write(*,*) 'i', i
          call solveAll_ptr()
          call runDDStep(dt)
      end do
      
      end subroutine runDDCycle          
************************************************************************
      function getMaxDislX(mnumfe) result(dislnum)
      
C     input variables
      integer :: mnumfe
      
C     output variables
      integer :: dislnum
      
C     local variables 
      integer :: i
      real(dp) :: posn(2), posncurr(2)
      
      posn(1) = -huge(dp)
      dislnum = 0
      do i = 1, disl(mnumfe)%ndisl
          if (disl(mnumfe)%list(i)%active) then
              posncurr = disl(mnumfe)%list(i)%posn
              if (posncurr(1) > posn(1)) then
                  dislnum = i
                  posn = posncurr
              end if
          end if
      end do
      
      end function getMaxDislX
************************************************************************
      subroutine applyShearDisplacements(gname,shearstrainincr)
      
      implicit none
      
C     input variables
      character(len=*) :: gname
      real(dp) :: shearstrainincr
      
C     local variables
      integer :: i
      integer :: gnum
      real(dp) :: halfshearstrainincr
      real(dp) :: xpos, ypos
      real(dp) :: xdispold, ydispold
      real(dp) :: xdisp, ydisp
      
      halfshearstrainincr = 0.5_dp*shearstrainincr
      gnum = getGroupNum(gname)
      do i = 1, nodes%nnodes
C     only do operation for atoms in group
      if (groups(gnum)%maskall(i)) then
          xpos = nodes%posn(1,i)
          ypos = nodes%posn(2,i)
          xdispold = nodes%posn(4,i)
          ydispold = nodes%posn(5,i)
          xdisp = ypos*halfshearstrainincr
          ydisp = xpos*halfshearstrainincr
          nodes%posn(1,i) = xpos + xdisp
          nodes%posn(2,i) = ypos + ydisp
          nodes%posn(4,i) = xdispold + xdisp
          nodes%posn(5,i) = ydispold + ydisp
      end if
      end do
      
      end subroutine applyShearDisplacements
************************************************************************
      end program