      program atomistic_1
      
      use mod_types, only: dp
      use mod_io, only: writeOutput, initAtomistic
      use mod_integrate, only: loopVerlet, velocityVerletSubAll
      use mod_kdispfield, only: applyKDispIso
      use mod_materials, only: materials
      use mod_misc, only: misc
      use mod_damping, only: damping, getDampingForcesAll
      use mod_poten_main, only: getPotForcesAll_ptr
      use mod_neighbors, only: updateNeighbors
      use mod_math, only: isMultiple
      
      implicit none
      
      real(dp) :: KIstart, KIend, KIincr, KII
      real(dp) :: mu, nu
      real(dp) :: xc, yc
      integer :: i, j, itot
      integer :: gnum
      
      write(*,*) 'Yay'
      call initAtomistic()
      write(*,*) misc%simname
      
      KIstart = 0.0_dp
      KIend = 6.0_dp
      KIincr = 6.0_dp
      itot = nint((KIend-KIstart)/KIincr)
      mu = 10.633
      nu = materials(1)%nu
      xc = 0.25_dp
      yc = 0.25_dp
      KII = 0.0_dp
      
      call writeDump_ptr()
      write(*,*) misc%timestep
      call runNVE('all',5000)
      call writeDump_ptr()

C         write(*,*) 'i', i
C         call applyKDispIso(KIincr,KII,mu,nu,xc,yc,gnum)
C         call updateNeighbors(.false.)

      contains
************************************************************************
      subroutine runNVE(gname,increments)

      implicit none

C     local variables
      character(len=*) :: gname
      integer :: increments

C     initialize all forces
      call getPotForcesAll_ptr()
      call getDampingForcesAll(damping%gname,damping%gamma,damping%flag)
      
C     run verlet loop
      do i = 1, increments
          if (isMultiple(misc%incrementcurr,100)) then
              write(*,*) 'i', i
          end if    
          call updateNeighbors()
          call velocityVerletSubAll(misc%timestep,gname,
     &             damping%gname,damping%flag,damping%gamma)
          misc%incrementcurr = misc%incrementcurr + 1
          call writeDump_ptr()
          call writeRestart_ptr()
      end do
      
      end subroutine
************************************************************************      
      end program