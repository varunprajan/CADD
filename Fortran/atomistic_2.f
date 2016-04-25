      program atomistic_2
      
      use mod_types, only: dp
      use mod_io, only: writeOutput, initSimulation, writeRestart_ptr,
     &                  writeDump_ptr
      use mod_integrate, only: loopVerlet, velocityVerletSubAll
      use mod_misc, only: updateMiscIncrementCurr
      use mod_damping, only: normaldamping, getDampingForcesAll
      use mod_poten_assign, only: getPotForcesAll_ptr
      use mod_neighbors, only: updateNeighborsCheck, 
     &   updateNeighborsNoCheck, updateNeighIncrementCurr
      use mod_math, only: isMultiple
      
      implicit none
      
      call initSimulation('atomistic_2','atomistic')
      call writeDump_ptr()
      call writeRestart_ptr()
C     K displacement field has already been applied
      call runNVE('all',5000,0.02_dp)
      call writeRestart_ptr()

      contains
************************************************************************
      subroutine runNVE(gname,increments,timestep)

      implicit none

C     input variables
      character(len=*) :: gname
      integer :: increments
      real(dp) :: timestep
      
C     local variables
      integer :: i

C     initialize neighbors, all forces
      call updateNeighborsNoCheck()
      call getPotForcesAll_ptr()
      call getDampingForcesAll(normaldamping)
      
C     run verlet loop
      do i = 1, increments
          if (isMultiple(i-1,100)) then
              write(*,*) 'i', i
          end if
          call updateNeighIncrementCurr(1)
          call velocityVerletSubAll(timestep,gname,normaldamping)
          call updateNeighborsCheck()
          
          call updateMiscIncrementCurr(1)
          call writeOutput()
      end do
      
      end subroutine
************************************************************************      
      end program