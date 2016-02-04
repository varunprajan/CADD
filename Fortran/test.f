      program test

      use mod_types, only: dp
      use mod_neighbors, only: updateNeighbors, neighbors
      use mod_misc, only: misc
      use mod_utils, only: prettyPrintIntMat, prettyPrintRealMat
      use mod_integrate, only: velocityVerlet
      use mod_io, only: readAll, writeAll, writeDump
      use mod_init, only: initAll
      implicit none
      
      character(len=:), allocatable :: simname
      integer :: i
      
      call initSimulation()
      
C     loop
C     do i = 1,misc%increments
      do i = 1,1
          write(*,*) i
          misc%incrementcurr = i
          
          write(*,*) 'YayStart'
          
C         reneighbor
          write(*,*) neighbors%every
          if (mod(i,neighbors%every)==0) then
              write(*,*) 'Inside'
              call updateNeighbors()
          end if
          
          write(*,*) 'YayNeighbors'
          
C         velocity verlet (main MD routine)
          call velocityVerlet(misc%timestep)
          
          write(*,*) 'YayVV'
          
C         dump
          if (mod(i,misc%dumpincrement)==0) then
              call writeDump()
          end if
          
          write(*,*) 'YayDump'
          
C         restart
          if (mod(i,misc%restartincrement)==0) then
              call writeAll()
          end if
          
          write(*,*) 'YayRestart'
      end do    

C     call getForcesAll()
C     call prettyPrintRealMat(nodes%realatomforces,'forces')
      
      
      end program