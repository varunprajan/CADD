      module mod_io
      
C     Purpose: Contains routines for reading/initializing and writing data
C     from/to files. Can be used to write "dump" or "restart" files.

C     Possible extensions: ?
      
      use mod_types, only: dp
      use mod_nodes, only: nodes, initNodeData, writeNodeData
      use mod_misc, only: misc, initMiscData, writeMiscData
      use mod_damping, only: initDampingData, writeDampingData
      use mod_fe_elements, only: initFEElementData, writeFEElementData
      use mod_materials, only: initMaterialData, writeMaterialData
      use mod_potentials, only: initPotentialData, writePotentialData
      use mod_poten_main, only: assignPotStyle
      use mod_interactions, only: initInteractionData,
     &                            writeInteractionData
      use mod_neighbors, only: initNeighborData, writeNeighborData
      use mod_disl, only: initDislData, writeDislData, disl,
     &                    initImageDislData, writeImageDislData
      use mod_disl_detect_pass, only: initDetectionData,
     &                                writeDetectionData
      use mod_groups, only: initGroupData, writeGroupData
      use mod_math, only: isMultiple
      
      implicit none
      
      private :: suffixnodes, suffixmaterials, suffixpotentials,
     &           suffixmisc, suffixdisl, suffiximagedisl,
     &           suffixinteractions, suffixneighbors, suffixdamping,
     &           suffixfeelements
      public :: writeDump, initAll, writeAll, initAtomistic
      
C     module variables
      character(len=6) :: suffixnodes = '_nodes'
      character(len=10) :: suffixmaterials = '_materials'
      character(len=11) :: suffixpotentials = '_potentials'
      character(len=5) :: suffixmisc = '_misc'
      character(len=5) :: suffixdisl = '_disl'
      character(len=10) :: suffiximagedisl = '_imagedisl'
      character(len=13) :: suffixinteractions = '_interactions'
      character(len=10) :: suffixneighbors = '_neighbors'
      character(len=8) :: suffixdamping = '_damping'
      character(len=11) :: suffixfeelements = '_feelements'
      character(len=10) :: suffixdetection = '_detection'
      character(len=7) :: suffixgroups = '_groups'
      
      contains
************************************************************************
      subroutine readSimname()
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file='main')
      read(iunit,*) misc%simname
      close(iunit)
      
      end subroutine readSimname
************************************************************************
      subroutine writeOutput()
      
      if (isMultiple(misc%incrementcurr,
     &               misc%dumpincrement,.false.)) then
          call writeDump()
      end if
      if (isMultiple(misc%incrementcurr,
     &               misc%restartincrement,.false.)) then
          call writeAll()
      end if
      
      end subroutine writeOutput
************************************************************************        
      subroutine writeDump()
      
C     TODO: What else is required for dump files?
C     TODO: Need to add back in tilde (DD) displacements to FE nodes
      
C     local variables
      integer :: iunit, j, i, nrow
      character(len=15) :: incrementstr
      character(len=:), allocatable :: pref
      
      write (incrementstr,'(I0)') misc%incrementcurr
      pref = trim(misc%simname)
      open(newunit=iunit,file=pref//'.'//trim(incrementstr)//'.dump')
      
C     nodes
      write(iunit,*) 'positions:'
      do i = 1,nodes%nnodes
          write(iunit,*) (nodes%posn(j,i), j=1,3)
      end do
      write(iunit,*) 'end'
      write(iunit,*) 'velocities:'
      do i = 1,nodes%nnodes
          write(iunit,*) (nodes%posn(j,i), j=6,7)
      end do
      write(iunit,*) 'end'
      write(iunit,*) 'types:'
      nrow = size(nodes%types,1)
      do i = 1,nodes%nnodes
          write(iunit,*) (nodes%types(j,i), j=1,nrow)
      end do
      write(iunit,*) 'end'
      
C     TODO: need to fix dislocation print-out!
CC     dislocations
CC     presumably, we don't need the material types?
C     write(iunit,*) 'Dislocation positions:'
C     nrow = size(disl%posn,2)
C     do i = 1,disl%ndisl
C        write(iunit,*) (disl%posn(j,i), j=1,nrow)
C     end do
C     write(iunit,*) 'end'
      
      close(iunit)
      
      end subroutine writeDump
************************************************************************        
      subroutine initAll()

C     local variables      
      character(len=:), allocatable :: pref

      call readSimname()
      pref = trim(misc%simname)
      call initNodeData(pref//suffixnodes)
      call initMaterialData(pref//suffixmaterials)
      call initPotentialData(pref//suffixpotentials)
      call initMiscData(pref//suffixmisc)
      call assignPotStyle()
      call initDampingData(pref//suffixdamping)
      call initGroupData(pref//suffixgroups)
      call initFEElementData(pref//suffixfeelements)
      call initDislData(pref//suffixdisl)
      call initImageDislData(pref//suffiximagedisl)
      call initInteractionData(pref//suffixinteractions)
      call initNeighborData(pref//suffixneighbors)
      call initDetectionData(pref//suffixdetection)
      
      end subroutine initAll
************************************************************************        
      subroutine initAtomistic()

C     local variables      
      character(len=:), allocatable :: pref

      call readSimname()
      pref = trim(misc%simname)
      call initNodeData(pref//suffixnodes)
      call initMaterialData(pref//suffixmaterials)
      call initPotentialData(pref//suffixpotentials)
      call initMiscData(pref//suffixmisc)
      call assignPotStyle()
      call initDampingData(pref//suffixdamping)
      call initGroupData(pref//suffixgroups)
      call initInteractionData(pref//suffixinteractions)
      call initNeighborData(pref//suffixneighbors)
      
      end subroutine initAtomistic
************************************************************************
      subroutine writeAll()
      
C     can be used to generate a set of "restart" files, which can
C     be read in (using initAll(filepref)) to run simulation
C     In other words, writeAll is a mirror image of initAll's read
C     statements

C     local variables
      character(len=15) :: incrementstr
      character(len=:), allocatable :: suff
      character(len=:), allocatable :: pref
      
      write (incrementstr,'(I0)') misc%incrementcurr
      suff = '.'//trim(incrementstr)//'.restart'
      pref = trim(misc%simname)
      
      call writeNodeData(pref//suffixnodes//suff)
      call writeMaterialData(pref//suffixmaterials//suff)
      call writePotentialData(pref//suffixpotentials//suff)
      call writeMiscData(pref//suffixmisc//suff)
      call writeDampingData(pref//suffixdamping//suff)
      call writeGroupDAta(pref//suffixgroups//suff)
      call writeFEElementData(pref//suffixfeelements//suff)
      call writeDislData(pref//suffixdisl//suff)
      call writeImageDislData(pref//suffixdisl//suff)
      call writeInteractionData(pref//suffixinteractions//suff)
      call writeNeighborData(pref//suffixneighbors//suff)
      call writeDetectionData(pref//suffixdetection//suff)
      
      end subroutine writeAll
************************************************************************  
      end module