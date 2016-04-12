      module mod_dump

C     Notes/TODO: Needs to be fixed for multi-material, and for 3D (currently, only 2D positions outputted)
C     For the latter, have to distinguish between atomistic simtypes and not.

      use mod_types, only: dp
      use mod_misc, only: misc
      use mod_nodes, only: nodes
      use mod_disl_try, only: disl, sources, obstacles
      use mod_fe_main_2d_assign, only: updateFENodalPosnAll_ptr
      use mod_fe_elements, only: feelements
      use mod_compute, only: compute, getCentroAtoms
      use mod_utils, only: writeMatTranspose, writeVec
      use mod_slip_sys, only: slipsys
      use mod_pad_atoms, only: updatePad
      
      implicit none
      
      private
      public :: getDumpFilename, writeDumpCADD,
     & writeDumpCADDNoDisl, writeDumpFE, writeDumpDD,
     & writeDumpAtomistic, writeDumpNodesElementsChunkSub,
     & writeDumpNodesDefElementsChunk, writeDumpNodesUndefElementsChunk,
     & writeDumpNodesUndefChunk,
     & writeDumpDDChunk, writeDumpNodesDefChunk, writeDumpNodesDef,
     & writeDumpNodesUndef, writeDumpNodesDisp, writeDumpNodesTypes,
     & writeDumpFEElements, writeDumpDisl, writeDumpSources,
     & writeDumpObstacles, writeDumpCompute, writeDumpCentro

      contains
************************************************************************
      function getDumpFilename() result(filename)

C     Inputs: None

C     Outputs: filename --- name of dump file

C     Purpose: Generate name of dump file from increment number:
C     e.g. if incrementcurr = 10, simname = test -> filename = 'test.10.dump'     
      
      implicit none
      
      character(len=15) :: incrementstr
      character(len=:), allocatable :: pref
      character(len=:), allocatable :: filename
      
      write (incrementstr,'(I0)') misc%incrementcurr
      pref = trim(misc%simname)
      filename = pref//'.'//trim(incrementstr)//'.dump'
      
      end function getDumpFilename
************************************************************************     
      subroutine writeDumpCADD()

C     Inputs: None

C     Outputs: None

C     Purpose: Write data for dump for a cadd simulation:
C     a) Nodal positions (deformed), displacements, and node types
C     b) FE elements (connectivity)
C     c) Dislocation, obstacle, and source positions and dislocation attributes (slip system, sign)
C     d) Compute data

      implicit none
      
C     local variables
      integer :: iunit
      character(len=:), allocatable :: filename
      
      filename = getDumpFilename()
      open(newunit=iunit,file=filename)
      
      call writeDumpNodesDefElementsChunk(iunit) ! deformed positions and displacements are most convenient
      call writeDumpDDChunk(iunit)
      call writeDumpCompute(iunit)
      
      close(iunit)
      
      end subroutine writeDumpCADD
************************************************************************
      subroutine writeDumpCADDNoDisl()

C     Inputs: None

C     Outputs: None

C     Purpose: Write data for dump for a cadd simulation with no dislocations/sources/obstacles:
C     a) Nodal positions (deformed), displacements, and node types
C     b) FE elements (connectivity)
C     c) Compute data
      
C     local variables
      integer :: iunit
      character(len=:), allocatable :: filename
      
      filename = getDumpFilename()
      open(newunit=iunit,file=filename)
      
      call writeDumpNodesDefElementsChunk(iunit) ! deformed positions and displacements are most convenient
      call writeDumpCompute(iunit)

      close(iunit)
      
      end subroutine writeDumpCADDNoDisl
************************************************************************ 
      subroutine writeDumpFE()

C     Inputs: None

C     Outputs: None

C     Purpose: Write data for dump for a pure fe simulation (i.e. with no dislocations/sources/obstacles):
C     a) Nodal positions (undeformed), displacements, and node types
C     b) FE elements (connectivity)
C     c) Compute data

      implicit none

C     local variables
      integer :: iunit
      character(len=:), allocatable :: filename
      
      filename = getDumpFilename()
      open(newunit=iunit,file=filename)
      
      call writeDumpNodesUndefElementsChunk(iunit) ! undeformed positions and displacements are most convenient
      call writeDumpCompute(iunit)
      
      close(iunit)
      
      end subroutine writeDumpFE
************************************************************************
      subroutine writeDumpDD()

C     Inputs: None

C     Outputs: None

C     Purpose: Write data for dump for a dd simulation (i.e. with dislocations/sources/obstacles):
C     a) Nodal positions (undeformed), displacements, and node types
C     b) FE elements (connectivity)
C     c) Dislocation, obstacle, and source positions and dislocation attributes (slip system, sign)
C     d) Compute data

      implicit none

C     local variables
      integer :: iunit
      character(len=:), allocatable :: filename

      filename = getDumpFilename()
      open(newunit=iunit,file=filename)
      
      call writeDumpNodesUndefElementsChunk(iunit) ! undeformed positions and displacements are most convenient
      call writeDumpDDChunk(iunit)
      call writeDumpCompute(iunit)
      
      close(iunit)
      
      end subroutine writeDumpDD
************************************************************************
      subroutine writeDumpAtomistic()

C     Inputs: None

C     Outputs: None

C     Purpose: Write data for dump for a pure atomistic simulation (i.e. no fe/dd):
C     a) Nodal positions (deformed) and types
C     b) Compute data

      implicit none

C     local variables
      integer :: iunit
      character(len=:), allocatable :: filename

      filename = getDumpFilename()
      open(newunit=iunit,file=filename)
      
      call writeDumpNodesDefChunk(iunit) ! deformed positions are most convenient
      call writeDumpCompute(iunit)
      
      close(iunit)
      
      end subroutine writeDumpAtomistic
************************************************************************
      subroutine writeDumpNodesDefElementsChunk(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write a chunk of data for dump
C     a) Deformed nodal positions, displacements, and types
C     b) FE element connectivity

C     Notes: Note that all FE nodal positions must be updated first by adding
C     hat and tilde fields back in *for all nodes* (not just boundary nodes);
C     this can be very costly for a DD or CADD simulation, so frequent dumps should be avoided.
      
      implicit none
      
C     input variables
      integer :: iunit
      
      call writeDumpNodesElementsChunkSub(iunit,.true.)
      
      end subroutine writeDumpNodesDefElementsChunk
************************************************************************
      subroutine writeDumpNodesUndefElementsChunk(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write a chunk of data for dump
C     a) Deformed nodal positions, displacements, and types
C     b) FE element connectivity

C     Notes: Note that all FE nodal positions must be updated first by adding
C     hat and tilde fields back in *for all nodes* (not just boundary nodes);
C     this can be very costly for a DD or CADD simulation, so frequent dumps should be avoided.
      
      implicit none
      
C     input variables
      integer :: iunit
      
      call writeDumpNodesElementsChunkSub(iunit,.false.)
      
      end subroutine writeDumpNodesUndefElementsChunk
************************************************************************
      subroutine writeDumpNodesElementsChunkSub(iunit,defoption)

C     Inputs: iunit --- integer file specifier
C             defoption --- true if deformed positions are to be dumped; false if undeformed positions

C     Outputs: None

C     Purpose: Write a chunk of data for dump
C     a) Nodal positions (deformed or undeformed), displacements, and types
C     b) FE element connectivity

C     Notes: Note that all FE nodal positions must be updated first by adding
C     hat and tilde fields back in *for all nodes* (not just boundary nodes);
C     this can be very costly for a DD or CADD simulation, so frequent dumps should be avoided.
      
      implicit none
      
C     input variables
      integer :: iunit
      logical :: defoption
      
      call updateFENodalPosnAll_ptr()
      call updatePad()
      if (defoption) then
          call writeDumpNodesDefChunk(iunit)
      else
          call writeDumpNodesUndefChunk(iunit)
      end if
      call writeDumpNodesDisp(iunit)
      call writeDumpFEElements(iunit)
      
      end subroutine writeDumpNodesElementsChunkSub
************************************************************************
      subroutine writeDumpDDChunk(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write a chunk of data for dump
C     a) Dislocation positions and attributes; obstacle and source positions
      
      implicit none
      
C     input variables
      integer :: iunit
      
      call writeDumpDisl(iunit)
      call writeDumpSources(iunit)
      call writeDumpObstacles(iunit)
      call writeDumpSlipSys(iunit)
      
      end subroutine writeDumpDDChunk
************************************************************************
      subroutine writeDumpNodesDefChunk(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write a chunk of data for dump
C     a) Deformed nodal positions and types

      implicit none

C     input variables
      integer :: iunit
      
      call writeDumpNodesTypes(iunit)
      call writeDumpNodesDef(iunit)
      
      end subroutine writeDumpNodesDefChunk
************************************************************************
      subroutine writeDumpNodesUndefChunk(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write a chunk of data for dump
C     a) Uneformed nodal positions and types

      implicit none

C     input variables
      integer :: iunit
      
      call writeDumpNodesTypes(iunit)
      call writeDumpNodesUndef(iunit)
      
      end subroutine writeDumpNodesUndefChunk
************************************************************************
      subroutine writeDumpNodesDef(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write deformed node positions for dump
      
      implicit none
      
C     input variables
      integer :: iunit
      
C     local variables
      integer :: i, j
      
      write(iunit,*) 'deformed_positions:'
      do i = 1, nodes%nnodes
          write(iunit,*) (nodes%posn(j,i), j=1,3)
      end do
      write(iunit,*) 'end'
      
      end subroutine writeDumpNodesDef
************************************************************************
      subroutine writeDumpNodesUndef(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write undeformed node positions for dump
      
      implicit none
      
C     input variables
      integer :: iunit
      
C     local variables
      integer :: i, j
      real(dp) :: disp(3)
      
      disp = 0.0_dp
      
      write(iunit,*) 'undeformed_positions:'
      do i = 1, nodes%nnodes
          disp(1:2) = nodes%posn(4:5,i)
          write(iunit,*) (nodes%posn(j,i)-disp(j), j=1,2)
      end do
      write(iunit,*) 'end'
      
      end subroutine writeDumpNodesUndef
************************************************************************
      subroutine writeDumpNodesDisp(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write nodal displacements for dump

C     input variables
      integer :: iunit

C     local variables
      integer :: i, j
      
      write(iunit,*) 'displacements:'
      do i = 1, nodes%nnodes
          write(iunit,*) (nodes%posn(j,i), j=4,5)
      end do
      write(iunit,*) 'end'
      
      end subroutine writeDumpNodesDisp
************************************************************************
      subroutine writeDumpNodesTypes(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write nodal types for dump

      implicit none

C     input variables
      integer :: iunit

      write(iunit,*) 'types:'
      call writeMatTranspose(iunit,nodes%types) 
      write(iunit,*) 'end'
      
      end subroutine writeDumpNodesTypes
************************************************************************
      subroutine writeDumpFEElements(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write finite element data for dump
      
C     Notes/TODO: Doesn't really work for multimaterial...need separate arrays?

      implicit none

C     input variables
      integer :: iunit
      
C     local variables
      integer :: i

      write(iunit,*) 'fe_elements:'
      do i = 1, size(feelements)
          call writeMatTranspose(iunit,feelements(i)%connect)
      end do
      write(iunit,*) 'end'
      
      end subroutine writeDumpFEElements
************************************************************************
      subroutine writeDumpDisl(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write dislocation data for dump (positions, attributes (slip system, sign))
      
      implicit none
      
C     input variables
      integer :: iunit
      
C     local variables
      integer :: i, j, k
      
      write(iunit,*) 'dislocation_positions:'
      do k = 1, size(disl) ! material number is not important
          do i = 1, disl(k)%ndisl
              if (disl(k)%list(i)%active) then
                  write(iunit,*) (disl(k)%list(i)%posn(j), j=1,2)
              end if    
          end do
      end do
      write(iunit,*) 'end'    
      
      write(iunit,*) 'dislocation_attributes:'
      do k = 1, size(disl) ! material number is not important
          do i = 1, disl(k)%ndisl
              if (disl(k)%list(i)%active) then
                  write(iunit,*) disl(k)%list(i)%slipsys,
     &                           disl(k)%list(i)%sgn
              end if 
          end do
      end do
      write(iunit,*) 'end'
      
      end subroutine writeDumpDisl
************************************************************************
      subroutine writeDumpSources(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write source positions for dump

      implicit none
      
C     input variables
      integer :: iunit
      
C     local variables
      integer :: i, j, k
      
      write(iunit,*) 'source_positions:'
      do k = 1, size(sources) ! material number is not important
          do i = 1, size(sources(k)%list)
              write(iunit,*) (sources(k)%list(i)%posn(j), j=1,2)    
          end do
      end do
      write(iunit,*) 'end'      
      
      end subroutine writeDumpSources
************************************************************************
      subroutine writeDumpObstacles(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write obstacle positions for dump

      implicit none

C     input variables
      integer :: iunit
      
C     local variables
      integer :: i, j, k
      
      write(iunit,*) 'obstacle_positions:'
      do k = 1, size(obstacles) ! material number is not important
          do i = 1, size(obstacles(k)%list)
              write(iunit,*) (obstacles(k)%list(i)%posn(j), j=1,2)    
          end do
      end do
      write(iunit,*) 'end'
      
      end subroutine writeDumpObstacles
************************************************************************
      subroutine writeDumpSlipSys(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write slip system angles      
      
      implicit none

C     input variables
      integer :: iunit
      
C     local variables
      integer :: i

      write(iunit,*) 'slipsys_angles:'
      do i = 1, size(slipsys) ! material number is not important
          call writeVec(iunit,slipsys(i)%theta)
      end do
      write(iunit,*) 'end'
      
      end subroutine writeDumpSlipSys
************************************************************************
      subroutine writeDumpCompute(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: If the relevant compute is active, call it, and
C     then write compute data to file
      
      implicit none
      
C     input variables
      integer :: iunit
      
      if (compute%centro%active) then
          call getCentroAtoms()
          call writeDumpCentro(iunit)
      end if
      
      end subroutine writeDumpCompute
************************************************************************
      subroutine writeDumpCentro(iunit)

C     Inputs: iunit --- integer file specifier

C     Outputs: None

C     Purpose: Write centrosymmetry data for each atom to file.
C     Assumes getCentroAtoms has already been called.
      
      implicit none
      
C     input variables
      integer :: iunit
      
      write(iunit,*) 'centro:'
      call writeVec(iunit,compute%centro%res)
      write(iunit,*) 'end'
      
      end subroutine writeDumpCentro
************************************************************************
      end module mod_dump