      module mod_nodes

C     Purpose: Has routines for reading, initializing, and writing
C     node-related data (positions, velocities, displacements, atom types,
C     etc.)
C     
C     Possible extensions: ?

      use mod_types, only: dp
      use mod_utils, only: writeMatTransposeSize, readMatTransposeSize
      implicit none
      
      private
      public :: nodes, initNodeData, writeNodeData,
     & readNodeData, processNodeData, getXYAtomBounds, getXYBounds,
     & getXYAtomBoundsDef, getXYAtomBoundsUndef
      
      type nodedata
C     read-in
      integer, allocatable :: types(:,:)
      real(dp), allocatable :: posn(:,:)
C     processed
      integer :: nnodes
      integer :: natoms
      integer, allocatable :: atomlist(:)
      integer :: nfenodes
      integer, allocatable :: fenodelist(:)
      integer :: npadatoms
      integer, allocatable :: padatomlist(:)
      integer :: nrealatoms
      integer, allocatable :: realatomlist(:)
      real(dp), allocatable :: potforces(:,:)
      real(dp), allocatable :: dampforces(:,:)
      end type
      
C     module variables
      type(nodedata) :: nodes
      
      contains
************************************************************************
      subroutine initNodeData(nodefile)
C     
C     Inputs: nodefile - name of file containing node information
C     (e.g. example_fortran_nodes)

C     Outputs: None
C     
C     Purpose: Read, initialize data in "nodes" structure

      implicit none
      
C     input variables
      character(len=*) :: nodefile
      
      call readNodeData(nodefile)
      call processNodeData()
      
      end subroutine initNodeData
************************************************************************
      subroutine readNodeData(nodefile)

C     Inputs: nodefile - name of file containing node information
C     (e.g. example_fortran_nodes)

C     Outputs: None

C     Purpose: Read data into "nodes" from file

      implicit none
      
C     input variables
      character(len=*) :: nodefile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=nodefile)
      
      call readMatTransposeSize(iunit,nodes%posn)
      call readMatTransposeSize(iunit,nodes%types)
      nodes%nnodes = size(nodes%posn,2)
      
      close(iunit)
      
      end subroutine readNodeData
************************************************************************
      subroutine processNodeData()

C     Inputs: None

C     Outputs: None

C     Purpose: Create lists of atoms, FE nodes, etc., initialize other
C     arrays.

      implicit none
      
C     local variables
      integer :: i
      integer :: curr
      integer :: counteratom, counterpadatom, counterfe, counterrealatom
      integer :: atomlisttemp(nodes%nnodes)
      integer :: fenodelisttemp(nodes%nnodes)
      integer :: padatomlisttemp(nodes%nnodes)
      integer :: realatomlisttemp(nodes%nnodes)
      
      counteratom = 0
      counterpadatom = 0
      counterrealatom = 0
      counterfe = 0
      
      do i = 1, nodes%nnodes
          curr = nodes%types(2,i)
          if ((curr==-1).or.(curr==1).or.(curr==2)) then ! all atoms
              counteratom = counteratom + 1
              atomlisttemp(counteratom) = i
              if (curr==-1) then ! just pad atoms
                  counterpadatom = counterpadatom + 1
                  padatomlisttemp(counterpadatom) = i
              else
                  counterrealatom = counterrealatom + 1
                  realatomlisttemp(counterrealatom) = i
              end if    
          end if
          if ((curr==0).or.(curr==2)) then ! all nodes
              counterfe = counterfe + 1
              fenodelisttemp(counterfe) = i
          end if 
      end do
      nodes%natoms = counteratom
      nodes%npadatoms = counterpadatom
      nodes%nfenodes = counterfe
      nodes%nrealatoms = counterrealatom
      
C     build lists from temp lists (allocation automatically taken care of)
      nodes%atomlist = atomlisttemp(1:nodes%natoms)
      nodes%padatomlist = padatomlisttemp(1:nodes%npadatoms)
      nodes%fenodelist = fenodelisttemp(1:nodes%nfenodes)
      nodes%realatomlist = realatomlisttemp(1:nodes%nrealatoms)
      
C     misc. arrays
      allocate(nodes%potforces(2,nodes%natoms))
      nodes%potforces = 0.0_dp
      allocate(nodes%dampforces(2,nodes%natoms))
      nodes%dampforces = 0.0_dp
      
      end subroutine processNodeData
************************************************************************
      subroutine getXYAtomBoundsDef(xmin,xmax,ymin,ymax)

C     Inputs: None

C     Outputs: xmin, xmax, ymin, ymax --- max x and y positions for all atoms

C     Purpose: Get max, min x- and y-positions for all atoms
C     (including pad atoms) in *deformed* state

      implicit none

C     output variables
      real(dp) :: xmin, xmax
      real(dp) :: ymin, ymax
      
      call getXYAtomBounds(.false.,xmin,xmax,ymin,ymax)
      
      end subroutine getXYAtomBoundsDef
************************************************************************
      subroutine getXYAtomBoundsUndef(xmin,xmax,ymin,ymax)

C     Inputs: None

C     Outputs: xmin, xmax, ymin, ymax --- max x and y positions for all atoms

C     Purpose: Get max, min x- and y-positions for all atoms
C     (including pad atoms) in *undeformed* state

      implicit none

C     output variables
      real(dp) :: xmin, xmax
      real(dp) :: ymin, ymax
      
      call getXYAtomBounds(.true.,xmin,xmax,ymin,ymax)
      
      end subroutine getXYAtomBoundsUndef
************************************************************************
      subroutine getXYAtomBounds(undeformed,xmin,xmax,ymin,ymax)

C     Inputs: undeformed --- flag indicating whether undeformed positions are to be used

C     Outputs: xmin, xmax, ymin, ymax --- max x and y positions for all atoms

C     Purpose: Get max, min x- and y-positions for all atoms
C     (including pad atoms)

      implicit none

C     input variables
      logical :: undeformed

C     output variables
      real(dp) :: xmin, xmax
      real(dp) :: ymin, ymax
      
C     local variables
      integer :: i, atom
      real(dp) :: xcurr, ycurr

C     get x, y bounds by looping over all atoms
      xmin = huge(xmin)
      xmax = -huge(xmax)
      ymin = huge(ymin)
      ymax = -huge(ymax)
      do i = 1, nodes%natoms
          atom = nodes%atomlist(i)
          xcurr = nodes%posn(1,atom)
          ycurr = nodes%posn(2,atom)
          if (undeformed) then
              xcurr = xcurr - nodes%posn(4,atom)
              ycurr = ycurr - nodes%posn(5,atom)
          end if    
          xmin = min(xmin,xcurr)
          xmax = max(xmax,xcurr)
          ymin = min(ymin,ycurr)
          ymax = max(ymax,ycurr) 
      end do
      
      end subroutine getXYAtomBounds
************************************************************************
      subroutine getXYBounds(xmin,xmax,ymin,ymax)

C     Inputs: None

C     Outputs: None

C     Purpose: Get max, min x- and y-positions for all nodes in *undeformed* state

      implicit none

C     output variables
      real(dp) :: xmin, xmax
      real(dp) :: ymin, ymax
      
C     local variables
      integer :: i
      real(dp) :: xdef, ydef, xdisp, ydisp, xcurr, ycurr

      xmin = huge(xmin)
      xmax = -huge(xmax)
      ymin = huge(ymin)
      ymax = -huge(ymax)
      do i = 1, nodes%nnodes
          xdef = nodes%posn(1,i)
          ydef = nodes%posn(2,i)
          xdisp = nodes%posn(4,i)
          ydisp = nodes%posn(5,i)
          xcurr = xdef - xdisp
          ycurr = ydef - ydisp
          xmin = min(xmin,xcurr)
          xmax = max(xmax,xcurr)
          ymin = min(ymin,ycurr)
          ymax = max(ymax,ycurr) 
      end do
      
      end subroutine getXYBounds
************************************************************************
      subroutine writeNodeData(nodefile)

C     Inputs: nodefile --- name of file containing node information
C     (e.g. example_fortran_nodes)

C     Outputs: None

C     Purpose: Write node data back to file (inverse of readNodeData),
C     for the purpose of writing restart files

      implicit none
      
C     input variables
      character(len=*) :: nodefile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=nodefile)
      
      call writeMatTransposeSize(iunit,nodes%posn)
      call writeMatTransposeSize(iunit,nodes%types)
      
      close(iunit)
      
      end subroutine writeNodeData
************************************************************************
      
      end module