      module mod_nodes

C     Purpose: Has routines for reading, initializing, and writing
C     node-related data (positions, velocities, displacements, atom types,
C     etc.)
C     
C     Possible extensions: ?

      use mod_types, only: dp
      use mod_utils, only: writeRealMat, writeIntMat
      implicit none
      
      private
      public :: nodes, initNodeData, writeNodeData,
     & readNodeData, processNodeData
      
      type nodedata
C     read-in
      integer :: nnodes
      integer, allocatable :: types(:,:)
      real(dp), allocatable :: posn(:,:)
C     processed
      integer :: natoms
      integer, allocatable :: atomlist(:)
      integer :: nfenodes
      integer, allocatable :: fenodelist(:)
      integer :: npadatoms
      integer, allocatable :: padatomlist(:)
      real(dp), allocatable :: potforces(:,:)
      real(dp), allocatable :: dampforces(:,:)
      end type
      
C     module variables
      type(nodedata) :: nodes
      
      contains
      
************************************************************************
      subroutine initNodeData(nodefile)
      
C     Subroutine: initNodeData
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
      
C     Subroutine: readNodeData

C     Inputs: nodefile - name of file containing node information
C     (e.g. example_fortran_nodes)

C     Outputs: None

C     Purpose: Read data into "nodes" from file

      implicit none
      
C     input variables
      character(len=*) :: nodefile
      
C     local variables
      integer :: iunit, nrow, j, k
      
      open(newunit=iunit,file=nodefile)
      
C     read posn
      read(iunit,*) nodes%nnodes, nrow
      allocate(nodes%posn(nrow,nodes%nnodes))
      read(iunit,*) ((nodes%posn(k,j),k=1,nrow),j=1,nodes%nnodes)
C     read types
      read(iunit,*) nodes%nnodes, nrow
      allocate(nodes%types(nrow,nodes%nnodes))
      read(iunit,*) ((nodes%types(k,j),k=1,nrow),j=1,nodes%nnodes)
      
      close(iunit)
      
      end subroutine readNodeData
************************************************************************
      subroutine processNodeData()

C     Subroutine: processNodeData

C     Inputs: None

C     Outputs: None

C     Purpose: Create lists of atoms, FE nodes, etc., initialize other
C     arrays.

      implicit none
      
C     local variables
      integer :: i
      integer :: curr
      integer :: counteratom, counterpadatom, counterfe
      integer, allocatable :: atomlisttemp(:)
      integer, allocatable :: fenodelisttemp(:)
      integer, allocatable :: padatomlisttemp(:)
      
C     allocate, initialize counters/temp arrays
      allocate(atomlisttemp(nodes%nnodes))
      allocate(fenodelisttemp(nodes%nnodes))
      allocate(padatomlisttemp(nodes%nnodes))
      counteratom = 0
      counterpadatom = 0
      counterfe = 0
      
C     loop through all nodes
      do i = 1, nodes%nnodes
          curr = nodes%types(2,i)
          if ((curr==-1).or.(curr==1).or.(curr==2)) then ! all atoms
              counteratom = counteratom + 1
              atomlisttemp(counteratom) = i
              if (curr==-1) then ! just pad atoms
                  counterpadatom = counterpadatom + 1
                  padatomlisttemp(counterpadatom) = i
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
      allocate(nodes%atomlist(nodes%natoms))
      allocate(nodes%padatomlist(nodes%npadatoms))
      allocate(nodes%fenodelist(nodes%nfenodes))
      
C     build lists from temp lists
      nodes%atomlist = atomlisttemp(1:nodes%natoms)
      nodes%padatomlist = padatomlisttemp(1:nodes%npadatoms)
      nodes%fenodelist = fenodelisttemp(1:nodes%nfenodes)
      
C     deallocate temp lists
      deallocate(atomlisttemp)
      deallocate(padatomlisttemp)
      deallocate(fenodelisttemp)
      
C     allocate, initialize other arrays
      allocate(nodes%potforces(2,nodes%natoms))
      nodes%potforces = 0.0_dp
      allocate(nodes%dampforces(2,nodes%natoms))
      nodes%dampforces = 0.0_dp
      
      end subroutine processNodeData
************************************************************************
      subroutine writeNodeData(nodefile)
      
C     Subroutine: writeNodeData

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
      
      call writeRealMat(nodes%posn,iunit,.true.)
      call writeIntMat(nodes%types,iunit,.true.)
      
      close(iunit)
      
      end subroutine writeNodeData
************************************************************************
      
      end module