      module mod_interactions

C     Purpose: Has routines for read/writing/generating interactions,
C     which describe which potential should be used for the bonds
C     between atoms i and j. Currently, only interactions for two-body
C     potentials have been implemented.

C     Possible extensions: Three-body potentials. This is reasonably
C     straightforward; instead of using mat (2d array), one could use
C     mat3d (3d array) for three-body interactions, etc.

      use mod_types, only: dp
      use mod_utils, only: writeMatTransposeSize, readMatTransposeSize
      use mod_materials, only: nmaterials
      use mod_nodes, only: nodes
      implicit none
      
      private
      public :: initInteractionData, writeInteractionData,
     &  interactions, readInteractionData, processInteractionData
      
      type interactiondata 
C     read-in
      integer, allocatable :: table(:,:)
C     processed
      integer, allocatable :: mat(:,:)
      end type
      
C     module variables (global)      
      type(interactiondata) :: interactions
      
      contains
************************************************************************     
      subroutine initInteractionData(interactionfile)
      
C     Subroutine: initInteractionData

C     Inputs: interactionfile --- filename where interaction data is stored
C     (should be something like '[filepref]_interactions')

C     Outputs: None

C     Purpose: Read, initialize data in "interactions" structure, which
C     contains information about bonds between atoms of different types
      
      implicit none
      
C     input variables
      character(len=*) :: interactionfile
      
      call readInteractionData(interactionfile)
      call processInteractionData()
      
      end subroutine initInteractionData
************************************************************************     
      subroutine readInteractionData(interactionfile)
      
C     Subroutine: readInteractionData

C     Inputs: interactionfile --- filename where interaction data is stored
C     (should be something like '[filepref]_interactions')

C     Outputs: None

C     Purpose: Read interaction data from file into "interaction" structure
      
      implicit none
      
C     input variables
      character(len=*) :: interactionfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=interactionfile)
      
      call readMatTransposeSize(iunit,interactions%table)
      
      close(iunit)
      
      end subroutine readInteractionData      
************************************************************************
      subroutine processInteractionData()
      
C     Subroutine: processInteractionData

C     Inputs: None

C     Outputs: None

C     Purpose: Create interactions%mat, given interactions%table (see
C     description in "variables/arrays" documentation. Also checks
C     if any interaction is missing
      
C     local variables
      integer :: i, j, n
      integer :: mati, matj, potij
      
      allocate(interactions%mat(nmaterials,nmaterials))

      interactions%mat = -1
      n = size(interactions%table,2)
      do i = 1, n
          mati = interactions%table(1,i)
          matj = interactions%table(2,i)
          potij = interactions%table(3,i)
          interactions%mat(mati,matj) = potij
          interactions%mat(matj,mati) = potij
      end do
      
C     check for uninitialized entries
      do i = 1, nmaterials
          do j = 1, nmaterials
              if (interactions%mat(i,j)==-1) then
                  write(*,*) 'Uninitialized interaction between:'
                  write(*,*) 'Atom ',i,' and atom ',j
                  stop
              end if
          end do
      end do
      
      end subroutine processInteractionData
************************************************************************
      subroutine writeInteractionData(interactionfile)
      
C     Subroutine: writeInteractionData

C     Inputs: interactionfile --- filename where misc data is stored
C     (should be something like '[filepref]_interactions')

C     Outputs: None

C     Purpose: Write interaction data to file (essentially
C     inverse of readInteractionData). Useful in creating "restart" file
C     (i.e. resulting file should be able to be read in by
C     readInteractionData)
      
      implicit none
      
C     input variables
      character(len=*) :: interactionfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=interactionfile)
      
      call writeMatTransposeSize(iunit,interactions%table)
      
      close(iunit)
      
      end subroutine writeInteractionData      
************************************************************************

      end module