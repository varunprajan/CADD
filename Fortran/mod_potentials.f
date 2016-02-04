      module mod_potentials
      
C     Purpose: Stores information about potentials, in an array of
C     structures (each structure contains information about a different
C     potential).
      
C     Possible extensions: Add more fields to structure, so that other
C     types of potentials can be used (i.e. other than tabular pair
C     potentials).

C     Notes: pottable is stored in the conventional (n by 4) format, as
C     opposed to the transposed (4 by n) format, because columns of
C     pottable are accessed by forceLookup
      
      use mod_types, only: dp
      use mod_utils, only: writeMatSize, readMatSize
      implicit none
      
      private
      public :: potentials, npotentials,
     &          initPotentialData, writePotentialData, readPotentialData
      
      type potentialdata
C     read-in
      real(dp) :: forcecutoff
      character(len=20) :: pname
      real(dp), allocatable :: pottable(:,:)
C     processed
      real(dp) :: forcecutoffsq
      end type
      
C     module variables (global)
      type(potentialdata), allocatable :: potentials(:)
      integer :: npotentials
      
      contains
************************************************************************
      subroutine initPotentialData(potentialfile)

C     Subroutine: initPotentialData

C     Inputs: potentialfile - filename where potential data is stored
C     (should be something like '[filepref]_potentials')

C     Outputs: None

C     Purpose: Read, initialize data in "potentials" structure 

C     Notes: Currently, "processPotentialData" is not defined, but this
C     might be useful later
      
      implicit none

C     input variables
      character(len=*) :: potentialfile
      
      call readPotentialData(potentialfile)
      call processPotentialData()
      
      end subroutine initPotentialData      
************************************************************************      
      subroutine readPotentialData(potentialfile)
      
C     Subroutine: readPotentialData

C     Inputs: potentialfile - filename where potential data is stored
C     (should be something like '[filepref]_potentials')

C     Outputs: None

C     Purpose: Read potential data from file
      
      implicit none

C     input variables
      character(len=*) :: potentialfile
      
C     local variables
      integer :: iunit, i
      
      open(newunit=iunit,file=potentialfile)
      
      read(iunit,*) npotentials
      allocate(potentials(npotentials))
      do i = 1, npotentials
          read(iunit,*) potentials(i)%forcecutoff
          read(iunit,*) potentials(i)%pname
          call readMatSize(iunit,potentials(i)%pottable)
      end do
      
      close(iunit)
      
      end subroutine readPotentialData
************************************************************************
      subroutine processPotentialData()
      
C     Subroutine: processPotentialData

C     Inputs: None

C     Outputs: None

C     Purpose: Create auxiliary quantities (forcecutoff squared) 
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, npotentials
          potentials(i)%forcecutoffsq = (potentials(i)%forcecutoff)**2
      end do    
      
      end subroutine processPotentialData
************************************************************************      
      subroutine writePotentialData(potentialfile)
      
C     Subroutine: writePotentialData

C     Inputs: potentialfile --- filename where potential data is stored
C     (should be something like '[filepref]_potentials')

C     Outputs: None

C     Purpose: Write potential data to file (essentially
C     inverse of readPotentialData). Useful in creating "restart" file
C     (i.e. resulting file should be able to be read in by readPotentialData)
      
      implicit none

C     input variables
      character(len=*) :: potentialfile
      
C     local variables
      integer :: iunit, i
      
      open(newunit=iunit,file=potentialfile)
      
      write(iunit,*) npotentials
      do i = 1, npotentials
          write(iunit,*) potentials(i)%forcecutoff
          write(iunit,*) potentials(i)%pname
          call writeMatSize(iunit,potentials(i)%pottable)
          write(iunit,*) ''
      end do
      
      close(iunit)
      
      end subroutine writePotentialData
************************************************************************
      
      end module