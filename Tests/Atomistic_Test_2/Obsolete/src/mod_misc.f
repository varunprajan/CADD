      module mod_misc

C     Purpose: Reads/writes/stores miscellaneous information about
C     simulation, that doesn't neatly fit in one of the other structures
      
C     Possible extensions: Move some of this stuff to its own structure?
      
      use mod_types, only: dp
      implicit none
      
      private
      public :: misc, writeMiscData, initMiscData, readMiscData
      
      type miscdata
C     read-in
      integer :: dumpincrement
      integer :: incrementcurr
      integer :: increments
      character(len=60) :: potstyle
      integer :: restartincrement
      real(dp) :: timestep
C     processed
      character(len=60) :: simname
      end type

C     module variables      
      type(miscdata) :: misc
      
      contains
************************************************************************
      subroutine initMiscData(miscfile)
      
C     Subroutine: initMiscData

C     Inputs: miscfile --- filename where misc data is stored
C     (should be something like '[filepref]_misc')

C     Outputs: None

C     Purpose: Read, initialize data in "misc" structure, which is a
C     container for information that doesn't fit anywhere else.
      
      implicit none
      
C     input variables
      character(len=*) :: miscfile
      
      call readMiscData(miscfile)
      
      end subroutine initMiscData
************************************************************************
      subroutine readMiscData(miscfile)
      
C     Subroutine: readMiscData

C     Inputs: miscfile --- filename where misc data is stored
C     (should be something like '[filepref]_misc')

C     Outputs: None

C     Purpose: Read misc data from file into "misc" structure
      
      implicit none
      
C     input variables
      character(len=*) :: miscfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=miscfile)
      read(iunit,*) misc%dumpincrement
      read(iunit,*) misc%incrementcurr
      read(iunit,*) misc%increments
      read(iunit,*) misc%potstyle
      read(iunit,*) misc%restartincrement
      read(iunit,*) misc%timestep
      
      close(iunit)
      
      end subroutine readMiscData
************************************************************************
      subroutine writeMiscData(miscfile)
      
C     Subroutine: writeMiscData

C     Inputs: miscfile --- filename where misc data is stored
C     (should be something like '[filepref]_misc')

C     Outputs: None

C     Purpose: Write misc data to file (essentially
C     inverse of readMiscData). Useful in creating "restart" file
C     (i.e. resulting file should be able to be read in by readMiscData)
      
      implicit none
      
C     input variables
      character(len=*) :: miscfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=miscfile)
      write(iunit,*) misc%dumpincrement
      write(iunit,*) misc%incrementcurr
      write(iunit,*) misc%increments
      write(iunit,*) misc%potstyle
      write(iunit,*) misc%restartincrement
      write(iunit,*) misc%timestep

      close(iunit)
      
      end subroutine writeMiscData
************************************************************************
      end module