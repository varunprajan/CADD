      module mod_misc

C     Purpose: Reads/writes/stores miscellaneous information about
C     simulation, that doesn't neatly fit in one of the other structures
      
C     Possible extensions: Move some of this stuff to its own structure?
      
      use mod_types, only: dp
      use mod_math, only: intToLogical, logicalToInt
      implicit none
      
      private
      public :: misc, writeMiscData, initMiscData, readMiscData,
     &          updateMiscIncrementCurr
      
      type miscdata
C     read-in
      integer :: dumpincrement
      integer :: incrementcurr
      integer :: increments
      logical :: iscrackproblem
      character(len=60) :: potstyle
      integer :: restartincrement
      real(dp) :: timestep
C     processed
      character(len=60) :: simname
      character(len=20) :: simtype
      end type

C     module variables      
      type(miscdata) :: misc
      
      contains
************************************************************************
      subroutine initMiscData(miscfile)

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

C     Inputs: miscfile --- filename where misc data is stored
C     (should be something like '[filepref]_misc')

C     Outputs: None

C     Purpose: Read misc data from file into "misc" structure
      
      implicit none
      
C     input variables
      character(len=*) :: miscfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=miscfile)
      read(iunit,*) misc%dumpincrement
      read(iunit,*) misc%incrementcurr
      read(iunit,*) misc%increments
      read(iunit,*) temp
C     use explicit integer -> logical conversion
      misc%iscrackproblem = intToLogical(temp)
      read(iunit,*) misc%potstyle
      read(iunit,*) misc%restartincrement
      read(iunit,*) misc%timestep
      
      close(iunit)
      
      end subroutine readMiscData
************************************************************************
      subroutine writeMiscData(miscfile)

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
      integer :: temp
      
      open(newunit=iunit,file=miscfile)
      write(iunit,*) misc%dumpincrement
      write(iunit,*) misc%incrementcurr
      write(iunit,*) misc%increments
C     use explicit logical -> integer conversion
      temp = logicalToInt(misc%iscrackproblem)
      write(iunit,*) temp
      write(iunit,*) misc%potstyle
      write(iunit,*) misc%restartincrement
      write(iunit,*) misc%timestep

      close(iunit)
      
      end subroutine writeMiscData
************************************************************************
      subroutine updateMiscIncrementCurr(deltaincrement)

C     Inputs: deltaincrement --- change in incrementcurr (incrementcurrnew - incrementcurrold)

C     Outputs: None

C     Purpose: Update current increment...usually, usage is updateIncrementCurr(1), to advance by 1 step
      
      implicit none
      
C     input variables
      integer :: deltaincrement
      
      misc%incrementcurr = misc%incrementcurr + deltaincrement
      
      end subroutine updateMiscIncrementCurr
************************************************************************      
      end module