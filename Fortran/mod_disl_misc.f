      module mod_disl_misc
      
      use mod_types, only: dp
      use mod_math, only: intToLogical, logicalToInt
      use mod_utils, only: writeVecSize, readVecSize
      implicit none

      private
      public :: dislmisc, writeDislMiscData, initDislMiscData,
     & readDislMiscData
      
      type dislmiscdata
C     read-in
      logical :: gradientcorrection
      integer, allocatable :: nmaxdisl(:)
      integer, allocatable :: nmaxdislslip(:)
      integer, allocatable :: nmaxescapeddisl(:)
      integer, allocatable :: nmaxghostdisl(:)
      integer, allocatable :: nmaxobsslip(:)
      integer, allocatable :: nmaxsrcslip(:)
      end type

C     module variables      
      type(dislmiscdata) :: dislmisc
      
      contains
************************************************************************
      subroutine initDislMiscData(dislmiscfile)
      
C     Subroutine: initDislMiscData

C     Inputs: dislmiscfile --- filename where misc DD data is stored
C     (should be something like '[filepref]_dislmisc')

C     Outputs: None

C     Purpose: Read, initialize data in "dislmisc" structure, which is a
C     container for constants used in DD (max number of dislocations,
C     max velocity, etc.)
      
      implicit none
      
C     input variables
      character(len=*) :: dislmiscfile
      
      call readDislMiscData(dislmiscfile)
      
      end subroutine initDislMiscData
************************************************************************
      subroutine readDislMiscData(dislmiscfile)
      
C     Subroutine: readDislMiscData

C     Inputs: dislmiscfile --- filename where misc DD data is stored
C     (should be something like '[filepref]_dislmisc')

C     Outputs: None

C     Purpose: Read misc DD data from file into "dislmisc" structure
      
      implicit none
      
C     input variables
      character(len=*) :: dislmiscfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=dislmiscfile)

      read(iunit,*) temp
C     use explicit integer -> logical conversion
      dislmisc%gradientcorrection = intToLogical(temp)
      call readVecSize(iunit,dislmisc%nmaxdisl)
      call readVecSize(iunit,dislmisc%nmaxdislslip)
      call readVecSize(iunit,dislmisc%nmaxescapeddisl)
      call readVecSize(iunit,dislmisc%nmaxghostdisl)
      call readVecSize(iunit,dislmisc%nmaxobsslip)
      call readVecSize(iunit,dislmisc%nmaxsrcslip)
      
      close(iunit)
      
      end subroutine readDislMiscData
************************************************************************
      subroutine writeDislMiscData(dislmiscfile)
      
C     Subroutine: writeDislMiscData

C     Inputs: dislmiscfile --- filename where misc DD data is stored
C     (should be something like '[filepref]_dislmisc')

C     Outputs: None

C     Purpose: Write misc DD data to file (essentially
C     inverse of readDislMiscData). Useful in creating "restart" file
C     (i.e. resulting file should be able to be read in by readDislMiscData)
      
      implicit none
      
C     input variables
      character(len=*) :: dislmiscfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=dislmiscfile)

C     use explicit logical -> integer conversion
      temp = logicalToInt(dislmisc%gradientcorrection)
      write(iunit,*) temp
      call writeVecSize(iunit,dislmisc%nmaxdisl)
      call writeVecSize(iunit,dislmisc%nmaxdislslip)
      call writeVecSize(iunit,dislmisc%nmaxescapeddisl)
      call writeVecSize(iunit,dislmisc%nmaxghostdisl)
      call writeVecSize(iunit,dislmisc%nmaxobsslip)
      call writeVecSize(iunit,dislmisc%nmaxsrcslip)
      
      close(iunit)
      
      end subroutine writeDislMiscData
************************************************************************
      end module mod_disl_misc