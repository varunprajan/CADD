      module mod_disl
      
C     Purpose: Reads/writes/stores information about dislocations,
C     for each fe (i.e. continuum) material. Also contains routines for
C     adding, deleting, and updating dislocations
      
C     Possible extensions: Add more fields to structure, to include
C     other properties that are relevant (e.g. slip planes if obstacles/
C     sources are present). 
      
      use mod_types, only: dp
      implicit none
      
      private
      public :: initDislData, writeDislData, disl, readDislData,
     &          addDislocation, assignDislocation, addImageDislocation,
     &          deleteDislocation, checkTooManyDisl, imagedisl,
     &          readImageDislData, writeImageDislData, initImageDislData
     
      type disldata
C     read-in
      real(dp), allocatable :: posn(:,:)
      integer :: ndisl
      integer, allocatable :: element(:)
      real(dp), allocatable :: localpos(:,:)
      integer, allocatable :: sgn(:)
      end type
      
      type imagedisldata
      integer, allocatable :: mnumfe(:)
      real(dp), allocatable :: posn(:,:)
      integer, allocatable :: sgn(:)
      integer :: nimagedisl
      end type
      
C     module variables
      type(disldata), allocatable :: disl(:)
      type(imagedisldata) :: imagedisl

C     HARD-CODED CONSTANTS
C     max dislocations in one material
      integer, parameter :: nmaxdisl = 1000
C     max image dislocations
      integer, parameter :: nmaximagedisl = 1000
      
      contains
************************************************************************
      subroutine initDislData(dislfile)
      
C     Subroutine: initDislData

C     Inputs: dislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read, initialize data in "disl" structure, which holds
C     information about dislocation positions, orientations, etc. for
C     *each* continuum material

      implicit none
      
C     input variables
      character(len=*) :: dislfile
      
      call readDislData(dislfile)
      
      end subroutine initDislData
************************************************************************
      subroutine initImageDislData(imagedislfile)
      
C     Subroutine: initImageDislData

C     Inputs: imagedislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_imagedisl')

C     Outputs: None

C     Purpose: Read, initialize data in "imagedisl" structure, which holds
C     information about image dislocations (dislocations to account for
C     those that have escaped mesh)

      implicit none
      
C     input variables
      character(len=*) :: imagedislfile
      
      call readImageDislData(imagedislfile)
      
      end subroutine initImageDislData
************************************************************************
      subroutine readDislData(dislfile)
      
C     Subroutine: readDislData

C     Inputs: dislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read dislocation data (positions, etc.) from file,
C     for *each* continuum material, initialize/allocate structures/arrays

C     Notes: Dislocation arrays are allocated with # of rows equal to
C     nmaxdisl, which is a hard-coded constant. If ndisl > nmaxdisl,
C     an error is thrown. An alternative is to reallocate the array
C     every time ndisl > nmaxdisl...

      implicit none
      
C     input variables
      character(len=*) :: dislfile
      
C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: nrow, nfematerials
      
      open(newunit=iunit,file=dislfile)
      
      read(iunit,*) nfematerials
      allocate(disl(nfematerials))
C     this should be identical to nfematerials from mod_fe_elements!
      do i = 1, nfematerials
      
C         element
          read(iunit,*) disl(i)%ndisl
          call checkTooManyDisl(disl(i)%ndisl,nmaxdisl)
          allocate(disl(i)%element(nmaxdisl))
          disl(i)%element = 0
          read(iunit,*) (disl(i)%element(j),j=1,disl(i)%ndisl)

C         local position w/in element          
          read(iunit,*) disl(i)%ndisl, nrow
          allocate(disl(i)%localpos(nrow,nmaxdisl))
          disl(i)%localpos = 0.0_dp
          read(iunit,*) ((disl(i)%localpos(k,j),k=1,nrow),
     &                                          j=1,disl(i)%ndisl)

C         position           
          read(iunit,*) disl(i)%ndisl, nrow
          allocate(disl(i)%posn(nrow,nmaxdisl))
          disl(i)%posn = 0.0_dp
          read(iunit,*) ((disl(i)%posn(k,j),k=1,nrow),j=1,disl(i)%ndisl)

C         dislocation sign
          read(iunit,*) disl(i)%ndisl          
          allocate(disl(i)%sgn(nmaxdisl))
          disl(i)%sgn = 0
          read(iunit,*) (disl(i)%sgn(j),j=1,disl(i)%ndisl)
      end do
      
      close(iunit)
      
      end subroutine readDislData
************************************************************************
      subroutine readImageDislData(imagedislfile)
      
C     Subroutine: readImageDislData

C     Inputs: imagedislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_imagedisl')

C     Outputs: None

C     Purpose: Read image dislocation data (positions, etc.) from file,
C     initialize/allocate structures/arrays

C     Notes: Image dislocation array is allocated with # of rows equal to
C     nmaximagedisl, which is a hard-coded constant. If nimagedisl > nmaximagedisl,
C     an error is thrown. An alternative is to reallocate the array
C     every time nimagedisl > nmaximagedisl...

      implicit none
      
C     input variables
      character(len=*) :: imagedislfile
      
C     local variables
      integer :: iunit
      integer :: j, k
      integer :: nrow
      
      open(newunit=iunit,file=imagedislfile)
      
C     fe material number
      read(iunit,*) imagedisl%nimagedisl
      call checkTooManyDisl(imagedisl%nimagedisl,nmaximagedisl)
      allocate(imagedisl%mnumfe(nmaximagedisl))
      imagedisl%mnumfe = 0
      read(iunit,*) (imagedisl%mnumfe(j),j=1,imagedisl%nimagedisl)

C     position      
      read(iunit,*) imagedisl%nimagedisl, nrow
      allocate(imagedisl%posn(nrow,nmaximagedisl))
      imagedisl%posn = 0.0_dp
      read(iunit,*) ((imagedisl%posn(k,j),k=1,nrow),
     &                                    j=1,imagedisl%nimagedisl)

C     dislocation sign      
      read(iunit,*) imagedisl%nimagedisl
      allocate(imagedisl%sgn(nmaximagedisl))
      imagedisl%sgn = 0
      read(iunit,*) (imagedisl%sgn(j),j=1,imagedisl%nimagedisl)
      
      close(iunit)
      
      end subroutine readImageDislData
************************************************************************
      subroutine addDislocation(mnumfe,element,x,y,theta,bsgn,r,s)

C     Subroutine: addDislocation

C     Inputs: mnumfe --- material number of dislocation
C             element --- element number of dislocation
C             x, y --- global coordinates of dislocation
C             theta --- orientation of dislocation (radians)
C             bsgn --- sign of dislocation (+1 or -1)
C             r, s --- local coordinates of dislocation in element

C     Outputs: None
      
C     Purpose: Adds dislocation to disl structure. First determines
C     empty slot for dislocation in disl structure. Then, assigns
C     dislocation its attributes. Finally, updates ndisl if necessary.
      
      implicit none
      
C     input variables
      integer :: mnumfe, element
      real(dp) :: x, y, theta
      integer :: bsgn
      real(dp) :: r, s
      
C     local variables
      integer :: i
      integer :: dislnum
      
      dislnum = nmaxdisl + 1
      do i = 1, nmaxdisl
          if (disl(mnumfe)%element(i) == 0) then
              dislnum = i
              exit
          end if
      end do
      
      call checkTooManyDisl(dislnum,nmaxdisl)
      
      call assignDislocation(mnumfe,element,x,y,theta,bsgn,r,s,dislnum)
      
      if (dislnum > disl(mnumfe)%ndisl) then
          disl(mnumfe)%ndisl = dislnum
      end if
      
      end subroutine addDislocation
************************************************************************
      subroutine addImageDislocation(x,y,theta,bsgn,mnumfe)

C     Subroutine: addImageDislocation

C     Inputs: x, y --- global coordinates of dislocation
C             theta --- orientation of dislocation (radians)
C             bsgn --- sign of dislocation (+1 or -1)
C             mnumfe --- material number of dislocation

C     Outputs: None
      
C     Purpose: Adds image dislocation to image disl structure. Updates
C     nimagedisl, checks against nmaximagedisl.
      
      implicit none
      
C     input variables
      real(dp) :: x, y
      real(dp) :: theta
      integer :: bsgn
      integer :: mnumfe
      
      imagedisl%nimagedisl = imagedisl%nimagedisl + 1
      call checkTooManyDisl(imagedisl%nimagedisl,nmaximagedisl)
      imagedisl%posn(1,imagedisl%nimagedisl) = x
      imagedisl%posn(2,imagedisl%nimagedisl) = y
      imagedisl%posn(3,imagedisl%nimagedisl) = theta
      imagedisl%sgn(imagedisl%nimagedisl) = bsgn
      imagedisl%mnumfe(imagedisl%nimagedisl) = mnumfe
      
      end subroutine addImageDislocation
************************************************************************
      subroutine assignDislocation(mnumfe,element,x,y,theta,bsgn,
     &                             r,s,dislnum)

C     Subroutine: assignDislocation

C     Inputs: mnumfe --- material number that dislocation belongs to
C             element --- element that dislocation belongs to
C             x, y --- position of dislocation
C             theta --- orientation of dislocation (radians)
C             bsgn --- sign of dislocation (+1 or -1)
C             r, s --- local coordinates of dislocation within element
C             dislnum --- number of dislocation

C     Outputs: None
      
C     Purpose: Assign dislocation with number dislnum to disl structure.
C     Helper routine for addDislocation, updateDislocation
      
      implicit none
      
C     input variables
      integer :: mnumfe, element
      real(dp) :: x, y, theta
      integer :: bsgn
      real(dp) :: r, s
      integer :: dislnum
      
      disl(mnumfe)%posn(1,dislnum) = x
      disl(mnumfe)%posn(2,dislnum) = y
      disl(mnumfe)%posn(3,dislnum) = theta
      disl(mnumfe)%element(dislnum) = element
      disl(mnumfe)%localpos(1,dislnum) = r
      disl(mnumfe)%localpos(2,dislnum) = s
      disl(mnumfe)%sgn(dislnum) = bsgn
      
      end subroutine assignDislocation
************************************************************************
      subroutine deleteDislocation(mnumfe,dislnum)

C     Subroutine: deleteDislocation

C     Inputs: mnumfe --- material number that dislocation belongs to
C             dislnum --- number of dislocation

C     Outputs: None
      
C     Purpose: Remove information about dislocation from disl structure,
C     update ndisl if needed.
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: dislnum
      
      disl(mnumfe)%posn(:,dislnum) = 0.0_dp
      disl(mnumfe)%element(dislnum) = 0
      disl(mnumfe)%localpos(:,dislnum) = 0.0_dp
      disl(mnumfe)%sgn(dislnum) = 0
      
      if (dislnum == disl(mnumfe)%ndisl) then
          disl(mnumfe)%ndisl = dislnum - 1
      end if
      
      end subroutine deleteDislocation
************************************************************************
      subroutine checkTooManyDisl(dislnum,dislmax)

C     Subroutine: checkTooManyDisl

C     Inputs: dislnum --- number of dislocation

C     Outputs: None
      
C     Purpose: Check to see if dislnum exceeds nmaxdisl. If so,
C     throws error and prints (hopefully) helpful message

      implicit none
      
C     input variables
      integer :: dislnum
      integer :: dislmax
      
      if (dislnum > dislmax) then
          write(*,*) 'Number of dislocations is too large'
          write(*,*) 'Increase nmaxdisl or nmaximagedisl'
          stop          
      end if
      
      end subroutine checkTooManyDisl
************************************************************************
      subroutine writeDislData(dislfile)
      
C     Subroutine: writeDislData

C     Inputs: dislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Write dislocation data to file (essentially
C     inverse of readDislData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: dislfile
      
C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: nrow, nfematerials
      
      open(newunit=iunit,file=dislfile)
      
      nfematerials = size(disl)
      write(iunit,*) nfematerials
      do i = 1, nfematerials
          write(iunit,*) disl(i)%ndisl
          do j = 1, disl(i)%ndisl
                write(iunit,*) disl(i)%element(j)
          end do
      
          nrow = size(disl(i)%localpos,1)
          write(iunit,*) disl(i)%ndisl, nrow
          do j = 1, disl(i)%ndisl
                write(iunit,*) (disl(i)%localpos(k,j), k=1,nrow)
          end do
      
          nrow = size(disl(i)%posn,1)
          write(iunit,*) disl(i)%ndisl, nrow
          do j = 1, disl(i)%ndisl
                write(iunit,*) (disl(i)%posn(k,j), k=1,nrow)
          end do
          
          write(iunit,*) disl(i)%ndisl
          do j = 1, disl(i)%ndisl
                write(iunit,*) disl(i)%sgn(j)
          end do
      end do
      
      close(iunit)
      
      end subroutine writeDislData
************************************************************************
      subroutine writeImageDislData(imagedislfile)
      
C     Subroutine: writeImageDislData

C     Inputs: imagedislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_imagedisl')

C     Outputs: None

C     Purpose: Write image dislocation data to file (essentially
C     inverse of readImageDislData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: imagedislfile
      
C     local variables
      integer :: iunit
      integer :: j, k
      integer :: nrow
      
      open(newunit=iunit,file=imagedislfile)
      
      write(iunit,*) imagedisl%nimagedisl
      do j = 1, imagedisl%nimagedisl
            write(iunit,*) imagedisl%mnumfe(j)
      end do

      nrow = size(imagedisl%posn,1)
      write(iunit,*) imagedisl%nimagedisl, nrow
      do j = 1, imagedisl%nimagedisl
            write(iunit,*) (imagedisl%posn(k,j), k=1,nrow)
      end do
      
      write(iunit,*) imagedisl%nimagedisl
      do j = 1, imagedisl%nimagedisl
            write(iunit,*) imagedisl%sgn(j)
      end do
      
      close(iunit)
      
      end subroutine writeImageDislData
************************************************************************   
      
      end module mod_disl