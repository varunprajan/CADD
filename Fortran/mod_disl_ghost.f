      module mod_disl_ghost

C     Purpose: Reads/writes/stores information about ghost dislocations,
C     designed to compensate for dislocations that have passed from/to
C     the atomistic region to/from the continuum region.

C     Long note/TODO: The elasticity problem is difficult to solve exactly in
C     FE if displacement discontinuities are present. The mesh may become
C     highly distorted (in the case of a small element), or the discontinuity
C     may not be resolved (in the case of a large element).

C     So, we instead solve an auxiliary problem, with the displacement discontinuity
C     subtracted out.

C     There are three steps:
C     1) the displacement and the traction fields should be subtracted out;
C     2) the FE problem (which is now "smooth") should be solved;
C     and 3) then the displacement and traction should be added back in.
C     This procedure is permissible if the field satisfies equilibrium
C     and if everything is linear elastic.

C     A discontinuity may arise when one member of a dislocation dipole
C     passes into the atomistic region, leaving the cut of the other
C     uncompensated, because the slip plane intersects the crack plane.
C     Or, it may arise when a dislocation is emitted from the crack tip,
C     since there is no corresponding oppositely-signed dislocation.

C     The tricky part is that a ghost dislocation is only needed when
C     the cut is uncompensated...is this always the case? It seems so,
C     but I need to check more carefully.

C     Also, note that the position of the dislocation within atomistic region 
C     doesn't matter (since the elasticity problem has a unique solution).
C     The only purpose is to cancel/compensate for the singular parts of the solution.

      use mod_types, only: dp
      use mod_disl_try, only: checkTooManyDisl
      use mod_disl_misc, only: dislmisc
      implicit none
      
      private
      public :: ghostdisl, readGhostDislData, writeGhostDislData,
     &  initGhostDislData, addGhostDislocation
          
      type ghostdislt
C     read-in
      integer :: cut
      real(dp) :: posn(2)
      integer :: sgn
      integer :: slipsys
      end type
        
      type ghostdisldata
C     read-in
      type(ghostdislt), allocatable :: list(:)
C     processed
      integer :: nghostdisl
      end type

C     module variables      
      type(ghostdisldata), allocatable :: ghostdisl(:)
      
      contains
************************************************************************
      subroutine initGhostDislData(ghostdislfile)
      
C     Subroutine: initGhostDislData

C     Inputs: ghostdislfile --- filename where ghost dislocation data is stored
C     (should be something like '[filepref]_ghostdisl')

C     Outputs: None

C     Purpose: Read, initialize data in "ghostdisl" structure, which holds
C     information about ghost dislocations (dislocations within atomistic
C     region that account for missing branch cuts of
C     dislocations that have been passed in or out)

      implicit none
      
C     input variables
      character(len=*) :: ghostdislfile
      
      call readGhostDislData(ghostdislfile)
      
      end subroutine initGhostDislData
************************************************************************
      subroutine readGhostDislData(ghostdislfile)
      
C     Subroutine: readGhostDislData

C     Inputs: ghostdislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_ghostdisl')

C     Outputs: None

C     Purpose: Read ghost dislocation data (positions, etc.) from file,
C     initialize/allocate structures/arrays

C     Notes: Dislocation array is allocated with # of rows equal to
C     nmaxghostdisl, which is a hard-coded constant. If nghostdisl > nmaxghostdisl,
C     an error is thrown. An alternative is to reallocate the array
C     every time nghostdisl > nmaxghostdisl...

      implicit none
      
C     input variables
      character(len=*) :: ghostdislfile
      
C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: m, n
      integer :: nfematerials
      
      open(newunit=iunit,file=ghostdislfile)

      read(iunit,*) nfematerials
      allocate(ghostdisl(nfematerials))
      
C     this should be identical to nfematerials from mod_fe_elements!
      do i = 1, nfematerials
          allocate(ghostdisl(i)%list(dislmisc%nmaxghostdisl(i)))

          read(iunit,*) m          
          do j = 1, m
              read(iunit,*) ghostdisl(i)%list(j)%cut
          end do
          
          read(iunit,*) m, n
          do j = 1, m
              read(iunit,*) (ghostdisl(i)%list(j)%posn(k), k = 1,n)
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) ghostdisl(i)%list(j)%sgn
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) ghostdisl(i)%list(j)%slipsys
          end do
          
          ghostdisl(i)%nghostdisl = m
      end do    
      
      close(iunit)
      
      end subroutine readghostDislData
************************************************************************
      subroutine writeGhostDislData(ghostdislfile)
      
C     Subroutine: writeGhostDislData

C     Inputs: ghostdislfile --- filename where data for ghost dislocations is stored
C     (should be something like '[filepref]_ghostdisl')

C     Outputs: None

C     Purpose: Write ghost dislocation data to file (essentially
C     inverse of readGhostDislData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: ghostdislfile
      
C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: m, n
      integer :: nfematerials
      
      open(newunit=iunit,file=ghostdislfile)

      nfematerials = size(ghostdisl)
      write(iunit,*) nfematerials
      do i = 1, nfematerials
          m = ghostdisl(i)%nghostdisl

          write(iunit,*) m
          do j = 1, m
              write(iunit,*) ghostdisl(i)%list(j)%cut
          end do
          write(iunit,*) ''
          
          n = size(ghostdisl(i)%list(1)%posn)
          write(iunit,*) m, n
          do j = 1, m
              write(iunit,*) (ghostdisl(i)%list(j)%posn(k), k = 1,n)
          end do
          write(iunit,*) ''

          write(iunit,*) m          
          do j = 1, m
              write(iunit,*) ghostdisl(i)%list(j)%sgn
          end do
          write(iunit,*) ''

          write(iunit,*) m
          do j = 1, m
              write(iunit,*) ghostdisl(i)%list(j)%slipsys
          end do
          write(iunit,*) ''
          
      end do
      
      close(iunit)
      
      end subroutine writeGhostDislData
************************************************************************
      subroutine addGhostDislocation(x,y,isys,bsgn,bcut,mnumfe)

C     Subroutine: addGhostDislocation

C     Inputs: x, y --- global coordinates of dislocation
C             isys --- slip system for dislocation (in mnumfe)
C             bsgn --- sign of dislocation (+1 or -1)
C             bcut --- branch cut of dislocation (0 if to the left, 1 if to the right)
C             mnumfe --- material number of dislocation

C     Outputs: None
      
C     Purpose: Adds ghost dislocation to ghost disl structure. Updates
C     nghostdisl, checks against nmaxghostdisl.
      
      implicit none
      
C     input variables
      real(dp) :: x, y
      integer :: isys
      integer :: bsgn, bcut
      integer :: mnumfe
      integer :: dislnum
      integer :: nmax
      
      dislnum = ghostdisl(mnumfe)%nghostdisl + 1
      nmax = size(ghostdisl(mnumfe)%list)
      call checkTooManyDisl(dislnum,nmax,'ghostdisl')
      ghostdisl(mnumfe)%list(dislnum)%cut = bcut
      ghostdisl(mnumfe)%list(dislnum)%posn(1) = x
      ghostdisl(mnumfe)%list(dislnum)%posn(2) = y
      ghostdisl(mnumfe)%list(dislnum)%slipsys = isys
      ghostdisl(mnumfe)%list(dislnum)%sgn = bsgn
      ghostdisl(mnumfe)%nghostdisl = dislnum
      
      end subroutine addGhostDislocation
************************************************************************
      end module mod_disl_ghost
