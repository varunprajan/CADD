      module mod_delaunay

C     Purpose: Computes delaunay triangulation (2D) of a set of points,
C     and includes some auxiliary functions (e.g. circumradius of triangle, etc.)

C     Notes/TODO: Only 2D, currently

      use mod_types, only: dp
      use mod_math, only: thirdconst, getCircumradiusSqforTriangle
      use mod_materials, only: materials
      use mod_neighbors, only: neighbors
      use mod_nodes, only: nodes

      implicit none

      type delaunaydata
      real(dp) :: circumradiussqcutoff
      integer, allocatable :: nodenums(:)
      real(dp), allocatable :: xy(:,:)
      integer :: numtri
      integer, allocatable :: trivert(:,:)
      integer, allocatable :: trineigh(:,:)
      logical, allocatable :: trigood(:)
      end type

      private :: diaedg, i4_modp, i4_sign, perm_check, perm_inverse,
     &           r82vec_permute, r82vec_sort_heap_index_a, swapec,
     &           vbedg, i4_wrap, lrline
      public :: genDelaunay, dtris2, identifyLargeTri, regenDelaunay,
     &          getTriNodes, getTriCenter, CIRCUMSQFACHEX,
     &          setDelaunayPos, genBadTriangles, setDelaunayPosDef,
     &          setDelaunayPosUndef
     
C     HARD-CODED CONSTANTS
      real(dp), parameter :: CIRCUMSQFACHEX = 2.0_dp ! see identifyLargeTri
                                                     ! ratio of circumradius**2 in largest dislocated triangle to circumradius**2 in equilibrium triangle
                                                     ! (must be less than 3, because otherwise edge triangles would be counted as "good")

      contains
************************************************************************
      subroutine regenDelaunay(delaunay,regen)

C     In/out: delaunay --- structure containing information about delaunay triangulation
C             regen --- parameter indicating whether delaunay triangulation needs to be regenerated
C                       (after regeneration, it is set to false)

C     Outputs: None
      
C     Purpose: Update delaunay%xy using new (deformed) node positions. If needed,
C     regenerate triangles/neighbors.

      implicit none
      
C     in/out variables
      type(delaunaydata) :: delaunay
      logical :: regen
      
      call setDelaunayPosDef(delaunay) ! update xy coordinates
      if (regen) then
          call genDelaunay(delaunay) ! regenerate triangulation
          call genBadTriangles(delaunay) ! regenerate bad triangles
          regen = .false. ! reset
      end if 
      
      end subroutine regenDelaunay
************************************************************************
      subroutine setDelaunayPosDef(delaunay)

C     Inputs: None

C     In/out: delaunay --- structure containing information about delaunay triangulation

C     Outputs: None
      
C     Purpose: Set delaunay%xy using new (deformed) positions of nodes listed in delaunay%nodenums.
C     Allocate if necessary.

      implicit none
      
C     in/out variables
      type(delaunaydata) :: delaunay
      
      call setDelaunayPos(delaunay,.false.)
      
      end subroutine setDelaunayPosDef
************************************************************************
      subroutine setDelaunayPosUndef(delaunay)

C     Inputs: None

C     In/out: delaunay --- structure containing information about delaunay triangulation

C     Outputs: None
      
C     Purpose: Set delaunay%xy using undeformed positions of nodes listed in delaunay%nodenums.
C     Allocate if necessary.

      implicit none
      
C     in/out variables
      type(delaunaydata) :: delaunay
      
      call setDelaunayPos(delaunay,.true.)
      
      end subroutine setDelaunayPosUndef
************************************************************************
      subroutine setDelaunayPos(delaunay,undeformed)

C     Inputs: undeformed --- flag indicating whether undeformed positions are to be used

C     In/out: delaunay --- structure containing information about delaunay triangulation

C     Outputs: None
      
C     Purpose: Set delaunay%xy using either undeformed or deformed positions of nodes listed in delaunay%nodenums.
C     Allocate if necessary.
      
      implicit none
      
C     input variables
      logical :: undeformed
      
C     in/out variables
      type(delaunaydata) :: delaunay
      
C     local variables
      integer :: i, npoints
      integer :: node
      real(dp) :: posn(2)
      
C     allocate xy
      if (allocated(delaunay%xy)) then
          deallocate(delaunay%xy)
      end if
      npoints = size(delaunay%nodenums)
      allocate(delaunay%xy(2,npoints))
      
C     update xy
      do i = 1, npoints
          node = delaunay%nodenums(i)
          posn = nodes%posn(1:2,node)
          if (undeformed) then
              posn = posn - nodes%posn(4:5,node)
          end if    
          delaunay%xy(:,i) = posn ! current position
      end do
      
      end subroutine setDelaunayPos
************************************************************************
      subroutine genDelaunay(delaunay)

C     In/out: delaunay --- structure containing information about delaunay triangulation

C     Outputs: None
      
C     Purpose: Generate triangles and neighbors from xy points
      
      implicit none
      
C     in/out variables
      type(delaunaydata) :: delaunay
      
C     local variables
      integer :: numpoints
      
C     deallocate, if necessary
      if (allocated(delaunay%trivert)) then
          deallocate(delaunay%trivert)
      end if
      if (allocated(delaunay%trineigh)) then
          deallocate(delaunay%trineigh)
      end if
      
C     run main routine
      numpoints = size(delaunay%xy,2)
      call dtris2(numpoints,delaunay%xy,
     &            delaunay%trivert,delaunay%trineigh)
      delaunay%numtri = size(delaunay%trivert,2)
      
      end subroutine genDelaunay
************************************************************************      
      subroutine genBadTriangles(delaunay)

C     In/out: delaunay --- structure containing information about delaunay triangulation

C     Outputs: None
      
C     Purpose: Generate list of bad triangles

      implicit none
      
C     in/out variables
      type(delaunaydata) :: delaunay
      
C     deallocate, if necessary
      if (allocated(delaunay%trigood)) then
          deallocate(delaunay%trigood)
      end if
      
C     find bad triangles
      allocate(delaunay%trigood(delaunay%numtri))
      call identifyLargeTri(delaunay) 
      
      end subroutine genBadTriangles      
************************************************************************
      function getTriNodes(delaunay,trinum) result(trinodes)

C     Inputs: delaunay --- structure containing information about delaunay triangulation 
C             trinum --- index of triangle in triangulation

C     Outputs: trinodes --- array, 2 by 3, of coordinates of nodes/vertices of triangle
      
C     Purpose: Get coordinates of nodes/vertices of triangle     
      
      implicit none
      
C     input variables
      type(delaunaydata) :: delaunay
      integer :: trinum
      
C     output variables
      real(dp) :: trinodes(2,3)
      
C     local variables
      integer :: i, v
      
      do i = 1, 3
          v = delaunay%trivert(i,trinum)
          trinodes(:,i) = delaunay%xy(:,v)
      end do
      
      end function getTriNodes
************************************************************************
      function getTriCenter(delaunay,trinum) result(tricenter)

C     Inputs: delaunay --- structure containing information about delaunay triangulation 
C             trinum --- index of triangle in triangulation

C     Outputs: tricenter --- vector, 2 by 1, coordinates of triangle center
      
C     Purpose: Get the center of a triangle (presumed to be location
C     of dislocation, if it is found in that triangle)
      
      implicit none
      
C     input variables
      type(delaunaydata) :: delaunay
      integer :: trinum
    
C     output variables
      real(dp) :: tricenter(2)
      
      tricenter = thirdconst*sum(getTriNodes(delaunay,trinum),2)
      
      end function getTriCenter
************************************************************************
      subroutine identifyLargeTri(delaunay)

C     In/out: delaunay --- structure containing information about delaunay triangulation 

C     Outputs: None
      
C     Purpose: Loop through triangles, labelling each
C     according to whether it is too large or not: i.e. a bad triangle is
C     one whose circumradius exceeds a cutoff
      
      implicit none
      
C     input variables
      type(delaunaydata) :: delaunay
      
C     local variables
      integer :: i
      real(dp) :: p(2,3)
      real(dp) :: circumradiussq
      
      do i = 1, delaunay%numtri
          p = getTriNodes(delaunay,i)
          circumradiussq = getCircumradiusSqForTriangle(p)    
          delaunay%trigood(i) = (circumradiussq <
     &                           delaunay%circumradiussqcutoff)
      end do
      
      end subroutine identifyLargeTri
************************************************************************      
      subroutine dtris2(point_num,point_xy,tri_vert_final,
     &                  tri_nabe_final)
!     
!     ! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!     
!     Discussion:
!     
!     The routine constructs the Delaunay triangulation of a set of 2D vertices
!     using an incremental approach and diagonal edge swaps.  Vertices are
!     first sorted in lexicographically increasing (X,Y) order, and
!     then are inserted one at a time from outside the convex hull.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     25 August 2001
!     
!     Author:
!     
!     Original FORTRAN77 version by Barry Joe.
!     FORTRAN90 version by John Burkardt.
!     
!     Reference:
!     
!     Barry Joe,
!     GEOMPACK - a software package for the generation of meshes
!     using geometric algorithms, 
!     Advances in Engineering Software,
!     Volume 13, pages 325-331, 1991.
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) POINT_NUM, the number of vertices.
!     
!     Input/output, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates 
!     of the vertices.  On output, the vertices have been sorted into 
!     dictionary order.
!     
!     Output, integer ( kind = 4 ) TRI_NUM, the number of triangles in the 
!     triangulation; TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the 
!     number of boundary vertices.
!     
!     Output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make up
!     each triangle.  The elements are indices of POINT_XY.  The vertices of the 
!     triangles are in counter clockwise order.
!     
!     Output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
!     list.  Positive elements are indices of TIL; negative elements are used 
!     for links of a counter clockwise linked list of boundary edges; 
!     LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
!     to the neighbor along edge from vertex J to J+1 (mod 3).
!     
      implicit none

C     input variables
      integer :: point_num
      real(dp) :: point_xy(2,point_num)
      
C     output variables
      integer, allocatable :: tri_vert_final(:,:)
      integer, allocatable :: tri_nabe_final(:,:)
      
C     local variables
      integer :: tri_vert(3,2*point_num)
      integer :: tri_nabe(3,2*point_num)
      integer :: tri_num
      real(dp) :: cmax
      integer :: e, i, ierr, j, k, l, ledg, lr, ltri
      integer :: m, m1, m2, n, redg, rtri, t, top
      integer :: indx(point_num)
      integer :: stack(point_num)
      real(dp) :: tol

      tol = 100.0_dp*epsilon(tol)

      ierr = 0
!     
!     Sort the vertices by increasing (x,y).
      call r82vec_sort_heap_index_a ( point_num, point_xy, indx )
      call r82vec_permute ( point_num, indx, point_xy )
!     
!     Make sure that the data points are "reasonably" distinct.
      m1 = 1
      do i = 2, point_num
          m = m1
          m1 = i
          k = 0
          do j = 1, 2
              cmax = max (abs(point_xy(j,m)),abs(point_xy(j,m1)))
              if (tol*(cmax+1.0_dp) <
     &                 abs(point_xy(j,m)-point_xy(j,m1))) then
                  k = j
                  exit
              end if
          end do
          if ( k == 0 ) then
              write ( *, '(a)' ) 'DTRIS2 - Fatal error 1!'
              stop
          end if
      end do
!     
!     Starting from points M1 and M2, search for a third point M that
!     makes a "healthy" triangle (M1,M2,M)
!     
      m1 = 1
      m2 = 2
      j = 3
      do
          if ( point_num < j ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'DTRIS2 - Fatal error 2!'
              ierr = 225
              return
          end if
          m = j
          lr = lrline (point_xy(1,m),point_xy(2,m),point_xy(1,m1),
     &         point_xy(2,m1), point_xy(1,m2), point_xy(2,m2),0.0_dp)
          if ( lr /= 0 ) then
              exit
          end if
          j = j + 1
      end do
!     
!     Set up the triangle information for (M1,M2,M), and for any other
!     triangles you created because points were collinear with M1, M2.
!     
      tri_num = j - 2
      if ( lr == -1 ) then
          tri_vert(1,1) = m1
          tri_vert(2,1) = m2
          tri_vert(3,1) = m
          tri_nabe(3,1) = -3
          do i = 2, tri_num
              m1 = m2
              m2 = i+1
              tri_vert(1,i) = m1
              tri_vert(2,i) = m2
              tri_vert(3,i) = m
              tri_nabe(1,i-1) = -3 * i
              tri_nabe(2,i-1) = i
              tri_nabe(3,i) = i - 1
          end do
          tri_nabe(1,tri_num) = -3 * tri_num - 1
          tri_nabe(2,tri_num) = -5
          ledg = 2
          ltri = tri_num
      else
          tri_vert(1,1) = m2
          tri_vert(2,1) = m1
          tri_vert(3,1) = m
          tri_nabe(1,1) = -4
          do i = 2, tri_num
              m1 = m2
              m2 = i+1
              tri_vert(1,i) = m2
              tri_vert(2,i) = m1
              tri_vert(3,i) = m
              tri_nabe(3,i-1) = i
              tri_nabe(1,i) = -3 * i - 3
              tri_nabe(2,i) = i - 1
          end do
          tri_nabe(3,tri_num) = -3 * tri_num
          tri_nabe(2,1) = -3 * tri_num - 2
          ledg = 2
          ltri = 1
      end if
!     
!     Insert the vertices one at a time from outside the convex hull,
!     determine visible boundary edges, and apply diagonal edge swaps until
!     Delaunay triangulation of vertices (so far) is obtained.
!     
      top = 0
      do i = j+1, point_num
          m = i
          m1 = tri_vert(ledg,ltri)
          if ( ledg <= 2 ) then
              m2 = tri_vert(ledg+1,ltri)
          else
              m2 = tri_vert(1,ltri)
          end if
          
          lr = lrline(point_xy(1,m), point_xy(2,m), point_xy(1,m1),
     &         point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0_dp)
     
          if ( 0 < lr ) then
              rtri = ltri
              redg = ledg
              ltri = 0
          else
              l = -tri_nabe(ledg,ltri)
              rtri = l / 3
              redg = mod(l,3) + 1
          end if

          call vbedg (point_xy(1,m),point_xy(2,m),point_num,point_xy,
     &                 tri_num,tri_vert,tri_nabe,ltri,ledg,rtri,redg)
          n = tri_num + 1
          l = -tri_nabe(ledg,ltri)
          
          do
              t = l / 3
              e = mod ( l, 3 ) + 1
              l = -tri_nabe(e,t)
              m2 = tri_vert(e,t)
              if ( e <= 2 ) then
                  m1 = tri_vert(e+1,t)
              else
                  m1 = tri_vert(1,t)
              end if
              tri_num = tri_num + 1
              tri_nabe(e,t) = tri_num
              tri_vert(1,tri_num) = m1
              tri_vert(2,tri_num) = m2
              tri_vert(3,tri_num) = m
              tri_nabe(1,tri_num) = t
              tri_nabe(2,tri_num) = tri_num - 1
              tri_nabe(3,tri_num) = tri_num + 1
              top = top + 1
              if (point_num < top) then
                  ierr = 8
                  write ( *, '(a)' ) 'DTRIS2 - Fatal error 3!'
                  stop
              end if
              stack(top) = tri_num
              if ( t == rtri .and. e == redg ) then
                  exit
              end if
          end do
          tri_nabe(ledg,ltri) = -3 * n - 1
          tri_nabe(2,n) = -3 * tri_num - 2
          tri_nabe(3,tri_num) = -l
          ltri = n
          ledg = 2
          call swapec(m,top,ltri,ledg,point_num,point_xy,tri_num,
     &                tri_vert,tri_nabe,stack,ierr)
          if (ierr/=0) then
              write ( *, '(a)' ) 'DTRIS2 - Fatal error 4!'
              stop
          end if
      end do
!     
!     Now account for the sorting that we did.
!     
      do i = 1, 3
          do j = 1, tri_num
              tri_vert(i,j) = indx (tri_vert(i,j))
          end do
      end do

      call perm_inverse(point_num,indx)

      call r82vec_permute(point_num,indx,point_xy)
      
C     Allocate final result
      tri_vert_final = tri_vert(:,1:tri_num)
      tri_nabe_final = tri_nabe(:,1:tri_num)

      return
      end
************************************************************************
      function diaedg (x0,y0,x1,y1,x2,y2,x3,y3) result(res)
!     
!     ! DIAEDG chooses a diagonal edge.
!     
!     Discussion:
!     
!     The routine determines whether 0--2 or 1--3 is the diagonal edge
!     that should be chosen, based on the circumcircle criterion, where
!     (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!     quadrilateral in counterclockwise order.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     19 February 2001
!     
!     Author:
!     
!     Original FORTRAN77 version by Barry Joe.
!     FORTRAN90 version by John Burkardt.
!     
!     Reference:
!     
!     Barry Joe,
!     GEOMPACK - a software package for the generation of meshes
!     using geometric algorithms, 
!     Advances in Engineering Software,
!     Volume 13, pages 325-331, 1991.
!     
!     Parameters:
!     
!     Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!     coordinates of the vertices of a quadrilateral, given in
!     counter clockwise order.
!     
!     Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
!     +1, if diagonal edge 02 is chosen;
!     -1, if diagonal edge 13 is chosen;
!         0, if the four vertices are cocircular.

      implicit none

!     input variables
      real(dp) :: x0, x1, x2, x3
      real(dp) :: y0, y1, y2, y3
      
C     output variables
      integer :: res

!     local variables
      real(dp) :: ca, cb
      real(dp) :: dx10, dx12, dx30, dx32
      real(dp) :: dy10, dy12, dy30, dy32
      real(dp) :: s
      real(dp) :: tol, tola, tolb

      tol = 100.0_dp*epsilon(tol)

      dx10 = x1 - x0
      dy10 = y1 - y0
      dx12 = x1 - x2
      dy12 = y1 - y2
      dx30 = x3 - x0
      dy30 = y3 - y0
      dx32 = x3 - x2
      dy32 = y3 - y2

      tola = tol*max(abs(dx10),abs(dy10),abs(dx30),abs(dy30))
      tolb = tol*max(abs(dx12),abs(dy12),abs(dx32),abs(dy32))

      ca = dx10 * dx30 + dy10 * dy30
      cb = dx12 * dx32 + dy12 * dy32

      if ( tola < ca .and. tolb < cb ) then
          res = -1
      else if ( ca < -tola .and. cb < -tolb ) then
          res = 1
      else
          tola = max(tola,tolb)
          s = ( dx10 * dy30 - dx30 * dy10 ) * cb +
     &        ( dx32 * dy12 - dx12 * dy32 ) * ca
          if ( tola < s ) then
              res = -1
          else if ( s < -tola ) then
              res = 1
          else
              res = 0
          end if
      end if

      end function diaedg
************************************************************************
      function i4_modp(i,j) result(res)
!     
!     ! I4_MODP returns the nonnegative remainder of I4 division.
!     
!     Discussion:
!     
!     If
!          NREM = I4_MODP ( I, J )
!          NMULT = ( I - NREM ) / J
!     then
!          I = J * NMULT + NREM
!     where NREM is always nonnegative.
!     
!     The MOD function computes a result with the same sign as the
!     quantity being divided.  Thus, suppose you had an angle A,
!     and you wanted to ensure that it was between 0 and 360.
!     Then mod(A,360) would do, if A was positive, but if A
!     was negative, your result would be between -360 and 0.
!     
!     On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!     
!     An I4 is an integer ( kind = 4 ) value.
!     
!     Example:
!     
!            I     J     MOD I4_MODP    Factorization
!     
!          107    50       7       7    107 =  2 *  50 + 7
!          107   -50       7       7    107 = -2 * -50 + 7
!         -107    50      -7      43   -107 = -3 *  50 + 43
!         -107   -50      -7      43   -107 =  3 * -50 + 43
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     02 March 1999
!     
!     Author:
!     
!     John Burkardt
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) I, the number to be divided.
!     
!     Input, integer ( kind = 4 ) J, the number that divides I.
!     
!     Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!     divided by J.
!     
      implicit none

C     input variables
      integer :: i, j
      
C     output variables
      integer :: res

      res = mod(i,j)

      if (res < 0) then
          res = res + abs(j)
      end if

      end function i4_modp
************************************************************************
      function i4_sign(x) result(res)
!     
!     ! I4_SIGN evaluates the sign of an I4.
!     
!     Discussion:
!     
!     An I4 is an integer ( kind = 4 ) value.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     27 March 2004
!     
!     Author:
!     
!     John Burkardt
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) X, the number whose sign is desired.
!     
!     Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
!     
      implicit none

C     input variables
      integer :: x
      
C     output variables
      integer :: res

      if ( x < 0 ) then
          res = -1
      else
          res = 1
      end if

      end function i4_sign
************************************************************************
      function i4_wrap (ival,ilo,ihi) result(res)
!     
!     ! I4_WRAP forces an I4 to lie between given limits by wrapping.
!     
!     Discussion:
!     
!     An I4 is an integer ( kind = 4 ) value.
!     
!     Example:
!     
!     ILO = 4, IHI = 8
!     
!     I  Value
!     
!     -2     8
!     -1     4
!         0     5
!         1     6
!         2     7
!         3     8
!         4     4
!         5     5
!         6     6
!         7     7
!         8     8
!         9     4
!     10     5
!     11     6
!     12     7
!     13     8
!     14     4
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     19 August 2003
!     
!     Author:
!     
!     John Burkardt
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) IVAL, a value.
!     
!     Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!     
!     Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!     
      implicit none

C     input variables
      integer :: ival
      integer :: ilo, ihi
      
C     output variables
      integer :: res

C     local variables
      integer :: jhi, jlo, wide

      jlo = min(ilo,ihi)
      jhi = max(ilo,ihi)
      wide = jhi - jlo + 1

      if (wide == 1) then
          res = jlo
      else
          res = jlo + i4_modp(ival-jlo,wide)
      end if

      end function
************************************************************************      
      function lrline (xu,yu,xv1,yv1,xv2,yv2,dv) result(res)
!     
!     ! LRLINE determines if a point is left of, right or, or on a directed line.
!     
!     Discussion:
!     
!     The directed line is parallel to, and at a signed distance DV from
!     a directed base line from (XV1,YV1) to (XV2,YV2).
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     14 July 2001
!     
!     Author:
!     
!     Original FORTRAN77 version by Barry Joe.
!     FORTRAN90 version by John Burkardt.
!     
!     Reference:
!     
!     Barry Joe,
!     GEOMPACK - a software package for the generation of meshes
!     using geometric algorithms,
!     Advances in Engineering Software,
!     Volume 13, pages 325-331, 1991.
!     
!     Parameters:
!     
!     Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
!     position relative to the directed line is to be determined.
!     
!     Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
!     that determine the directed base line.
!     
!     Input, real ( kind = 8 ) DV, the signed distance of the directed line
!     from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!     DV is positive for a line to the left of the base line.
!     
!     Output, integer ( kind = 4 ) LRLINE, the result:
!     +1, the point is to the right of the directed line;
!         0, the point is on the directed line;
!     -1, the point is to the left of the directed line.
!     
      implicit none

C     input variables
      real(dp) :: xu, yu
      real(dp) :: xv1, yv1, xv2, yv2
      real(dp) :: dv
      
C     output variables
      integer :: res
      
C     local variables
      real ( kind = 8 ) dx
      real ( kind = 8 ) dxu
      real ( kind = 8 ) dy
      real ( kind = 8 ) dyu
      real ( kind = 8 ) t
      real ( kind = 8 ) tolabs
      real(dp) :: tol
      
      tol = 100.0_dp*epsilon(tol)

      dx = xv2 - xv1
      dy = yv2 - yv1
      dxu = xu - xv1
      dyu = yu - yv1

      tolabs = tol*max(abs(dx),abs(dy),abs(dxu),abs(dyu),abs(dv))
      t = dy*dxu - dx*dyu + dv*sqrt(dx*dx + dy*dy)

      if ( tolabs < t ) then
          res = 1
      else if ( -tolabs <= t ) then
          res = 0
      else
          res = -1
      end if

      end function lrline
************************************************************************
      subroutine perm_check ( n, p, base, ierror )
!     
!     ! PERM_CHECK checks that a vector represents a permutation.
!     
!     Discussion:
!     
!     The routine verifies that each of the integers from BASE to
!     to BASE+N-1 occurs among the N entries of the permutation.
!     
!     Set the input quantity BASE to 0, if P is a 0-based permutation,
!     or to 1 if P is a 1-based permutation.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     31 October 2008
!     
!     Author:
!     
!     John Burkardt
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) N, the number of entries.
!     
!     Input, integer ( kind = 4 ) P(N), the array to check.
!     
!     Input, integer ( kind = 4 ) BASE, the index base.
!     
!     Output, integer ( kind = 4 ) IERROR, error flag.
!     0, the array represents a permutation.
!     nonzero, the array does not represent a permutation.  The smallest
!     missing value is equal to IERROR.
!     
      implicit none

C     input variables
      integer :: n, base
      integer :: p(n)

C     output variables
      integer :: ierror
      
C     local variables
      integer :: ffind, seek

      ierror = 0

      do seek = base, base + n - 1
          ierror = 1
          do ffind = 1, n
              if ( p(ffind) == seek ) then
                    ierror = 0
                    exit
              end if
          end do
          if ( ierror /= 0 ) then
              write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
              stop
          end if
      end do
      
      end subroutine perm_check
************************************************************************
      subroutine perm_inverse ( n, p )
!     
!     ! PERM_INVERSE inverts a permutation "in place".
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     02 January 2006
!     
!     Author:
!     
!     John Burkardt
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) N, the number of objects being permuted.
!     
!     Input/output, integer ( kind = 4 ) P(N), the permutation, in standard 
!     index form.  On output, P describes the inverse permutation
!     
      implicit none

C     input variables
      integer :: n
      
C     output variables
      integer :: p(n)
      
C     local variables
      integer, parameter :: base = 1
      integer :: i, i0, i1, i2, ierror, is

      call perm_check ( n, p, base, ierror )
      if ( ierror /= 0 ) then
          write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
          stop
      end if
      
      is = 1
      do i = 1, n
          i1 = p(i)
          do while ( i < i1 )
              i2 = p(i1)
              p(i1) = -i2
              i1 = i2
          end do
          is = -i4_sign(p(i))
          p(i) = is*abs(p(i))
      end do

      do i = 1, n
          i1 = - p(i)
          if ( 0 <= i1 ) then
              i0 = i
              do
                  i2 = p(i1)
                  p(i1) = i0
                  if ( i2 < 0 ) then
                      exit
                  end if
                  i0 = i1
                  i1 = i2
              end do
          end if
      end do

      end subroutine perm_inverse
************************************************************************
      subroutine r82vec_permute ( n, p, a )

!     ! R82VEC_PERMUTE permutes an R82VEC in place.
!     
!     Discussion:
!     
!     An R82VEC is an array of pairs of R8 values.
!     
!     The same logic can be used to permute an array of objects of any 
!     arithmetic type, or an array of objects of any complexity.  The only
!     temporary storage required is enough to store a single object.  The number
!     of data movements made is N + the number of cycles of order 2 or more,
!     which is never more than N + N/2.
!     
!     Example:
!     
!     Input:
!     
!          N = 5
!          P = (   2,    4,    5,    1,    3 )
!          A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!              (11.0, 22.0, 33.0, 44.0, 55.0 )
!     
!     Output:
!     
!          A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!                 ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     13 March 2005
!     
!     Author:
!     
!     John Burkardt
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) N, the number of objects.
!     
!     Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!     that the I-th element of the output array should be the J-th
!     element of the input array.  
!     
!     Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!     
      implicit none

      integer, parameter :: dim_num = 2

C     input variables
      integer :: n
      integer :: p(n)
      
C     in/out variables
      real(dp) :: A(dim_num,n)

C     local variables
      real(dp) :: a_temp(dim_num)
      integer, parameter :: base = 1
      integer :: ierror, iget, iput, istart

      call perm_check ( n, p, base, ierror )

      if ( ierror /= 0 ) then
          write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
          stop
      end if
!     
!     Search for the next element of the permutation that has not been used.
!     
      do istart = 1, n

      if ( p(istart) < 0 ) then

          cycle

      else if ( p(istart) == istart ) then

          p(istart) = - p(istart)
          cycle

      else

          a_temp(1:dim_num) = a(1:dim_num,istart)
          iget = istart
!     
!     Copy the new value into the vacated entry.
!     
          do

            iput = iget
            iget = p(iget)

            p(iput) = - p(iput)

            if ( iget < 1 .or. n < iget ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
              write ( *, '(a)' ) 'A permutation index is out of range.'
              write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
              stop
            end if

            if ( iget == istart ) then
              a(1:dim_num,iput) = a_temp(1:dim_num)
              exit
            end if

            a(1:dim_num,iput) = a(1:dim_num,iget)

          end do

      end if

      end do
!     
!     Restore the signs of the entries.
!     
      p(1:n) = - p(1:n)

      end subroutine r82vec_permute
************************************************************************
      subroutine r82vec_sort_heap_index_a ( n, a, indx )
!     
!     ! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
!     
!     Discussion:
!     
!     An R82VEC is an array of R82's.
!     
!     The sorting is not actually carried out.  Rather an index array is
!     created which defines the sorting.  This array may be used to sort
!     or index the array, or to sort or index related arrays keyed on the
!     original array.
!     
!     Once the index array is computed, the sorting can be carried out
!     "implicitly:
!     
!          A(1:2,INDX(1:N)) is sorted,
!     
!     or explicitly, by the call
!     
!          call r82vec_permute ( n, indx, a )
!     
!     after which A(1:2,I), I = 1 to N is sorted.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     08 December 2004
!     
!     Author:
!     
!     John Burkardt
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) N, the number of entries in the array.
!     
!     Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!     
!     Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!     I-th element of the sorted array is A(1:2,INDX(I)).
!     
      implicit none

      integer, parameter :: dim_num = 2

C     input variables
      integer :: n
      real(dp) :: a(dim_num,n)

C     output variables
      integer :: indx(n)
      real(dp) :: aval(dim_num)
      integer :: i, ir, j, l
      integer :: indxt

      if ( n < 1 ) then
          return
      end if

      do i = 1, n
          indx(i) = i
      end do

      if ( n == 1 ) then
          return
      end if

      l = n / 2 + 1
      ir = n
      do
          if ( 1 < l ) then
              l = l - 1
              indxt = indx(l)
              aval(1:dim_num) = a(1:dim_num,indxt)
          else
              indxt = indx(ir)
              aval(1:dim_num) = a(1:dim_num,indxt)
              indx(ir) = indx(1)
              ir = ir - 1
              if ( ir == 1 ) then
                  indx(1) = indxt
                  exit
              end if
          end if

          i = l
          j = l + l
          
          do while ( j <= ir )
              if ( j < ir ) then
                  if (   a(1,indx(j)) <  a(1,indx(j+1)) .or.
     &                ( a(1,indx(j)) == a(1,indx(j+1)) .and.
     &                  a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
                      j = j + 1
                  end if
              end if

              if (   aval(1) <  a(1,indx(j)) .or.
     &              ( aval(1) == a(1,indx(j)) .and.
     &                aval(2) <  a(2,indx(j)) ) ) then
                  indx(i) = indx(j)
                  i = j
                  j = j + j
              else
                  j = ir + 1
              end if

          end do

          indx(i) = indxt

      end do

      end subroutine r82vec_sort_heap_index_a
************************************************************************      
      subroutine swapec ( i, top, btri, bedg, point_num, point_xy,
     &         tri_num, tri_vert, tri_nabe, stack, ierr )
!     
!     ! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!     
!     Discussion:
!     
!     The routine swaps diagonal edges in a 2D triangulation, based on
!     the empty circumcircle criterion, until all triangles are Delaunay,
!     given that I is the index of the new vertex added to the triangulation.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     14 July 2001
!     
!     Author:
!     
!     Original FORTRAN77 version by Barry Joe.
!     FORTRAN90 version by John Burkardt.
!     
!     Reference:
!     
!     Barry Joe,
!     GEOMPACK - a software package for the generation of meshes
!     using geometric algorithms, 
!     Advances in Engineering Software,
!     Volume 13, pages 325-331, 1991.
!     
!     Parameters:
!     
!     Input, integer ( kind = 4 ) I, the index of the new vertex.
!     
!     Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
!     On output, TOP is zero.
!     
!     Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are 
!     the triangle and edge indices of a boundary edge whose updated indices
!     must be recorded.  On output, these may be updated because of swaps.
!     
!     Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!     
!     Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
!     of the points.
!     
!     Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!     
!     Input/output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle 
!     incidence list.  May be updated on output because of swaps.
!     
!     Input/output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle 
!     neighbor list; negative values are used for links of the counter-clockwise 
!     linked list of boundary edges;  May be updated on output because of swaps.
!          LINK = -(3*I + J-1) where I, J = triangle, edge index.
!     
!     Workspace, integer ( kind = 4 ) STACK(MAXST); on input, entries 1 through
!     TOP contain the indices of initial triangles (involving vertex I)
!     put in stack; the edges opposite I should be in interior;  entries
!     TOP+1 through MAXST are used as a stack.
!     
!     Output, integer ( kind = 4 ) IERR is set to 8 for abnormal return.
!     
      implicit none

C     input variables
      integer :: i
      integer :: point_num
      integer :: tri_num
      real(dp) :: point_xy(2,point_num)
      
C     in/out variables
      integer :: top
      integer :: btri, bedg
      integer :: tri_vert(3,tri_num)
      integer :: tri_nabe(3,tri_num)
      
C     output variables
      integer :: ierr
      
C     local variables
      integer :: a, b, c, e, ee, em1, ep1
      integer :: f, fm1, fp1, l, r, s, swap, t, tt, u
      real(dp) :: x, y
      integer :: stack(point_num)
!     
!     Determine whether triangles in stack are Delaunay, and swap
!     diagonal edge of convex quadrilateral if not.
!     
      x = point_xy(1,i)
      y = point_xy(2,i)

      do
      if ( top <= 0 ) then
          exit
      end if

      t = stack(top)
      top = top - 1

      if ( tri_vert(1,t) == i ) then
          e = 2
          b = tri_vert(3,t)
      else if ( tri_vert(2,t) == i ) then
          e = 3
          b = tri_vert(1,t)
      else
          e = 1
          b = tri_vert(2,t)
      end if

      a = tri_vert(e,t)
      u = tri_nabe(e,t)

      if ( tri_nabe(1,u) == t ) then
          f = 1
          c = tri_vert(3,u)
      else if ( tri_nabe(2,u) == t ) then
          f = 2
          c = tri_vert(1,u)
      else
          f = 3
          c = tri_vert(2,u)
      end if

      swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a),point_xy(1,c),
     &               point_xy(2,c), point_xy(1,b), point_xy(2,b) )

      if ( swap == 1 ) then

          em1 = i4_wrap ( e - 1, 1, 3 )
          ep1 = i4_wrap ( e + 1, 1, 3 )
          fm1 = i4_wrap ( f - 1, 1, 3 )
          fp1 = i4_wrap ( f + 1, 1, 3 )

          tri_vert(ep1,t) = c
          tri_vert(fp1,u) = i
          r = tri_nabe(ep1,t)
          s = tri_nabe(fp1,u)
          tri_nabe(ep1,t) = u
          tri_nabe(fp1,u) = t
          tri_nabe(e,t) = s
          tri_nabe(f,u) = r

          if ( 0 < tri_nabe(fm1,u) ) then
            top = top + 1
            stack(top) = u
          end if

          if ( 0 < s ) then

            if ( tri_nabe(1,s) == u ) then
              tri_nabe(1,s) = t
            else if ( tri_nabe(2,s) == u ) then
              tri_nabe(2,s) = t
            else
              tri_nabe(3,s) = t
            end if

            top = top + 1

            if ( point_num < top ) then
              ierr = 8
              return
            end if

            stack(top) = t

          else

            if ( u == btri .and. fp1 == bedg ) then
              btri = t
              bedg = e
            end if

            l = - ( 3 * t + e - 1 )
            tt = t
            ee = em1

            do while ( 0 < tri_nabe(ee,tt) )

              tt = tri_nabe(ee,tt)

              if ( tri_vert(1,tt) == a ) then
                ee = 3
              else if ( tri_vert(2,tt) == a ) then
                ee = 1
              else
                ee = 2
              end if

            end do

            tri_nabe(ee,tt) = l

          end if

          if ( 0 < r ) then

            if ( tri_nabe(1,r) == t ) then
              tri_nabe(1,r) = u
            else if ( tri_nabe(2,r) == t ) then
              tri_nabe(2,r) = u
            else
              tri_nabe(3,r) = u
            end if

          else

            if ( t == btri .and. ep1 == bedg ) then
              btri = u
              bedg = f
            end if

            l = - ( 3 * u + f - 1 )
            tt = u
            ee = fm1

            do while ( 0 < tri_nabe(ee,tt) )

              tt = tri_nabe(ee,tt)

              if ( tri_vert(1,tt) == b ) then
                ee = 3
              else if ( tri_vert(2,tt) == b ) then
                ee = 1
              else
                ee = 2
              end if

            end do

            tri_nabe(ee,tt) = l

          end if

      end if
      end do

      end subroutine swapec
************************************************************************
      subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert,
     &  tri_nabe, ltri, ledg, rtri, redg )
!     
!     ! VBEDG determines which boundary edges are visible to a point.
!     
!     Discussion:
!     
!     The point (X,Y) is assumed to be outside the convex hull of the
!     region covered by the 2D triangulation.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     25 August 2001
!     
!     Author:
!     
!     Original FORTRAN77 version by Barry Joe.
!     FORTRAN90 version by John Burkardt.
!     
!     Reference:
!     
!     Barry Joe,
!     GEOMPACK - a software package for the generation of meshes
!     using geometric algorithms, 
!     Advances in Engineering Software,
!     Volume 13, pages 325-331, 1991.
!     
!     Parameters:
!     
!     Input, real ( kind = 8 ) X, Y, the coordinates of a point outside
!     the convex hull of the current triangulation.
!     
!     Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!     
!     Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates 
!     of the vertices.
!     
!     Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!     
!     Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle incidence 
!     list.
!     
!     Input, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
!     list; negative values are used for links of a counter clockwise linked 
!     list of boundary edges;
!          LINK = -(3*I + J-1) where I, J = triangle, edge index.
!     
!     Input/output, integer ( kind = 4 ) LTRI, LEDG.  If LTRI /= 0 then these 
!     values are assumed to be already computed and are not changed, else they 
!     are updated.  On output, LTRI is the index of boundary triangle to the 
!     left of the leftmost boundary triangle visible from (X,Y), and LEDG is 
!     the boundary edge of triangle LTRI to the left of the leftmost boundary
!     edge visible from (X,Y).  1 <= LEDG <= 3.
!     
!     Input/output, integer ( kind = 4 ) RTRI.  On input, the index of the 
!     boundary triangle to begin the search at.  On output, the index of the 
!     rightmost boundary triangle visible from (X,Y).
!     
!     Input/output, integer ( kind = 4 ) REDG, the edge of triangle RTRI that 
!     is visible from (X,Y).  1 <= REDG <= 3.
!     
      implicit none

C     input variables
      integer :: point_num
      real(dp) :: x, y
      real(dp) :: point_xy(2,point_num)
      integer :: tri_num
      integer :: tri_vert(3,tri_num)
      integer :: tri_nabe(3,tri_num)

C     in/out variables
      integer :: ltri, ledg, rtri, redg
      
C     local variables
      integer :: a, b, e, l, lr, t
      logical :: ldone
!     
!     Find the rightmost visible boundary edge using links, then possibly
!     leftmost visible boundary edge using triangle neighbor information.
!     
      if ( ltri == 0 ) then
          ldone = .false.
          ltri = rtri
          ledg = redg
      else
          ldone = .true.
      end if

      do

      l = -tri_nabe(redg,rtri)
      t = l / 3
      e = mod ( l, 3 ) + 1
      a = tri_vert(e,t)

      if ( e <= 2 ) then
          b = tri_vert(e+1,t)
      else
          b = tri_vert(1,t)
      end if

      lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b),
     &     point_xy(2,b), 0.0_dp )

      if ( lr <= 0 ) then
          exit
      end if

      rtri = t
      redg = e

      end do

      if ( ldone ) then
          return
      end if

      t = ltri
      e = ledg

      do

      b = tri_vert(e,t)
      e = i4_wrap ( e-1, 1, 3 )

      do while ( 0 < tri_nabe(e,t) )

          t = tri_nabe(e,t)

          if ( tri_vert(1,t) == b ) then
            e = 3
          else if ( tri_vert(2,t) == b ) then
            e = 1
          else
            e = 2
          end if

      end do

      a = tri_vert(e,t)

      lr = lrline(x,y,point_xy(1,a), point_xy(2,a), point_xy(1,b),
     &      point_xy(2,b), 0.0_dp)

      if ( lr <= 0 ) then
          exit
      end if

      end do

      ltri = t
      ledg = e

      end subroutine vbedg
************************************************************************      
      end module
