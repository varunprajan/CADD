      module mod_fe_el_2d
      
C     Purpose: Collection of subroutines related to a single finite element
C     (stiffness matrix, shape functions, Jacobian, gauss points, etc.)
      
C     Possible extensions: Different element types?      
      
      use mod_types, only: dp
      use mod_math, only: invertMat2, getDeterminant2
      implicit none
      
      private
      public :: getK_2d, getB_2d, getBAlt_2d, getElTypeNum, felib,
     &          initFELibrary, getBSub_2d, getJ_2d, findInCPE4,
     &          findInCPE3, checkCPE4, checkCPE3, getN_CPE3_2d,
     &          getN_CPE4_2d, getdN_CPE3_2d, getdN_CPE4_2d,
     &          getGauss_CPE4_2d, getGauss_CPE3_2d
     
      type felibrarydata
C     processed
      character(len=20) :: elname
C     n items
      integer :: nedges
      integer :: nedgenodes
      integer :: nelnodes
      integer :: neldof
      integer :: nelip
      integer :: nedgeip
C     gauss
      real(dp), allocatable :: rgauss(:)
      real(dp), allocatable :: sgauss(:)
      real(dp), allocatable :: wgauss(:)
      real(dp), allocatable :: sedgegauss(:)
      real(dp), allocatable :: wedgegauss(:,:)
C     functions
      procedure(getN_2d), pointer, nopass :: getN_2d_ptr => NULL()
      procedure(getdN_2d), pointer, nopass :: getdN_2d_ptr => NULL()
      procedure(findInElement), pointer, nopass ::
     &                                   findInElement_ptr => NULL()
      procedure(checkElement), pointer, nopass ::
     &                                   checkElement_ptr => NULL()
      end type
      
      integer, parameter :: nfeeltypes = 2
      type(felibrarydata) :: felib(nfeeltypes)
      
      contains
************************************************************************
      subroutine initFELibrary()

C     Subroutine: initFELibrary()

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize data, procedures for each finite element type.
C     Currently, only CPE4 and CPE3 have been implemented.

      implicit none
      
C     local variables
      integer :: eltypenum
      
C     CPE4
      eltypenum = 1
      felib(eltypenum)%elname = 'CPE4'
      felib(eltypenum)%nedges = 4
      felib(eltypenum)%nedgenodes = 2
      felib(eltypenum)%nelnodes = 4
      felib(eltypenum)%neldof = 8
      felib(eltypenum)%nelip = 4
      felib(eltypenum)%nedgeip = 2
      call getGauss_CPE4_2d(felib(eltypenum)%rgauss,
     &                  felib(eltypenum)%sgauss,felib(eltypenum)%wgauss)
      call getGaussEdge_2d(felib(eltypenum)%sedgegauss,
     &                     felib(eltypenum)%wedgegauss)
      felib(eltypenum)%getN_2d_ptr => getN_CPE4_2d
      felib(eltypenum)%getdN_2d_ptr => getdN_CPE4_2d
      felib(eltypenum)%findInElement_ptr => findInCPE4
      felib(eltypenum)%checkElement_ptr => checkCPE4
     
C     CPE3
      eltypenum = 2
      felib(eltypenum)%elname = 'CPE3'
      felib(eltypenum)%nedges = 3
      felib(eltypenum)%nedgenodes = 2
      felib(eltypenum)%nelnodes = 3
      felib(eltypenum)%neldof = 6
      felib(eltypenum)%nelip = 1
      felib(eltypenum)%nedgeip = 2
      call getGauss_CPE3_2d(felib(eltypenum)%rgauss,
     &                  felib(eltypenum)%sgauss,felib(eltypenum)%wgauss)
      call getGaussEdge_2d(felib(eltypenum)%sedgegauss,
     &                     felib(eltypenum)%wedgegauss)
      felib(eltypenum)%getN_2d_ptr => getN_CPE3_2d
      felib(eltypenum)%getdN_2d_ptr => getdN_CPE3_2d
      felib(eltypenum)%findInElement_ptr => findInCPE3
      felib(eltypenum)%checkElement_ptr => checkCPE3
      
      end subroutine initFELibrary
************************************************************************
      function getElTypeNum(elname) result(eltypenum)

C     Function: getElTypeNum

C     Inputs: elname --- name of element type, following ABAQUS convention

C     Outputs: eltypenum --- number of the element type in felib structure
                      
C     Purpose: Convert element type name into element type number; return
C     error if name is not found

C     input variables
      character(len=*) :: elname
      
C     output variables
      integer :: eltypenum
      
C     local variables
      integer :: i

      do i = 1, nfeeltypes
          if (trim(felib(i)%elname) == trim(elname)) then
              eltypenum = i
              return
          end if
      end do
      
C     if we've gotten here, we didn't find the element
      write(*,*) 'Unknown element type'
      write(*,*) 'Currently defined elements:'
      do i = 1, nfeeltypes
          write(*,*) felib(i)%elname
      end do    
      stop
      
      end function getElTypeNum
************************************************************************      
      function getK_2d(posn,dsde,eltypenum) result(K)

C     Subroutine: getK_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             dsde --- stiffness matrix for element, 3 by 3
C             eltypenum --- number of element in element library

C     Outputs: K --- stiffness matrix for element, neldof by neldof

C     Purpose: Evaluate stiffness matrix for element (using Gaussian quadrature)
      
      implicit none
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: dsde(3,3)
      integer :: eltypenum
      
C     output variables
      real(dp) :: K(felib(eltypenum)%neldof,felib(eltypenum)%neldof)
      
C     local variables
      integer :: i
      real(dp) :: B(3,felib(eltypenum)%neldof)
      real(dp) :: J(2,2)
      real(dp) :: detJ
      
      K = 0.0_dp
      do i = 1, felib(eltypenum)%nelip
          call getB_2d(posn,felib(eltypenum)%rgauss(i),
     &                      felib(eltypenum)%sgauss(i),eltypenum,B,J)
          detJ = getDeterminant2(J)
          K = K + felib(eltypenum)%wgauss(i)*
     &            matmul(matmul(transpose(B),dsde),B)*detJ
      end do
      
      end function getK_2d
************************************************************************
      subroutine getB_2d(posn,r,s,eltypenum,B,J)

C     Subroutine: getB_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             r, s --- local (element-level) coordinates
C             eltypenum --- number of element in element library

C     Outputs: B --- matrix involving shape function derivatives, 3 by neldof,
C                    which relates nodal displacements to element strains 
C              J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate strain-displacement matrix for an element at point
C              r, s
      
      implicit none
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: r, s
      integer :: eltypenum

C     output variables
      real(dp) :: B(3,felib(eltypenum)%neldof)
      real(dp) :: J(2,2)

C     local variables
      real(dp) :: dNuglob(2,felib(eltypenum)%neldof)
      real(dp) :: dNvglob(2,felib(eltypenum)%neldof)
      
      call getBSub_2d(posn,r,s,eltypenum,dNuglob,dNvglob,J)
     
      B(1,:) = dNuglob(1,:)
      B(2,:) = dNvglob(2,:)
      B(3,:) = dNuglob(2,:) + dNvglob(1,:)
      
      end subroutine getB_2d
************************************************************************
      subroutine getBAlt_2d(posn,r,s,eltypenum,Balt,J)

C     Subroutine: getBAlt_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             r, s --- local (element-level) coordinates
C             eltypenum --- number of element in element library

C     Outputs: Balt --- matrix involving shape function derivatives, 4 by neldof
C              J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate dU/dX - displacement matrix for an element at point r, s
      
      implicit none
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: r, s
      integer :: eltypenum

C     output variables
      real(dp):: Balt(4,felib(eltypenum)%neldof)
      real(dp) :: J(2,2)

C     local variables
      real(dp) :: dNuglob(2,felib(eltypenum)%neldof)
      real(dp) :: dNvglob(2,felib(eltypenum)%neldof)
      
      call getBSub_2d(posn,r,s,eltypenum,dNuglob,dNvglob,J)

      Balt(1:2,:) = dNuglob
      Balt(3:4,:) = dNvglob
      
      end subroutine getBAlt_2d
************************************************************************
      subroutine getBSub_2d(posn,r,s,eltypenum,dNuglob,dNvglob,J)

C     Subroutine: getBSub_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             r, s --- local (element-level) coordinates
C             eltypenum --- number of element in element library

C     Outputs: dNuglob - [du/dx; du/dy], 2 by neldof
C              dNvglob - [dv/dx, dv/dy], 2 by neldof
C              J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate global derivatives for an element at point r, s
      
      implicit none
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: r, s
      integer :: eltypenum

C     output variables
      real(dp) :: dNuglob(2,felib(eltypenum)%neldof)
      real(dp) :: dNvglob(2,felib(eltypenum)%neldof)
      real(dp) :: J(2,2)

C     local variables
      integer :: i
      real(dp) :: invJ(2,2)
      real(dp), allocatable :: dNdr(:), dNds(:)
      real(dp) :: dNuloc(2,felib(eltypenum)%neldof)
      real(dp) :: dNvloc(2,felib(eltypenum)%neldof)
      
      call felib(eltypenum)%getdN_2d_ptr(r,s,dNdr,dNds)
      J = getJ_2d(posn,dNdr,dNds)
      invJ = invertMat2(J)
      dNuloc = 0.0_dp
      dNvloc = 0.0_dp
      do i = 1, felib(eltypenum)%nelnodes
          dNuloc(1,2*i-1) = dNdr(i)
          dNvloc(1,2*i) = dNdr(i)
          dNuloc(2,2*i-1) = dNds(i)
          dNvloc(2,2*i) = dNds(i)
      end do
      
      dNuglob = matmul(invJ,dNuloc)
      dNvglob = matmul(invJ,dNvloc)
      
      end subroutine getBSub_2d
************************************************************************
      function getJ_2d(posn,dNdr,dNds) result(J)

C     Subroutine: getJ_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             dNdr, dNds --- vector of shape functions derivatives,
C                       of length nelnodes, evaluated at r, s

C     Outputs: J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate Jacobian for an element (tacitly, at a particular
C              point r, s)

      implicit none
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: dNdr(:), dNds(:)
      
C     output variables
      real(dp) :: J(2,2)
      
C     local variables
      integer :: i
      
      do i = 1, 2
          J(1,i) = dot_product(dNdr,posn(i,:))
          J(2,i) = dot_product(dNds,posn(i,:))
      end do

      end function getJ_2d
************************************************************************
      subroutine getdN_2d(r,s,dNdr,dNds)

C     Subroutine: getdN_2d

C     Purpose: Interface for getdN (shape function derivative evaluation)
C              for various element types (see below)

      implicit none
      
C     input variables
      real(dp) :: r, s
      
C     output variables
      real(dp), allocatable :: dNdr(:), dNds(:)
      
      end subroutine getdN_2d
************************************************************************
      subroutine getdN_CPE4_2d(r,s,dNdr,dNds)

C     Subroutine: getdN_CPE4_2d

C     Inputs: r, s --- local (element-level) coordinates

C     Outputs: dNdr, dNds --- vector of shape functions derivatives, evaluated at r, s

C     Purpose: Evaluate all derivatives of all shape functions at r, s for element CPE4

      implicit none

C     input variables
      real(dp) :: r, s
      
C     output variables
      real(dp), allocatable :: dNdr(:), dNds(:)     
     
C     local variables
      real(dp) :: temp1, temp2
      
      temp1 = 0.25_dp*(s - 1.0_dp)
      temp2 = 0.25_dp*(s + 1.0_dp)
      dNdr = [temp1,-temp1,temp2,-temp2]

      temp1 = 0.25_dp*(r - 1.0_dp)
      temp2 = 0.25_dp*(r + 1.0_dp)
      dNds = [temp1,-temp2,temp2,-temp1]
      
      end subroutine getdN_CPE4_2d
************************************************************************
      subroutine getdN_CPE3_2d(r,s,dNdr,dNds)

C     Subroutine: getdN_CPE3_2d

C     Inputs: r, s --- local (element-level) coordinates

C     Outputs: dNdr, dNds --- vector of shape functions derivatives, evaluated at r, s

C     Purpose: Evaluate all derivatives of all shape functions at r, s for element CPE3

      implicit none

C     input variables
      real(dp) :: r, s
      
C     output variables
      real(dp), allocatable :: dNdr(:), dNds(:)

      dNdr = [-1.0_dp,1.0_dp,0.0_dp]
      dNds = [-1.0_dp,0.0_dp,1.0_dp]  
      
      end subroutine getdN_CPE3_2d
************************************************************************
      function getN_2d(r,s) result(N)

C     Function: getN_2d

C     Purpose: Interface for getN (shape function evaluation) for various element types (see below)

      implicit none

C     input variables      
      real(dp) :: r, s
      
C     output variables
      real(dp), allocatable :: N(:)
      
      end function getN_2d
************************************************************************
      function getN_CPE4_2d(r,s) result(N)

C     Function: getN_CPE4_2d

C     Inputs: r, s --- local (element-level) coordinates

C     Outputs: N --- vector of shape functions evaluated at r, s

C     Purpose: Evaluate all shape functions at r, s for element CPE4

C     Notes: Node numbering follows ABAQUS ordering (for simple elements,
C     nodes are ordered counterclockwise)
      
      implicit none

C     input variables      
      real(dp) :: r, s
      
C     output variables
      real(dp), allocatable :: N(:)
      
C     local variables
      real(dp) :: fac1, fac2, fac3, fac4      
      
      fac1 = 1.0_dp - r
      fac2 = 1.0_dp + r
      fac3 = 1.0_dp - s
      fac4 = 1.0_dp + s
      N = 0.25_dp*[fac1*fac3,fac2*fac3,fac2*fac4,fac1*fac4]
      
      end function getN_CPE4_2d
************************************************************************      
      function getN_CPE3_2d(r,s) result(N)

C     Function: getN_CPE3_2d

C     Inputs: r, s --- local (element-level) coordinates

C     Outputs: N --- vector of shape functions evaluated at r, s

C     Purpose: Evaluate all shape functions at r, s for element CPE3

C     Notes: Node numbering follows ABAQUS ordering (for simple elements,
C     nodes are ordered counterclockwise)
      
      implicit none

C     input variables      
      real(dp) :: r, s
      
C     output variables
      real(dp), allocatable :: N(:)
      
      N = [1.0_dp - r - s, r, s]
      
      end function getN_CPE3_2d
************************************************************************      
      subroutine getGauss_CPE4_2d(rgauss,sgauss,wgauss)

C     Subroutine: getGauss_CPE4_2d

C     Inputs: None

C     Outputs: rgauss, sgauss --- vector of r, s-locations for Gaussian quadrature,
C                                 of length nelip
C              wgauss --- vector of weights for Gaussian quadrature, of length nelip

C     Purpose: Return locations, weights for Gaussian quadrature for *area* integral for element CPE4

      implicit none
      
C     output variables
      real(dp), allocatable :: rgauss(:), sgauss(:), wgauss(:)

      rgauss = [0.57735026919_dp,
     &     0.57735026919_dp,-0.57735026919_dp,-0.57735026919_dp]
      sgauss = [0.57735026919_dp,
     &     -0.57735026919_dp,0.57735026919_dp,-0.57735026919_dp]
      wgauss = [1.0_dp,1.0_dp,1.0_dp,1.0_dp]

      end subroutine getGauss_CPE4_2d
************************************************************************      
      subroutine getGauss_CPE3_2d(rgauss,sgauss,wgauss)

C     Subroutine: getGauss_CPE3_2d

C     Inputs: None

C     Outputs: rgauss, sgauss --- vector of r, s-locations for Gaussian quadrature,
C                                 of length nelip
C              wgauss --- vector of weights for Gaussian quadrature, of length nelip

C     Purpose: Return locations, weights for Gaussian quadrature for *area* integral for element CPE3

      implicit none
      
C     output variables
      real(dp), allocatable :: rgauss(:), sgauss(:), wgauss(:)

      rgauss = [0.33333333333_dp]
      sgauss = [0.33333333333_dp]
      wgauss = [0.5_dp]
      
      end subroutine getGauss_CPE3_2d
************************************************************************
      subroutine getGaussEdge_2d(sgauss,wgauss)

C     Subroutine: getGaussEdge_2d

C     Outputs: sgauss --- vector of s-locations (local coordinate along edge)
C                         for Gaussian quadrature
C              wgauss --- array of weights for Gaussian quadrature

C     Purpose: Return locations, weights for Gaussian quadrature for *line* integral
C              (e.g. for calculation of nodal forces from traction)

C     Notes: For CPE4 *or* CPE3

C     Notes: Is accurate for quadratic traction...is this level of accuracy
C     necessary? This was the algorithm used in the original DD code...

C     Algorithm: Uses concept of work-equivalent forces. So, W = 
C     integral(t(s)*u(s),s) = sum(F_i*u_i), where t is the traction on the
C     edge, u is the displacement, and F_i are the nodal forces. For these
C     to be equal for quadratic t(s) and linear u(s) (since the shape function
C     is linear for CPE3 and CPE4), we can evaluate t, u at the Gauss points
C     and use the equation: [F] = L*[wgauss]*[t(sgauss_1) t(sgauss_2)]
      
C     output variables
      real(dp), allocatable :: sgauss(:)  
      real(dp), allocatable :: wgauss(:,:)

      sgauss = [0.21132486540518713_dp, 0.78867513459481287_dp]
      allocate(wgauss(2,2))
      wgauss(1,1) = 0.39433756729740643_dp
      wgauss(1,2) = 0.10566243270259354_dp
      wgauss(2,1) = 0.10566243270259354_dp
      wgauss(2,2) = 0.39433756729740643_dp    
      
      end subroutine getGaussEdge_2d
************************************************************************
      subroutine findInElement(posn,xp,yp,r,s)

C     Subroutine: findInElement

C     Purpose: Interface for findInCPE3, etc. (local coordinates of point within element)
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: xp, yp      

C     in/out variables
      real(dp) :: r, s
      
      end subroutine findInElement
************************************************************************     
      subroutine findInCPE3(posn,xp,yp,r,s)

C     Subroutine: findInCPE3

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             xp, yp --- coordinates of point to search for
C             r, s --- local coordinates of point in element

C     Outputs: r, s --- updated local coordinates of point in element

C     Purpose: Determine local coordinates of point in CPE3 element, using
C              procedure described in Allievi and Bermejo, JCP, 1997
C              (if point does not lie in element, r, s will be outside bounds)
      
      implicit none
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: xp, yp
      
C     in/out variables
      real(dp) :: r, s
      
C     local variables
      real(dp) :: diffx2, diffx3
      real(dp) :: diffy2, diffy3
      real(dp) :: diffxp, diffyp
      real(dp), allocatable :: N(:)
      real(dp) :: xold, yold
      real(dp) :: invdet

      N = getN_CPE3_2d(r,s)
      xold = dot_product(posn(:,1),N)
      yold = dot_product(posn(:,2),N)      
      diffx2 = posn(2,1) - posn(1,1)
      diffx3 = posn(3,1) - posn(1,1)
      diffy2 = posn(2,2) - posn(1,2)
      diffy3 = posn(3,2) - posn(1,2)
      diffxp = xp - xold
      diffyp = yp - yold
      invdet = 1.0_dp/(diffx2*diffy3 - diffx3*diffy2)
      r = r + invdet*(diffy3*diffxp - diffx3*diffyp)
      s = s + invdet*(diffx2*diffyp - diffy2*diffxp)
      
      end subroutine findInCPE3
************************************************************************
      subroutine findInCPE4(posn,xp,yp,r,s)

C     Subroutine: findInCPE4

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             xp, yp --- coordinates of point to search for
C             r, s --- local coordinates of point in element

C     Outputs: r, s --- updated local coordinates of point in element

C     Purpose: Determine local coordinates of point in CPE4 element, using
C              procedure described in Allievi and Bermejo, JCP, 1997
C              (if point does not lie in element, r, s will be outside bounds)

      implicit none
      
C     input variables
      real(dp) :: posn(:,:)
      real(dp) :: xp, yp
      
C     in/out variables
      real(dp) :: r, s
      
C     local variables
      real(dp) :: diffx12, diffx34
      real(dp) :: diffy12, diffy34
      real(dp) :: diffxp, diffyp
      real(dp), allocatable :: N(:)
      real(dp) :: xold, yold
      real(dp) :: invdet
      real(dp) :: a1, a2, a3, b1, b2, b3
      real(dp) :: mat11, mat12, mat21, mat22
      
      N = getN_CPE4_2d(r,s)
      xold = dot_product(posn(:,1),N)
      yold = dot_product(posn(:,2),N)
      diffxp = xp - xold
      diffyp = yp - yold
      diffx12 = posn(1,1) - posn(2,1)
      diffx34 = posn(3,1) - posn(4,1)
      diffy12 = posn(1,2) - posn(2,2)
      diffy34 = posn(3,2) - posn(4,2)
      a1 = 0.25_dp*(diffx34 - diffx12)
      a2 = 0.25_dp*(-posn(1,1) - posn(2,1) + posn(3,1)  + posn(4,1))
      a3 = 0.25_dp*(diffx34 + diffx12)
      b1 = 0.25_dp*(diffy34 - diffy12)
      b2 = 0.25_dp*(-posn(1,2) - posn(2,2) + posn(3,2)  + posn(4,2))
      b3 = 0.25_dp*(diffy34 + diffy12)
      invdet = 1.0_dp/(a1*b2-a2*b1 + (a1*b3-a3*b1)*r + (a3*b2-a2*b3)*s)
      mat11 = b2 + b3*r
      mat12 = a2 + a3*r
      mat21 = b1 + b3*s
      mat22 = a1 + a3*s
      r = r + invdet*(mat11*diffxp - mat12*diffyp)
      s = s + invdet*(mat22*diffyp - mat21*diffxp)
      
      end subroutine findInCPE4
************************************************************************
      function checkElement(edge,r,s) result(check)

C     Function: checkElement

C     Purpose: Interface for checkCPE4, etc. (determine whether a point
C     lies on the correct side of an edge of an element)
      
C     input variables
      integer :: edge
      real(dp) :: r, s
      
C     output variables
      logical :: check
      
      end function checkElement
************************************************************************
      function checkCPE4(edge,r,s) result(check)

C     Function: checkCPE4

C     Inputs: edge --- number of edge of element
C             r, s --- local coordinates of point in element

C     Outputs: check --- logical value, giving whether the local coordinates
C                        lie on the correct side of the edge
                      
C     Purpose: Determine whether local coordinates from findinCPE4
C              lie on the correct side of edge (if checkCPE4(edge,r,s)
C              is true for all edges, then point lies inside element)
      
C     input variables
      integer :: edge
      real(dp) :: r, s
      
C     output variables
      logical :: check
      
      if (edge == 1) then
          check = (s >= -1)
      else if (edge == 2) then
          check = (r <= 1)
      else if (edge == 3) then
          check = (s <= 1)
      else
          check = (r >= -1)
      end if
      
      end function checkCPE4  
************************************************************************
      function checkCPE3(edge,r,s) result(check)

C     Function: checkCPE3

C     Inputs: edge --- number of edge of element
C             r, s --- local coordinates of point in element

C     Outputs: check --- logical value, giving whether the local coordinates
C                        lie on the correct side of the edge
                      
C     Purpose: Determine whether local coordinates from findinCPE3
C              lie on the correct side of edge (if checkCPE3(edge,r,s)
C              is true for all edges, then point lies inside element)
      
C     input variables
      integer :: edge
      real(dp) :: r, s
      
C     output variables
      logical :: check
      
      if (edge == 1) then
          check = (s >= 0)
      else if (edge == 2) then
          check = (r + s <= 1)
      else
          check = (r >= 0)
      end if
      
      end function checkCPE3
************************************************************************
      end module