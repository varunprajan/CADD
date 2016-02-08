      module mod_fe_el_2d
      
C     Purpose: Collection of subroutines related to a single finite element
C     (stiffness matrix, shape functions, Jacobian, gauss points, etc.)
      
C     Possible extensions: Different element types?      
      
      use mod_types, only: dp
      use mod_math, only: invertmat2
      use mod_utils, only: prettyPrintRealMat
      implicit none
      
      private
      public :: getK_2d, getGaussEdge_2d, getB_2d, getN_2d,
     &          getdN_2d, getGauss_2d, getJ_2d, getBAlt_2d, getBSub_2d
      
      contains
************************************************************************      
      function getK_2d(posn,dsde,neldof,nelnodes,nelip,elname,
     &                 rgauss,sgauss,wgauss) result(K)

C     Subroutine: getK_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             dsde --- stiffness matrix for element, 3 by 3
C             neldof --- number of degrees of freedom per element (should be 2*nelnodes)
C             nelnodes --- number of nodes per element
C             nelip --- number of integration points per element
C             elname --- name of element, character string, follows ABAQUS convention

C     Outputs: K --- stiffness matrix for element, neldof by neldof

C     Purpose: Evaluate stiffness matrix for element (using Gaussian quadrature)
      
      implicit none
      
C     input variables
      real(dp) :: posn(2,nelnodes)
      real(dp) :: dsde(3,3)
      integer :: neldof, nelnodes, nelip
      character(len=*) :: elname
      real(dp) :: rgauss(nelip), sgauss(nelip), wgauss(nelip)
      
C     output variables
      real(dp) :: K(neldof,neldof)
      
C     local variables
      integer :: i
      real(dp) :: B(3,neldof)
      real(dp) :: J(2,2)
      real(dp) :: detJ
      
      K = 0
      do i = 1, nelip
          call getB_2d(posn,rgauss(i),sgauss(i),
     &                 neldof,nelnodes,elname,B,J)
          detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
          K = K + wgauss(i)*matmul(matmul(transpose(B),dsde),B)*detJ
      end do
      
      end function getK_2d
************************************************************************
      subroutine getB_2d(posn,r,s,neldof,nelnodes,elname,B,J)

C     Subroutine: getB_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             r, s --- local (element-level) coordinates
C             neldof --- number of degrees of freedom per element (should be 2*nelnodes)
C             nelnodes --- number of nodes per element
C             elname --- name of element, character string, follows ABAQUS convention

C     Outputs: B --- matrix involving shape function derivatives, 3 by neldof,
C                    which relates nodal displacements to element strains 
C              J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate strain-displacement matrix for an element at point
C              r, s
      
      implicit none
      
C     input variables
      real(dp) :: posn(2,nelnodes)
      real(dp) :: r, s
      integer :: neldof, nelnodes
      character(len=*) :: elname

C     output variables
      real(dp) :: B(3,neldof)
      real(dp) :: J(2,2)

C     local variables
      real(dp) :: dNuglob(2,neldof), dNvglob(2,neldof)
      
      call getBSub_2d(posn,r,s,neldof,nelnodes,elname,dNuglob,dNvglob,J)
      
      B(1,:) = dNuglob(1,:)
      B(2,:) = dNvglob(2,:)
      B(3,:) = dNuglob(2,:) + dNvglob(1,:)
      
      end subroutine getB_2d
************************************************************************
      subroutine getBAlt_2d(posn,r,s,neldof,nelnodes,elname,Balt,J)

C     Subroutine: getBAlt_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             r, s --- local (element-level) coordinates
C             neldof --- number of degrees of freedom per element (should be 2*nelnodes)
C             nelnodes --- number of nodes per element
C             elname --- name of element, character string, follows ABAQUS convention

C     Outputs: Balt --- matrix involving shape function derivatives, 4 by neldof
C              J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate dU/dX - displacement matrix for an element at point r, s
      
      implicit none
      
C     input variables
      real(dp) :: posn(2,nelnodes)
      real(dp) :: r, s
      integer :: neldof, nelnodes
      character(len=*) :: elname

C     output variables
      real(dp) :: Balt(4,neldof)
      real(dp) :: J(2,2)

C     local variables
      real(dp) :: dNuglob(2,neldof), dNvglob(2,neldof)
      
      call getBSub_2d(posn,r,s,neldof,nelnodes,elname,dNuglob,dNvglob,J)
      
      Balt(1:2,:) = dNuglob
      Balt(3:4,:) = dNvglob
      
      end subroutine getBAlt_2d
************************************************************************
      subroutine getBSub_2d(posn,r,s,neldof,nelnodes,elname,
     &                      dNuglob,dNvglob,J)

C     Subroutine: getBSub_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             r, s --- local (element-level) coordinates
C             neldof --- number of degrees of freedom per element (should be 2*nelnodes)
C             nelnodes --- number of nodes per element
C             elname --- name of element, character string, follows ABAQUS convention

C     Outputs: dNuglob - [du/dx; du/dy], 2 by neldof
C              dNvglob - [dv/dx, dv/dy], 2 by neldof
C              J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate global derivatives for an element at point r, s
      
      implicit none
      
C     input variables
      real(dp) :: posn(2,nelnodes)
      real(dp) :: r, s
      integer :: neldof, nelnodes
      character(len=*) :: elname

C     output variables
      real(dp) :: dNuglob(2,neldof), dNvglob(2,neldof)
      real(dp) :: J(2,2)

C     local variables
      integer :: i
      real(dp) :: invJ(2,2)
      real(dp) :: dNdr(nelnodes), dNds(nelnodes)
      real(dp) :: dNuloc(2,neldof), dNvloc(2,neldof)
      
      call getdN_2d(r,s,elname,nelnodes,dNdr,dNds)
      J = getJ_2d(posn,dNdr,dNds,nelnodes)
      invJ = invertmat2(J)
      dNuloc = 0
      dNvloc = 0
      do i = 1, nelnodes
          dNuloc(1,2*i-1) = dNdr(i)
          dNuloc(2,2*i-1) = dNds(i)
          dNvloc(1,2*i) = dNdr(i)
          dNvloc(2,2*i) = dNds(i)
      end do
      
      dNuglob = matmul(invJ,dNuloc)
      dNvglob = matmul(invJ,dNvloc)
      
      end subroutine getBSub_2d
************************************************************************
      function getJ_2d(posn,dNdr,dNds,nelnodes) result(J)

C     Subroutine: getJ_2d

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             dNdr, dNds --- vector of shape functions derivatives,
C                       of length nelnodes, evaluated at r, s
C             nelnodes --- number of nodes per element

C     Outputs: J --- jacobian matrix, 2 by 2, for the element
C                    (tacitly evaluated at r, s through dNdr, dNds)

C     Purpose: Evaluate Jacobian for an element (tacitly, at a particular
C              point r, s)

      implicit none
      
C     input variables
      integer :: nelnodes
      real(dp) :: posn(2,nelnodes)
      real(dp) :: dNdr(nelnodes), dNds(nelnodes)
      
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
      subroutine getdN_2d(r,s,elname,nelnodes,dNdr,dNds)

C     Subroutine: getdN_2d

C     Inputs: r, s --- local (element-level) coordinates
C             elname --- name of element, character string, follows ABAQUS convention
C             nelnodes --- number of nodes per element

C     Outputs: dNdr, dNds --- vector of shape functions derivatives,
C                       of length nelnodes, evaluated at r, s

C     Purpose: Evaluate shape function derivatives at r, s

C     Notes: Node numbering follows ABAQUS ordering (for simple elements,
C     nodes are ordered counterclockwise)

      implicit none
      
C     input variables
      real(dp) :: r, s
      character(len=*) :: elname
      integer :: nelnodes
      
C     output variables
      real(dp) :: dNdr(nelnodes), dNds(nelnodes)
      
      if (elname == 'CPE4') then
          dNdr(1) = 0.25_dp*(s - 1.0_dp)
          dNdr(2) = -dNdr(1)
          dNdr(3) = 0.25_dp*(s + 1.0_dp)
          dNdr(4) = -dNdr(3)
      
          dNds(1) = 0.25_dp*(r - 1.0_dp)
          dNds(3) = 0.25_dp*(r + 1.0_dp)
          dNds(2) = -dNds(3)
          dNds(4) = -dNds(1)
      else if (elname == 'CPE3') then
          dNdr(1) = -1.0_dp
          dNdr(2) = 1.0_dp
          dNdr(3) = 0.0_dp
      
          dNds(1) = -1.0_dp
          dNds(2) = 0.0_dp
          dNds(3) = 1.0_dp
      end if    
      
      end subroutine getdN_2d
************************************************************************
      function getN_2d(r,s,elname,nelnodes) result(N)

C     Function: getN_2d

C     Inputs: r, s --- local (element-level) coordinates
C             elname --- name of element, character string, follows ABAQUS convention
C             nelnodes --- number of nodes per element

C     Outputs: N --- vector of shape functions of length nelnodes evaluated
C                    at r, s

C     Purpose: Evaluate all shape functions at r, s

C     Notes: Node numbering follows ABAQUS ordering (for simple elements,
C     nodes are ordered counterclockwise)

      implicit none

C     input variables      
      real(dp) :: r, s
      character(len=*) :: elname
      integer :: nelnodes
      
C     output variables
      real(dp) :: N(nelnodes)
      
C     local variables
      real(dp) :: fac1, fac2, fac3, fac4
      
      if (elname == 'CPE4') then
          fac1 = 1.0_dp - r
          fac2 = 1.0_dp + r
          fac3 = 1.0_dp - s
          fac4 = 1.0_dp + s
          N(1) = 0.25_dp*fac1*fac3
          N(2) = 0.25_dp*fac2*fac3
          N(3) = 0.25_dp*fac2*fac4
          N(4) = 0.25_dp*fac1*fac4
      else if (elname == 'CPE3') then
          N(1) = 1.0_dp - r - s
          N(2) = r
          N(3) = s
      end if
      
      end function getN_2d
************************************************************************
      subroutine getGauss_2d(elname,nelip,rgauss,sgauss,wgauss)

C     Subroutine: getGauss_2d

C     Inputs: elname --- name of element, character string, follows ABAQUS convention
C             nelip --- number of integration points per element

C     Outputs: rgauss, sgauss --- vector of r, s-locations for Gaussian quadrature,
C                                 of length nelip
C              wgauss --- vector of weights for Gaussian quadrature, of length nelip

C     Purpose: Return locations, weights for Gaussian quadrature for *area* integral

      implicit none

C     input variables
      character(len=*) :: elname
      integer :: nelip
      
C     output variables
      real(dp) :: rgauss(nelip), sgauss(nelip), wgauss(nelip)

      if (elname == 'CPE4') then
          rgauss = [0.57735026919_dp,
     &     0.57735026919_dp,-0.57735026919_dp,-0.57735026919_dp]
          sgauss = [0.57735026919_dp,
     &     -0.57735026919_dp,0.57735026919_dp,-0.57735026919_dp]
          wgauss = [1.0_dp,1.0_dp,1.0_dp,1.0_dp]
      else if (elname == 'CPE3') then
          rgauss = [0.33333333333_dp]
          sgauss = [0.33333333333_dp]
          wgauss = [0.5_dp]
      end if
      
      end subroutine getGauss_2d
************************************************************************
      subroutine getGaussEdge_2d(elname,nedgenodes,nedgeip,
     &                           sgauss,wgauss)

C     Subroutine: getGaussEdge_2d

C     Inputs: elname --- name of element, character string, follows ABAQUS convention
C             nedgenodes --- number of nodes on edge of element
C             nedgeip --- number of integration points on edge of element

C     Outputs: sgauss --- vector of s-locations (local coordinate along edge)
C                         for Gaussian quadrature, of length nedgeip
C              wgauss --- array of weights for Gaussian quadrature,
C                         nedgenodes by nedgeip

C     Purpose: Return locations, weights for Gaussian quadrature for *line* integral
C              (e.g. for calculation of nodal forces from traction)

C     Notes: Is accurate for quadratic traction...is this level of accuracy
C     necessary? This was the algorithm used in the original DD code...

C     Algorithm: Uses concept of work-equivalent forces. So, W = 
C     integral(t(s)*u(s),s) = sum(F_i*u_i), where t is the traction on the
C     edge, u is the displacement, and F_i are the nodal forces. For these
C     to be equal for quadratic t(s) and linear u(s) (since the shape function
C     is linear for CPE3 and CPE4), we can evaluate t, u at the Gauss points
C     and use the equation: [F] = L*[wgauss]*[t(sgauss_1) t(sgauss_2)]
      
C     input variables
      character(len=*) :: elname
      integer :: nedgenodes
      integer :: nedgeip
      
C     output variables
      real(dp) :: sgauss(nedgeip)  
      real(dp) :: wgauss(nedgenodes,nedgeip)

      sgauss = [0.21132486540518713_dp, 0.78867513459481287_dp]
      if ((elname == 'CPE3').or.
     &    (elname == 'CPE4')) then
          wgauss(1,1) = 0.39433756729740643_dp
          wgauss(1,2) = 0.10566243270259354_dp
          wgauss(2,1) = 0.10566243270259354_dp
          wgauss(2,2) = 0.39433756729740643_dp
      end if    
      
      end subroutine getGaussEdge_2d
************************************************************************
      end module