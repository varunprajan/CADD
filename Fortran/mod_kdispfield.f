      module mod_kdispfield
      
C     Purpose: Has routines for apply K-(displacement) fields
C     to a group of nodes in the body. Currently just isotropic K-field
C     for plane strain has been implemented.

C     Possible extensions: Anisotropic K-fields
      
      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_math, only: piconst
      use mod_groups, only: groups, getGroupNum
      use mod_utils, only: prettyPrintVec
      implicit none
      
      private
      public :: applyKDispIsoSub, applyKDispIsoSet, applyKDispIsoIncr
      
      contains

************************************************************************
      subroutine applyKDispIsoSet(KI,KII,mu,nu,xc,yc,gname)

C     Inputs: KI - mode I stress intensity factor
C             KII - mode II stress intensity factor
C             mu - shear modulus (2d, isotropic)
C             nu - Poisson's ratio (2d, isotropic)
C             xc - x-coordinate of crack
C             yc - y-coordinate of crack
C             gname - group name
C                    (displacements applied only to nodes in group)

C     Outputs: None
      
C     Purpose: Set group of nodes to have plane strain K field in isotropic body.
C     I.e. x_{curr,new} = x_{ref} + u

      implicit none
      
C     input variables
      real(dp) :: KI, KII
      real(dp) :: mu, nu
      real(dp) :: xc, yc
      character(len=*) :: gname      
      
      call applyKDispIsoSub(KI,KII,mu,nu,xc,yc,gname,.false.)
      
      end subroutine
************************************************************************
      subroutine applyKDispIsoIncr(KIincr,KIIincr,mu,nu,xc,yc,gname)

C     Inputs: KIincr - mode I stress intensity factor increment
C             KIIincr - mode II stress intensity factor increment
C             mu - shear modulus (2d, isotropic)
C             nu - Poisson's ratio (2d, isotropic)
C             xc - x-coordinate of crack
C             yc - y-coordinate of crack
C             gname - group name
C                    (displacements applied only to nodes in group)

C     Outputs: None
      
C     Purpose: Incrementally apply plane strain K field to group of nodes in isotropic body.
C     I.e. x_{curr,new} = x_{curr} + u

      implicit none
      
C     input variables
      real(dp) :: KIincr, KIIincr
      real(dp) :: mu, nu
      real(dp) :: xc, yc
      character(len=*) :: gname      
      
      call applyKDispIsoSub(KIincr,KIIincr,mu,nu,xc,yc,gname,.true.)
      
      end subroutine
************************************************************************
      subroutine applyKDispIsoSub(KI,KII,mu,nu,xc,yc,gname,incr)

C     Inputs: KI - mode I stress intensity factor
C             KII - mode II stress intensity factor
C             mu - shear modulus (2d, isotropic)
C             nu - Poisson's ratio (2d, isotropic)
C             xc - x-coordinate of crack
C             yc - y-coordinate of crack
C             gname - group name
C                    (displacements applied only to nodes in group)
C             incr - flag indicating whether displacements should be incrementally applied or not (see below)

C     Outputs: None
      
C     Purpose: Apply plane strain K field to group of nodes in isotropic body.
C     Applies field based on undeformed positions; x_def = x_undef + K*sqrt(2*pi*r)*f(theta); y_def = y_undef + K*sqrt(2*pi*r)*g(theta)
      
      implicit none
      
C     input variables
      real(dp) :: KI, KII
      real(dp) :: mu, nu
      real(dp) :: xc, yc
      character(len=*) :: gname
      logical :: incr
      
C     local variables
      integer :: gnum
      integer :: i
      real(dp) :: posnundef(2)
      real(dp) :: cost, sint
      real(dp) :: x, y, r, invr
      real(dp) :: prefac
      real(dp) :: sinthalfsq, sinthalf
      real(dp) :: costhalfsq, costhalf
      real(dp) :: uxInorm, uyInorm
      real(dp) :: uxIInorm, uyIInorm
      real(dp) :: unorm(2), u(2)
      
C     from Bower, Applied Mechanics of Solids, Chapter 9.3
C     valid only for plane strain
      
      gnum = getGroupNum(gname)
      prefac = sqrt(1.0_dp/(2.0_dp*piconst))/mu
      do i = 1, nodes%nnodes
C     only do operation for atoms in group
      if (groups(gnum)%maskall(i)) then
          posnundef = nodes%posn(1:2,i) - nodes%posn(4:5,i)
          x = posnundef(1) - xc
          y = posnundef(2) - yc
          r = sqrt(x*x + y*y)
          invr = 1.0_dp/r
          cost = x*invr
          sint = y*invr
C         assume -pi < theta < pi
          sinthalfsq = 0.5_dp*(1.0_dp - cost)
          sinthalf = sign(sqrt(sinthalfsq),sint) ! since sin is odd
          costhalfsq = 0.5_dp*(1.0_dp + cost) 
          costhalf = sqrt(costhalfsq) ! since cos is even
C         disp fields (again, see Bower) without prefactor
          uxInorm = KI*(1.0_dp - 2.0_dp*nu + sinthalfsq)*costhalf
          uyInorm = KI*(2.0_dp - 2.0_dp*nu - costhalfsq)*sinthalf
          uxIInorm = KII*(2.0_dp - 2.0_dp*nu + costhalfsq)*sinthalf
          uyIInorm = KII*(-1.0_dp + 2.0_dp*nu + sinthalfsq)*costhalf
          unorm(1) = uxInorm + uxIInorm
          unorm(2) = uyInorm + uyIInorm
          u = (sqrt(r)*prefac)*unorm
          if (incr) then
              nodes%posn(1:2,i) = nodes%posn(1:2,i) + u
              nodes%posn(4:5,i) = nodes%posn(4:5,i) + u
          else
              nodes%posn(1:2,i) = posnundef + u
              nodes%posn(4:5,i) = u
          end if
      end if
      end do
      
      end subroutine applyKDispIsoSub
************************************************************************
      end module