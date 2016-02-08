      module mod_kdispfield
      
C     Purpose: Has routines for apply K-(displacement) fields
C     to a group of nodes in the body. Currently just isotropic K-field
C     for plane strain has been implemented.

C     Possible extensions: Anisotropic K-fields
      
      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_math, only: piconst
      use mod_groups, only: groups
      implicit none
      
      private
      public :: applyKDispIso
      
      contains

************************************************************************
      subroutine applyKDispIso(KI,KII,mu,nu,xc,yc,gnum)
      
C     Subroutine: applyKDispIso

C     Inputs: KI - mode I stress intensity factor
C             KII - mode II stress intensity factor
C             mu - shear modulus (2d, isotropic)
C             nu - Poisson's ratio (2d, isotropic)
C             xc - x-coordinate of crack
C             yc - y-coordinate of crack
C             gnum - group number
C                    (displacements applied only to nodes in group)

C     Outputs: None
      
C     Purpose: Apply plane strain K field to group of nodes in isotropic body.
C     Update both positions and displacements.
      
C     input variables
      real(dp) :: KI, KII
      real(dp) :: mu, nu
      real(dp) :: xc, yc
      integer :: gnum
      
      integer :: i
      real(dp) :: nodepos(2)
      real(dp) :: cost, sint
      real(dp) :: x, y, r
      real(dp) :: prefac
      real(dp) :: sinthalfsq, sinthalf
      real(dp) :: costhalfsq, costhalf
      real(dp) :: uxInorm, uyInorm
      real(dp) :: uxIInorm, uyIInorm
      real(dp) :: unorm(2), u(2)
      
C     from Bower, Applied Mechanics of Solids, Chapter 9.3
C     valid only for plane strain
      
      prefac = sqrt(1.0_dp/(2.0_dp*piconst))/mu
      do i = 1, nodes%nnodes
C     only do operation for atoms in group
      if (groups(gnum)%maskall(i)) then
          nodepos = nodes%posn(1:2,i)
          x = nodepos(1) - xc
          y = nodepos(2) - yc
          r = sqrt(x*x + y*y)
          cost = x/r
          sint = y/r
C         assume -pi < theta < pi
          sinthalfsq = 0.5_dp*(1.0_dp - cost)
          sinthalf = sqrt(sinthalfsq)
C         since sin is odd
          if (sint < 0.0_dp) then
              sinthalf = -sinthalf
          end if    
          costhalfsq = 0.5_dp*(1.0_dp + cost)
C         since cos is even
          costhalf = sqrt(costhalfsq)
C         disp fields (again, see Bower) without prefactor
          uxInorm = KI*(1.0_dp - 2.0_dp*nu + sinthalfsq)*costhalf
          uyInorm = KI*(2.0_dp - 2.0_dp*nu - costhalfsq)*sinthalf
          uxIInorm = KII*(2.0_dp - 2.0_dp*nu + costhalfsq)*sinthalf
          uyIInorm = KII*(-1.0_dp + 2.0_dp*nu + sinthalfsq)*costhalf
          unorm(1) = uxInorm + uxIInorm
          unorm(2) = uyInorm + uyIInorm
          u = sqrt(r)*prefac*unorm
          nodes%posn(1:2,i) = nodepos + u
          nodes%posn(4:5,i) = nodes%posn(4:5,i) + u
      end if
      end do
      
      end subroutine applyKDispIso
************************************************************************      
      
      end module