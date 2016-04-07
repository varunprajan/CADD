      module mod_dd_integrate

C     Purpose: Integrates dd equation of motion: dx/dt = F_{pk}/b,
C     either using the basic forward Euler method or using the gradient
C     correction method of Chakravarthy and Curtin, 2011, MSMSE. 

C     Notes: See long note for dispFromPKOneMatVCorr describing
C     possible issues with the gradient correction method.
      
      use mod_types, only: dp
      use mod_disl_try, only: disl
      use mod_fe_elements, only: fematerials
      use mod_materials, only: materials
      use mod_disl_fields2, only: getPKTildeStressAll
      use mod_fe_main_2d, only: getFEStressAtPoint
      use mod_slip_sys, only: slipsys, resolveStress
      implicit none
      
      private
      public :: dispFromPK, dispFromPKOneMat, dispFromPKOneMatVCorr,
     &  velFromTau, velFromTauCorr, getResolvedStressOnDisl,
     &  assignDispFromPKOneMat, dispFromVel
     
      procedure(Dummy), pointer :: dispFromPKOneMat_ptr => NULL()
      
C     HARD-CODED CONSTANTS
      real(dp), parameter :: DPOSFAC = 1.0e-6_dp      
      
      contains
************************************************************************
      subroutine assignDispFromPKOneMat(gradientcorr)

C     Inputs: gradientcorr --- true if gradient correction method is to be used

C     Outputs: None

C     Purpose: Assigns pointer for subroutine to compute dislocation
C     displacements based on whether the gradient correction method
C     (Chakravarthy and Curtin, 2011, MSMSE) is to be used
      
      implicit none
      
C     input variables
      logical :: gradientcorr
      
      select case (gradientcorr)
          case (.true.) ! use gradient correction
              dispFromPKOneMat_ptr => dispFromPKOneMatVCorr
          case (.false.) ! don't use gradient correction (i.e., just forward Euler)
              dispFromPKOneMat_ptr => dispFromPKOneMat
      end select
      
      end subroutine assignDispFromPKOneMat
************************************************************************
      subroutine dispFromPK(dt)

C     Inputs: dt --- time increment for DD update

C     Outputs: None

C     Purpose: Computes dislocation displacements for all materials
      
      implicit none
      
C     input variables
      real(dp) :: dt
      
C     local variables
      integer :: i
      
      do i = 1, size(disl)
          call dispFromPKOneMat_ptr(i,dt)
      end do
      
      end subroutine dispFromPK
************************************************************************
      subroutine dispFromPKOneMat(mnumfe,dt)

C     Inputs: mnumfe --- material in which to update dislocation positions
C             dt --- time increment for DD update

C     Outputs: None

C     Purpose: Compute dislocation displacements for one material
C              using Peach-Koehler force and mobility law. Uses hat and
C              tilde stresses to determine Peach-Koehler force.
      
      implicit none
      
C     input variables
      integer :: mnumfe
      real(dp) :: dt
      
C     local variables
      integer :: i
      integer :: mnum
      real(dp) :: tau
      real(dp) :: burgers, disldrag, bfacabs, bfac
      real(dp) :: v, vmax
      integer :: bsgn

      mnum = fematerials%list(mnumfe)
      burgers = materials(mnum)%burgers
      disldrag = materials(mnum)%disldrag
      vmax = materials(mnum)%dislvmax
      bfacabs = burgers/disldrag
      
      do i = 1, disl(mnumfe)%ndisl
      if (disl(mnumfe)%list(i)%active) then ! loop over active dislocations
          tau = getResolvedStressOnDisl(mnumfe,i) 
          bfac = sign(bfacabs,real(bsgn,dp))
          v = velFromTau(tau,bfac)
          disl(mnumfe)%list(i)%disp = dispFromVel(v,dt,vmax)
      end if    
      end do
      
      end subroutine dispfromPKOneMat
************************************************************************
      subroutine dispFromPKOneMatVCorr(mnumfe,dt)

C     Inputs: mnumfe --- material in which to update dislocation positions
C             dt --- time increment for DD update

C     Outputs: None

C     Purpose: Compute dislocation displacements for one material
C              using Peach-Koehler force and mobility law. Uses hat and
C              tilde stresses to determine Peach-Koehler force. Uses
C              velocity correction method described in Chakravarthy and Curtin, MSMSE, 2011

C     Long note: There are several problems.

C     1) It's completely unclear to me whether the algorithm presented in the paper
C     actually fixes the stability problem identify for the forward-Euler method
C     *in general*. (Granted, it works well in the specific example of the dislocation pile-up)
C     The method is not truly an implicit method (like backward-Euler), so I'm
C     not sure why it would be unconditionally stable. Looking at the literature on
C     integration algorithms, it's difficult to even say what type of algorithm
C     this is (semi-implicit, perhaps?).
C     
C     2) Even if the algorithm can be mathematically proven to fix the stability issue,
C     it appears to be written down incorrectly. If the gradients are used
C     to correct the dislocation velocity, the entire *matrix* (Ndisl by Ndisl)
C     should be used: the change in position of the ith disloctaion affects
C     not only the force on it, but also the force on every other dislocation
C     These "cross" terms don't seem to be taken into account. Instead,
C     it is (tacitly) assumed that their affect is negligible, so the 
C     matrix of gradients (the Jacobian) is diagonal. Otherwise,
C     one has to use something like: v(x(t+dt)) = F/B*[I - 1/B*dF_i/dx_j*dt]^{-1}, where I is the identity

C     3) Similarly, the FE fields (hat fields) are also dependent on the dislocation position
C     (for quadrilateral elements). I believe this is ignored in the original code
C     (presumably, because these fields vary slowly), and I have ignored it as well.

C     4) There seems to be a lot of stuff in the code that's not described
C     in the paper. For instance, the algorithm in velFromTauCorr and
C     the underrelaxation (which is a complete mess, and may not even
C     be used any more, etc.)
      
      implicit none
      
C     input variables
      integer :: mnumfe
      real(dp) :: dt
      
C     local variables
      integer :: i
      integer :: isys
      integer :: mnum
      real(dp) :: taubackwards, tauforwards, tau, dtauds
      real(dp) :: burgers, disldrag, bfacabs, bfac
      real(dp) :: v, vmax
      real(dp) :: lengthscale, ds, invds
      real(dp) :: posnsave(2), dpos(2)
      integer :: bsgn

      mnum = fematerials%list(mnumfe)
      burgers = materials(mnum)%burgers
      disldrag = materials(mnum)%disldrag
      vmax = materials(mnum)%dislvmax
      bfacabs = burgers/disldrag
      
      lengthscale = vmax*dt ! this sets length scale
      ds = lengthscale*DPOSFAC ! increment for computing derivative numerically
      invds = 1.0_dp/ds
      
      do i = 1, disl(mnumfe)%ndisl
      if (disl(mnumfe)%list(i)%active) then ! loop over active dislocations
          bsgn = disl(mnumfe)%list(i)%sgn
      
C         get new disl. position (for sake of computing derivative)
          posnsave = disl(mnumfe)%list(i)%posn
          isys = disl(mnumfe)%list(i)%slipsys
          dpos = 0.5_dp*ds*slipsys(mnumfe)%trig(:,isys) ! 1/2*displacement along slip plane
      
C         backwards
          disl(mnumfe)%list(i)%posn = posnsave - dpos
          taubackwards = getResolvedStressOnDisl(mnumfe,i)
          
C         forwards
          disl(mnumfe)%list(i)%posn = posnsave + dpos
          tauforwards = getResolvedStressOnDisl(mnumfe,i)                
          
C         compute derivative numerically
          dtauds = (tauforwards - taubackwards)*invds
          tau = 0.5_dp*(taubackwards + tauforwards)
          
C         flip if bsgn = -1
          bfac = sign(bfacabs,real(bsgn,dp))
          
C         reset position
          disl(mnumfe)%list(i)%posn = posnsave
          
C         compute velocity, displacement
          v = velFromTauCorr(tau,dtauds,bfac,dt)
          disl(mnumfe)%list(i)%disp = dispFromVel(v,dt,vmax)
      end if    
      end do
      
      end subroutine dispfromPKOneMatVCorr
************************************************************************
      subroutine Dummy(mnumfe,dt)
      
C     input variables
      integer :: mnumfe
      real(dp) :: dt
      
      end subroutine
************************************************************************
      function getResolvedStressOnDisl(mnumfe,i) result(tau)

C     Inputs: mnumfe --- fe material number for dislocation
C             i --- index of dislocation in disl(mnumfe)%list

C     Outputs: tau --- shear stress acting on dislocation resolved on slip plane

C     Purpose: Compute resolved shear stress for dislocation i
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: i
      
C     output variables
      real(dp) :: tau
      
C     local variables
      integer :: isys, element
      real(dp) :: r, s
      real(dp) :: stresstilde(3), stresshat(3), stress(3)
      
      isys = disl(mnumfe)%list(i)%slipsys
      element = disl(mnumfe)%list(i)%element
      r = disl(mnumfe)%list(i)%localpos(1)
      s = disl(mnumfe)%list(i)%localpos(2)
      stresstilde = getPKTildeStressAll(i,mnumfe)
      stresshat = getFEStressAtPoint(mnumfe,element,r,s)
      stress = stresshat + stresstilde
      tau = resolveStress(mnumfe,isys,stress)
      
      end function getResolvedStressOnDisl
************************************************************************
      function velFromTau(tau,bfac) result(v)

C     Inputs: tau --- resolved shear stress on dislocation
C             bfac --- = b/B, b = burgers (signed), B = mobility/drag coefficient

C     Outputs: v --- velocity of dislocation along slip plane

C     Purpose: Compute dislocation velocity using resolved shear stress,
C     tau, on dislocation    
      
      implicit none
      
C     input variables
      real(dp) :: tau
      real(dp) :: bfac
      
C     output variables
      real(dp) :: v
      
      v = tau*bfac
      
      end function velfromTau
************************************************************************
      function velFromTauCorr(tau,dtauds,bfac,dt) result(v)

C     Inputs: tau --- resolved shear stress on dislocation
C             dtauds --- derivative of resolved shear stress on dislocation
C                        w.r.t. position along slip plane (s)
C             bfac --- = b/B, b = burgers (signed), B = mobility/drag coefficient
C             bsgn --- sign of dislocation (+1 or -1)
C             dt --- time increment for DD update
C             vmax --- max dislocation velocity

C     Outputs: v --- dislocation velocity along slip plane

C     Purpose: Compute dislocation velocity using resolved shear stress,
C     tau, and derivative of shear stress, dtauds, on dislocation,
C     using velocity correction method (Equation 2, Chakravarthy and Curtin, MSMSE, 2011)

C     Notes: I simply copied the algorithm in the old CADD code (see mod_dd_slip, lines 980-989).
C     I don't understand why this algorithm works well, but it seems to work much better than the 
C     naive approach of not doing anything when v0 and vcorr have opposite signs
      
      implicit none
      
C     input variables
      real(dp) :: tau
      real(dp) :: dtauds
      real(dp) :: bfac
      real(dp) :: dt
      
C     output variables
      real(dp) :: v
      
C     local variables
      real(dp) :: v0, vcorr
      real(dp) :: corrfac
      
C     construct uncorrected, corrected, actual velocity
      v0 = tau*bfac
      corrfac = 1.0_dp - bfac*dt*dtauds
      vcorr = v0/corrfac
      if (corrfac < 0.0_dp) then ! v0 and vcorr have opposite signs
          if (abs(v0) > abs(vcorr)) then
              v = v0 + vcorr
          else
              v = v0
          end if
      else
          v = vcorr
      end if
      
      end function velFromTauCorr
************************************************************************
      function dispFromVel(v,dt,vmax) result(disp)

C     Inputs: v --- velocity of dislocation along slip plane
C             dt --- DD timestep
C             vmax --- max dislocation velocity along slip plane
C                      (either ad hoc underrelaxation
C                       or real, physical velocity limit related to speed of sound)

C     Outputs: disp --- displacement of dislocation along slip plane

C     Purpose: Compute dislocation displacement using dislocation velocity,
C     after enforcing vmax
      
      implicit none
      
C     input variables
      real(dp) :: v
      real(dp) :: dt
      real(dp) :: vmax
      
C     output variables
      real(dp) :: disp
       
      if (abs(v) > vmax) then ! cap velocity
          v = sign(vmax,v)
      end if
      disp = v*dt
      
      end function dispFromVel
************************************************************************
      end module mod_dd_integrate