      module mod_disl_fields2
      
C     Purpose: Computes fields (tilde fields, in van der Giessen and
C     Needleman's terminology) created by all dislocations acting in an
C     *infinite* medium. Also computes fields from "image" dislocations,
C     which have left the mesh, and contribute only plastic displacement
C     (slip) and no stress.
      
C     Possible extensions:
C     1) Currently, expressions for fields assume
C     isotropy. Although our toy potentials are isotropic, this is not
C     the case for many materials.

C     2) Calculates fields directly, which has a total cost of O(N^3) (?)
C     The multipole method (see Chakravarthy and Curtin, MSMSE, 2011)
C     can reduce this to O(N), but it is more complicated to implement
C     and not necessarily faster unless ndisl is large
C     (few hundred or more?). This is also important if the velocity
C     correction in Chakravarthy, 2011 is employed (since derivatives
C     of fields are required.)

C     3) Assumes all dislocations are in same material w/ same elastic
C     constants. This is fine, but extension to multiple grains is
C     tricky and requires extra bookkeeping. Basically, the DD
C     problem in *each* grain is treated independently (see O'Day and
C     Curtin, JAM, 2004). So, the only relevant dislocations are
C     those in material i for a point (x,y) in material i.
      
      use mod_types, only: dp
      use mod_math, only: tolconst, piconst, rotateVec2d, rotateStress2d
      use mod_disl, only: disl, imagedisl
      use mod_materials, only: materials
      use mod_fe_elements, only: fematerials
      implicit none
      
      private
      public :: getTildeDispAtPointAll, getTildeStressatPointAll,
     &            getPKTildeStressAll, getDispAtPointSub,
     &            getStressAtPointSub, adjustDxnDyn,
     &            getDispAtPoint, getStressAtPoint, getImageDispAtPoint,
     &            getImageDispAtPointSub, getDispAtPointHelper
      
      contains
************************************************************************
      function getTildeDispAtPointAll(posn,mnumfe) result(disp)
      
C     Subroutine: getTildeDispAtPointAll

C     Inputs: posn --- position at which tilde displacement field is sought
C                      (vector, length 2)
C             mnumfe --- fe material in which point lies

C     Outputs: disp --- tilde displacement at point (vector, length 2)

C     Purpose: Get displacement field from all dislocations *in material*
C              at point of interest (tilde field --- does not account for finite boundaries)
      
      implicit none
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      
C     output variables
      real(dp) :: disp(2)
      
C     local variables
      integer :: i
      integer :: mnum
      real(dp) :: dislpos(2)
      real(dp) :: theta
      integer :: bsgn
      real(dp) :: dispnew(2)
      
      disp = 0.0_dp
      
C     loop over dislocations in material
      mnum = fematerials%list(mnumfe)
      do i = 1, disl(mnumfe)%ndisl
C     an actual dislocation
      if (disl(mnumfe)%element(i) /= 0) then
          dislpos = disl(mnumfe)%posn(1:2,i)
          theta = disl(mnumfe)%posn(3,i)
          bsgn = disl(mnumfe)%sgn(i)
          dispnew = getDispAtPoint(posn,dislpos,theta,bsgn,mnum)
          disp = disp + dispnew
      end if    
      end do
      
C     TODO: check if it is correct to apply these fields to all materials!
C     loop over image dislocations
      do i = 1, imagedisl%nimagedisl
          dislpos = imagedisl%posn(1:2,i)
          theta = imagedisl%posn(3,i)
          bsgn = imagedisl%sgn(i)
          mnumfe = imagedisl%mnumfe(i)
          mnum = fematerials%list(mnumfe)
          dispnew = getImageDispAtPoint(posn,dislpos,theta,bsgn,mnum)
          disp = disp + dispnew
      end do
      
      end function getTildeDispAtPointAll
************************************************************************
      function getTildeStressAtPointAll(posn,mnumfe) result(stress)
      
C     Subroutine: getTildeStressAtPointAll

C     Inputs: posn --- position at which tilde displacement field is sought
C                      (vector, length 2)
C             mnumfe --- fe material in which point lies

C     Outputs: stress --- tilde stress at point (vector, length 3 --- sxx, syy, sxy)

C     Purpose: Get stress field from all dislocations *in material*
C              at point of interest (tilde field --- does not account for finite boundaries)
      
      implicit none
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      integer :: i
      integer :: mnum
      real(dp) :: dislpos(2)
      real(dp) :: theta
      integer :: bsgn
      real(dp) :: stressnew(3)
      
C     loop over dislocations in material      
      stress = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, disl(mnumfe)%ndisl
C     an actual dislocation
      if (disl(mnumfe)%element(i) /= 0) then
          dislpos = disl(mnumfe)%posn(1:2,i)
          theta = disl(mnumfe)%posn(3,i)
          bsgn = disl(mnumfe)%sgn(i)
          stressnew = getStressAtPoint(posn,dislpos,theta,bsgn,mnum)
          stress = stress + stressnew
      end if    
      end do
      
      end function getTildeStressAtPointAll
************************************************************************
      function getPKTildeStressAll(dislnum,mnumfe) result(stress)
      
C     Subroutine: getPKTildeStressAll

C     Inputs: dislnum --- number of dislocation for which Peach-Koehler stress is sought
C             mnumfe --- fe material in which point lies

C     Outputs: stress --- Peach-Koehler stress at point (vector, length 3 --- sxx, syy, sxy)

C     Purpose: Get contribution *of dislocations in material*
C              to Peach-Koehler stress (total PK stress also includes FE contribution)

C     Notes: See expression in parentheses in Equation 12 in vdG and Needleman, MSMSE, 1995.
      
      implicit none
      
C     input variables
      integer :: dislnum
      integer :: mnumfe
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      integer :: i
      integer :: mnum
      real(dp) :: posn(2), dislpos(2)
      real(dp) :: theta
      integer :: bsgn
      real(dp) :: stressnew(3)
      
C     loop over dislocations in material      
      stress = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, disl(mnumfe)%ndisl
C     an actual dislocation
      if (disl(mnumfe)%element(i) /= 0) then
C         no self force
          if (i /= dislnum) then
              dislpos = disl(mnumfe)%posn(1:2,i)
              theta = disl(mnumfe)%posn(3,i)
              bsgn = disl(mnumfe)%sgn(i)
              stressnew = getStressAtPoint(posn,dislpos,theta,bsgn,mnum)
              stress = stress + stressnew
          end if
      end if
      end do
      
      end function getPKTildeStressAll
************************************************************************
      function getDispAtPoint(posn,dislpos,theta,bsgn,mnum) result(disp)
      
C     Subroutine: getDispAtPoint

C     Inputs: posn - 2 by 1 vector of position of point at which displacement is sought
C             dislpos - 2 by 1 vector of position of dislocation
C             theta - angle of dislocation
C             bsgn - sign of dislocation (+1 or -1)
C             mnum - number of material (*not* fe material) that point and (edge) dislocation
C                    reside within

C     Outputs: disp --- displacement at point (vector, 2 by 1)

C     Purpose: Get tilde displacement from a *single* dislocation

      implicit none

C     input variables
      real(dp) :: posn(2), dislpos(2)
      real(dp) :: theta
      integer :: mnum
      integer :: bsgn

C     output variables
      real(dp) :: disp(2)
     
      disp = getDispAtPointHelper(posn,dislpos,theta,bsgn,mnum,.false.)
      
      end function getDispAtPoint
************************************************************************
      function getImageDispAtPoint(posn,dislpos,theta,bsgn,
     &                              mnum) result(disp)
      
C     Subroutine: getImageDispAtPoint

C     Inputs: posn - 2 by 1 vector of position of point at which displacement is sought
C             dislpos - 2 by 1 vector of position of dislocation
C             theta - angle of dislocation
C             bsgn - sign of dislocation (+1 or -1)
C             mnum - number of material (*not* fe material) that point and (edge) dislocation
C                    reside within

C     Outputs: disp --- displacement at point (vector, 2 by 1)

C     Purpose: Get (plastic) displacement from a *single* *image* dislocation

      implicit none

C     input variables
      real(dp) :: posn(2), dislpos(2)
      real(dp) :: theta
      integer :: mnum
      integer :: bsgn

C     output variables
      real(dp) :: disp(2)
     
      disp = getDispAtPointHelper(posn,dislpos,theta,bsgn,mnum,.true.)
     
      end function getImageDispAtPoint
************************************************************************
      function getDispAtPointHelper(posn,dislpos,theta,bsgn,
     &                              mnum,imageoption) result(disp)

C     Function: getDispAtPointHelper

C     Inputs: posn - 2 by 1 vector of position of point at which displacement is sought
C             dislpos - 2 by 1 vector of position of dislocation
C             theta - angle of dislocation
C             bsgn - sign of dislocation (+1 or -1)
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             imageoption --- flag indicating whether dislocation is an image
                        
C     Outputs: disp - 2 by 1 vector of displacements of point (in original
C                     coordinate system, *not* dislocation coordinate system)
     
C     Purpose: Compute displacements at point due to a single edge
C              dislocation. Dislocation can be either a real dislocation
C              (in which case, elasticity expressions apply) or image dislocation
C              (which contributes only plastic field).
C              Performs appropriate rotations so that
C              elasticity expressions can be used. Also adjusts for points
C              that are too close to the dislocation core.

C     Notes: 'n' in dxn, etc. refers to the 'new' (dislocation)
C            coordinate system; 'o' refers to the old (unrotated) system.
      
      implicit none
      
C     input variables
      real(dp) :: posn(2), dislpos(2)
      real(dp) :: theta
      integer :: mnum
      integer :: bsgn
      logical :: imageoption
      
C     output variables
      real(dp) :: disp(2)
      
C     local variables
      real(dp) :: dxo, dyo, dxn, dyn
      real(dp) :: uxo, uyo, uxn, uyn
      real(dp) :: cost, sint
      
      cost = cos(theta)
      sint = sin(theta)
      dxo = posn(1) - dislpos(1)
      dyo = posn(2) - dislpos(2)
      
C     rotate coordinates into new coordinate system
      call rotateVec2d(cost,sint,dxo,dyo,dxn,dyn)
      
C     adjust positions if too close to core
      call adjustDxnDyn(dxn,dyn,materials(mnum)%rcore)
      
C     get displacements from analytical solution
      if (imageoption) then
          call getImageDispAtPointSub(dxn,dyn,bsgn,mnum,uxn,uyn)
      else    
          call getDispAtPointSub(dxn,dyn,bsgn,mnum,uxn,uyn)
      end if    
      
C     rotate displacements back (by negative theta)
      call rotateVec2d(cost,-sint,uxn,uyn,uxo,uyo)
      
      disp(1) = uxo
      disp(2) = uyo
      
      end function getDispAtPointHelper
************************************************************************
      function getStressAtPoint(posn,dislpos,theta,bsgn,mnum)
     &                             result(stress)

C     Function: getStressAtPoint

C     Inputs: posn - 2 by 1 vector of position of point at which stress is sought
C             dislpos - 2 by 1 vector of position of dislocation
C             theta - angle of dislocation
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             bsgn - sign of dislocation (+1 or -1)
                        
C     Outputs: stress - 3 by 1 vector of (2d) stresses of point,
C                     [s11, s22, tau] (in original coordinate system,
C                     *not* dislocation coordinate system)
     
C     Purpose: Compute stresses at point due to a single edge
C              dislocation. Performs appropriate rotations so that
C              elasticity expressions can be used. Also adjusts for points
C              that are too close to the dislocation core.

C     Notes: 'n' in dxn, etc. refers to the 'new' (dislocation)
C            coordinate system; 'o' refers to the old (unrotated) system.

      implicit none
      
C     input variables
      real(dp) :: posn(2), dislpos(2)
      real(dp) :: theta
      integer :: bsgn
      integer :: mnum
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      real(dp) :: dxo, dyo, dxn, dyn
      real(dp) :: cost, sint
      real(dp) :: s11o, s22o, tauo
      real(dp) :: s11n, s22n, taun
      
      cost = cos(theta)
      sint = sin(theta)
      dxo = posn(1) - dislpos(1)
      dyo = posn(2) - dislpos(2)
      
C     rotate coordinates into new coordinate system
      call rotateVec2d(cost,sint,dxo,dyo,dxn,dyn)
      
C     adjust positions if too close to core
      call adjustDxnDyn(dxn,dyn,materials(mnum)%rcore)
      
C     get stresses
      call getStressAtPointSub(dxn,dyn,bsgn,mnum,s11n,s22n,taun)
      
C     rotate stresses back (by negative theta)
      call rotateStress2d(cost,-sint,s11n,s22n,taun,s11o,s22o,tauo)
      
      stress(1) = s11o
      stress(2) = s22o
      stress(3) = tauo
      
      end function getStressAtPoint
************************************************************************
      subroutine adjustDxnDyn(dxn,dyn,rcore)
     
C     Subroutine: adjustDxnDyn

C     Inputs: dxn, dyn - relative coordinates of point w.r.t. dislocation,
C                        in the dislocation coordinate system
C             rcore - radius of dislocation core
                        
C     Outputs: dxn, dyn - adjusted relative coordinates (moved outside
C                         of core)
     
C     Purpose: Adjust coordinates of point so it lies outside of disl.
C              core., to avoid spuriously large displacements/stresses

C     Notes/TODO: Should these fields simply be excluded instead?
     
      implicit none
      
C     notes
      
C     input variables
      real(dp) :: rcore
      
C     in/out variables
      real(dp) :: dxn, dyn
      
C     local variables
      real(dp) :: rsq
      real(dp) :: rcoresq, rcorefac, rfac
      
C     core radius
C     (Points too close to core are treated somewhat inconsistently in
C     original DD code --- compare routines dissig.f and dislp.f)
      rcoresq = rcore**2

      rsq = dxn*dxn + dyn*dyn      
C     If point and dislocation are basically coincident
      if (rsq < tolconst*rcoresq) then
          rcorefac = rcore/sqrt(2.0_dp)
          dxn = sign(rcorefac,dxn)
          dyn = sign(rcorefac,dyn)
C     For other points inside core radius
      else if (rsq < rcoresq) then
          rfac = rcore/sqrt(rsq)
          dxn = dxn*rfac
          dyn = dyn*rfac
      end if
      
      end subroutine adjustDxnDyn
************************************************************************
      subroutine getDispAtPointSub(dxn,dyn,bsgn,mnum,uxn,uyn)

C     Subroutine: getDispAtPointSub

C     Inputs: dxn, dyn - (adjusted) relative coordinates of point w.r.t.
C                         dislocation, in the dislocation coordinate system
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             bsgn - sign of dislocation (+1 or -1)
                        
C     Outputs: uxn, uyn - displacements of point due to (edge) dislocation,
C                         in the dislocation coordinate system
     
C     Purpose: Compute displacements at point due to a single edge
C              dislocation, using standard expressions from 2D elasticity
C              (see vdG and Needleman, MSMSE, 1995)

C     Notes: 'n' in dxn, etc. refers to the 'new' (dislocation)
C            coordinate system; 'o' refers to the old (unrotated) system.
C            Also, see long note below about error in original DD code.
      
      implicit none
      
C     input variables
      real(dp) :: dxn, dyn
      integer :: bsgn
      integer :: mnum
      
C     output variables
      real(dp) :: uxn, uyn
      
C     local variables
      real(dp) :: nu
      real(dp) :: burgers
      real(dp) :: rsq, invrsq
      real(dp) :: prefac, phi
      
C     Initialize
      nu = materials(mnum)%nu
      burgers = materials(mnum)%burgers
      
C     Auxiliary constants
      rsq = dxn*dxn + dyn*dyn
      invrsq = 1.0_dp/rsq
      prefac = 0.25_dp*bsgn*burgers/(piconst*(1.0_dp-nu))
      phi = atan2(dyn,dxn)
      
C     Displacement fields

C     Long note: The ux displacement field in vdG and Needleman are not exactly
C     correct. Using the bare arctan(dx/dy) leads to two discontinuities,
C     one at x2 = 0, x1 > 0, the other at x2 = 0, x1 < 0. In fact, there
C     should only be one discontinuity; the former one is spurious.
C     This can also be seen in the polar coordinate solution of Bower (
C     Section 5.3.4), although his sign convention appears to be different.
C     Using atan2 with two arguments, and the Hirth and Lothe solution for ux
C     (which is the VDG and Needleman solution, modified using the identity
C     tan(x) + tan(1/x) = const.), gives only one discontinuity,
C     and in the correct location. This appeared to have been fixed in
C     the old CADD code, but was not correct in the DD code used by O'Day,
C     Cleveringa, etc.

      uxn = prefac*(dxn*dyn*invrsq + 2.0_dp*(1.0_dp-nu)*phi)
      uyn = prefac*(dyn*dyn*invrsq - (0.5_dp-nu)*log(rsq/burgers**2))
      
      end subroutine getDispAtPointSub
************************************************************************
      subroutine getImageDispAtPointSub(dxn,dyn,bsgn,mnum,uxn,uyn)

C     Subroutine: getImageDispAtPointSub

C     Inputs: dxn, dyn - (adjusted) relative coordinates of point w.r.t.
C                         dislocation, in the dislocation coordinate system
C             mnum - number of material that image dislocation used to reside in
C             bsgn - sign of dislocation (+1 or -1)
                        
C     Outputs: uxn, uyn - displacements of point due to (edge) image dislocation,
C                         in the dislocation coordinate system
     
C     Purpose: Compute displacements at point due to a plastic slip field
C              from "image" dislocation

C     Notes: 'n' in dxn, etc. refers to the 'new' (dislocation)
C            coordinate system; 'o' refers to the old (unrotated) system.
C            Also, see long note below about error in original DD code.
      
      implicit none
      
C     input variables
      real(dp) :: dxn, dyn
      integer :: bsgn
      integer :: mnum
      
C     output variables
      real(dp) :: uxn, uyn
      
C     local variables
      real(dp) :: burgers

      burgers = materials(mnum)%burgers

C     Long note: I believe this is implemented incorrectly in the original
C     DD code. There, the plastic displacements from an image dislocation (i.e.
C     a dislocation that left the material) were b/4 and -b/4 for points
C     above and below the dislocation, respectively, regardless of whether
C     the point was to the right or left of the dislocation. This doesn't seem to
C     make much sense. And, in fact, they appear to be cancelling errors...
C     if one considers a dipole, where one dislocation exits the right side
C     and the other the left, only the right one contributes plastic slip
C     to the body. So, the slip is b/2 on top and -b/2 on bottom. The original
C     DD code gets this result, but only by adding b/4 from both dislocations.
C     This procedure results in errors when only one of the dislocations (not both)
C     have left the body.
      
      if (dxn < 0) then
          uxn = bsgn*burgers*sign(0.5_dp,dyn)
      else
          uxn = 0.0_dp
      end if
      uyn = 0.0_dp
      
      end subroutine getImageDispAtPointSub
************************************************************************
      subroutine getStressAtPointSub(dxn,dyn,bsgn,mnum,
     &                               sig11n,sig22n,taun)

C     Subroutine: getStressAtPointSub

C     Inputs: dxn, dyn - (adjusted) relative coordinates of point w.r.t.
C                         dislocation, in the dislocation coordinate system
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             bsgn - sign of dislocation (+1 or -1)
                        
C     Outputs: sig11n, sig22n, taun - stress at point due to (edge)
C                     dislocation, in the dislocation coordinate system
     
C     Purpose: Compute stresses at point due to a single edge
C              dislocation, using standard expressions from 2D elasticity
C              (see vdG and Needleman, MSMSE, 1995)

C     Notes: 'n' in dxn, etc. refers to the 'new' (dislocation)
C            coordinate system; 'o' refers to the old (unrotated) system.

      implicit none

C     input variables
      real(dp) :: dxn, dyn
      integer :: bsgn
      integer :: mnum
      
C     output variables
      real(dp) :: sig11n, sig22n, taun
      
C     local variables
      real(dp) :: mu, nu
      real(dp) :: burgers
      real(dp) :: dxnsq, dynsq, rsq, dxynsqdiff
      real(dp) :: prefac

C     Initialize
      mu = materials(mnum)%mu
      nu = materials(mnum)%nu
      burgers = materials(mnum)%burgers

C     Auxiliary constants
      dxnsq = dxn*dxn
      dynsq = dyn*dyn
      rsq = dxnsq + dynsq
      dxynsqdiff = dxnsq - dynsq
      prefac = 0.5_dp*bsgn*burgers*mu/(piconst*(1.0_dp-nu)*rsq**2)
          
C     Stress fields
      sig11n = -prefac*dyn*(3.0_dp*dxnsq + dynsq)
      sig22n = prefac*dyn*dxynsqdiff
      taun = prefac*dxn*dxynsqdiff
      
      end subroutine getStressAtPointSub
************************************************************************

      end module