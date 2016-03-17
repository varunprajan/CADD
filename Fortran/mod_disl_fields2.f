      module mod_disl_fields2
      
C     Purpose: Computes fields (tilde fields, in van der Giessen and
C     Needleman's terminology) created by all dislocations acting in an
C     *infinite* medium. Also computes fields from "escaped" dislocations,
C     which have left the mesh, and contribute only plastic displacement
C     (slip) and no stress.
      
C     Possible extensions:
C     1) Currently, expressions for fields assume
C     isotropy. Although our toy potentials are isotropic, this is not
C     the case for many materials.

C     2) Calculates fields directly, which has a total cost of O(N^2) (?)
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
      use mod_math, only: TOLCONST, rotateVec2d, rotateStress2d
      use mod_disl_try, only: disl, sources
      use mod_disl_ghost, only: ghostdisl
      use mod_slip_sys, only: slipsys
      use mod_materials, only: materials
      use mod_fe_elements, only: fematerials
      use mod_disl_escaped, only: getEscapedDispAtPointAll
      implicit none

      private
      public :: getTildeDispAtPointAll, getTildeStressatPointAll,
     &            getPKTildeStressAll, getDispAtPointSub,
     &            getStressAtPointSub, adjustDxnDyn,
     &            getDispAtPoint, getStressAtPoint,
     &          getGhostStressAtPointAll, getRealStressAtPointAll,
     &          getLatentStressAtPointAll, getLatentStressOnSourceAll,
     &          getTildeStressOnSourceAll
      
      contains
************************************************************************
      function getTildeDispAtPointAll(posn,mnumfe) result(disp)

C     Inputs: posn --- position at which tilde displacement field is sought
C                      (vector, length 2)
C             mnumfe --- fe material in which point lies

C     Outputs: disp --- tilde displacement at point (vector, length 2)

C     Purpose: Get displacement field from all dislocations (real, escaped, latent, gost)
C              *associated with fe material*
C              at point of interest (to get total displacement, FE contribution is also required)
      
      implicit none
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      
C     output variables
      real(dp) :: disp(2)
      
      disp = 0.0_dp
      disp = disp + getRealDispAtPointAll(posn,mnumfe)
      disp = disp + getLatentDispAtPointAll(posn,mnumfe)
      disp = disp + getGhostDispAtPointAll(posn,mnumfe)
      disp = disp + getEscapedDispAtPointAll(posn,mnumfe)
      
      end function getTildeDispAtPointAll
************************************************************************
      function getTildeStressAtPointAll(posn,mnumfe) result(stress)

C     Inputs: posn --- position at which tilde displacement field is sought
C                      (vector, length 2)
C             mnumfe --- fe material in which point lies

C     Outputs: stress --- tilde stress at point (vector, length 3 --- sxx, syy, sxy)

C     Purpose: Get stress field from all dislocations (real, latent)
C              *associated with fe material*
C              at point of interest (to get total stress, FE contribution is also required)
      
      implicit none
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      
C     output variables
      real(dp) :: stress(3)
      
      stress = 0.0_dp
      stress = stress + getRealStressAtPointAll(posn,mnumfe)
      stress = stress + getLatentStressAtPointAll(posn,mnumfe)
      stress = stress + getGhostStressAtPointAll(posn,mnumfe)
      
      end function getTildeStressAtPointAll
************************************************************************
      function getPKTildeStressAll(dislnum,mnumfe) result(stress)

C     Inputs: dislnum --- number of dislocation for which Peach-Koehler stress is sought
C             mnumfe --- fe material in which point lies

C     Outputs: stress --- Peach-Koehler stress at point (vector, length 3 --- sxx, syy, sxy)

C     Purpose: Get contribution of all dislocations (real, latent) 
C              *associated with fe material*
C              to Peach-Koehler stress (to get total PK stress, FE contribution is also required)

C     Notes: See expression in parentheses in Equation 12 in vdG and Needleman, MSMSE, 1995.
      
      implicit none
      
C     input variables
      integer :: dislnum
      integer :: mnumfe
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      real(dp) :: posn(2)
      
      posn = disl(mnumfe)%list(dislnum)%posn
      stress = 0.0_dp
      stress = stress + getRealPKStressAll(dislnum,mnumfe)
      stress = stress + getLatentStressAtPointAll(posn,mnumfe)
C     (no contributions from ghost or escaped dislocations, since they're not real)
      
      end function getPKTildeStressAll
************************************************************************
      function getTildeStressOnSourceAll(sourcenum,mnumfe)
     &                                                   result(stress)

C     Inputs: sourcenum --- number of source for which stress is sought
C             mnumfe --- fe material in which point lies

C     Outputs: stress --- tilde stress at point (vector, length 3 --- sxx, syy, sxy)

C     Purpose: Get contribution of all dislocations (real, latent) 
C              *associated with fe material*
C              to stress on source (to get total source stress, FE contribution is also required)

C     Notes: See expression in parentheses in Equation 12 in vdG and Needleman, MSMSE, 1995.
      
      implicit none
      
C     input variables
      integer :: sourcenum
      integer :: mnumfe
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      real(dp) :: posn(2)
      
      posn = sources(mnumfe)%list(sourcenum)%posn
      stress = 0.0_dp
      stress = stress + getRealStressAtPointAll(posn,mnumfe)
      stress = stress + getLatentStressOnSourceAll(sourcenum,mnumfe)
      stress = stress + getGhostStressAtPointAll(posn,mnumfe)
      
      end function getTildeStressOnSourceAll
************************************************************************
      function getRealDispAtPointAll(posn,mnumfe) result(disp)

C     Inputs: posn - 2 by 1 vector of position of point at which displacement is sought
C             mnumfe --- fe material in which point lies 

C     Outputs: disp --- displacement at point (vector, 2 by 1)

C     Purpose: Get displacement at point from *all* real dislocations in material
      
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
      integer :: isys
      real(dp) :: cost, sint
      integer :: bsgn, bcut
      real(dp) :: dispnew(2)

      disp = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, disl(mnumfe)%ndisl
C     an actual dislocation
      if (disl(mnumfe)%list(i)%active) then
          dislpos = disl(mnumfe)%list(i)%posn
          isys = disl(mnumfe)%list(i)%slipsys
          cost = slipsys(mnumfe)%trig(1,isys)
          sint = slipsys(mnumfe)%trig(2,isys)
          bsgn = disl(mnumfe)%list(i)%sgn
          bcut = disl(mnumfe)%list(i)%cut
          dispnew = getDispAtPoint(posn,dislpos,
     &                             cost,sint,bsgn,bcut,mnum,1.0_dp) ! no fudge
          disp = disp + dispnew
      end if    
      end do
      
      end function getRealDispAtPointAll
************************************************************************
      function getRealStressAtPointAll(posn,mnumfe) result(stress)

C     Inputs: posn - 2 by 1 vector of position of point at which stress is sought
C             mnumfe --- fe material in which point lies 

C     Outputs: stress --- stress at point (vector, 3 by 1)
 
C     Purpose: Get stress at point from *all* real dislocations in material
      
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
      integer :: isys
      real(dp) :: cost, sint
      integer :: bsgn
      real(dp) :: stressnew(3)
           
      stress = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, disl(mnumfe)%ndisl
C     an actual dislocation
      if (disl(mnumfe)%list(i)%active) then
          dislpos = disl(mnumfe)%list(i)%posn
          isys = disl(mnumfe)%list(i)%slipsys
          cost = slipsys(mnumfe)%trig(1,isys)
          sint = slipsys(mnumfe)%trig(2,isys)
          bsgn = disl(mnumfe)%list(i)%sgn
          stressnew = getStressAtPoint(posn,dislpos,cost,sint,bsgn,
     &                                 mnum,1.0_dp) ! no fudge
          stress = stress + stressnew
      end if    
      end do
      
      end function getRealStressAtPointAll
************************************************************************
      function getRealPKStressAll(dislnum,mnumfe) result(stress)

C     Inputs: dislnum --- number of dislocation within fe material
C             mnumfe --- fe material in which point lies 

C     Outputs: stress --- stress at point (vector, 3 by 1)

C     Purpose: Get PK stress at point from *all* real dislocations in material
      
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
      integer :: isys
      real(dp) :: cost, sint
      integer :: bsgn
      real(dp) :: stressnew(3)
      
      posn = disl(mnumfe)%list(dislnum)%posn
      stress = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, disl(mnumfe)%ndisl
C     an actual dislocation
      if (disl(mnumfe)%list(i)%active) then
C         no self force
          if (i /= dislnum) then
              dislpos = disl(mnumfe)%list(i)%posn
              isys = disl(mnumfe)%list(i)%slipsys
              cost = slipsys(mnumfe)%trig(1,isys)
              sint = slipsys(mnumfe)%trig(2,isys)
              bsgn = disl(mnumfe)%list(i)%sgn
              stressnew = getStressAtPoint(posn,dislpos,cost,sint,bsgn,
     &                                                     mnum,1.0_dp) ! no fudge
              stress = stress + stressnew
          end if
      end if
      end do
      
      end function getRealPKStressAll
************************************************************************
      function getGhostDispAtPointAll(posn,mnumfe) result(disp)

C     Inputs: posn - 2 by 1 vector of position of point at which displacement is sought
C             mnumfe --- fe material in which point lies 

C     Outputs: disp --- displacement at point (vector, 2 by 1)

C     Purpose: Get displacement at point from *all* ghost dislocations associated with material
      
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
      integer :: isys
      real(dp) :: cost, sint
      integer :: bsgn, bcut
      real(dp) :: dispnew(2)

      disp = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, ghostdisl(mnumfe)%nghostdisl
          dislpos = ghostdisl(mnumfe)%list(i)%posn
          isys = ghostdisl(mnumfe)%list(i)%slipsys
          cost = slipsys(mnumfe)%trig(1,isys)
          sint = slipsys(mnumfe)%trig(2,isys)
          bsgn = ghostdisl(mnumfe)%list(i)%sgn
          bcut = ghostdisl(mnumfe)%list(i)%cut
          dispnew = getDispAtPoint(posn,dislpos,
     &                             cost,sint,bsgn,bcut,mnum,1.0_dp) ! no fudge
          disp = disp + dispnew    
      end do
      
      end function getGhostDispAtPointAll
************************************************************************
      function getGhostStressAtPointAll(posn,mnumfe) result(stress)

C     Inputs: posn - 2 by 1 vector of position of point at which stress is sought
C             mnumfe --- fe material in which point lies 

C     Outputs: stress --- stress at point (vector, 3 by 1)
 
C     Purpose: Get stress at point from *all* ghost dislocations associated with material
      
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
      integer :: isys
      real(dp) :: cost, sint
      integer :: bsgn
      real(dp) :: stressnew(3)
           
      stress = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, ghostdisl(mnumfe)%nghostdisl
          dislpos = ghostdisl(mnumfe)%list(i)%posn
          isys = ghostdisl(mnumfe)%list(i)%slipsys
          cost = slipsys(mnumfe)%trig(1,isys)
          sint = slipsys(mnumfe)%trig(2,isys)
          bsgn = ghostdisl(mnumfe)%list(i)%sgn
          stressnew = getStressAtPoint(posn,dislpos,cost,sint,bsgn,
     &                                 mnum,1.0_dp) ! no fudge
          stress = stress + stressnew   
      end do
      
      end function getGhostStressAtPointAll
************************************************************************
      function getLatentDispAtPointAll(posn,mnumfe) result(disp)

C     Inputs: posn - 2 by 1 vector of position of point at which displacement is sought
C             mnumfe --- fe material in which point lies 

C     Outputs: disp --- displacement at point (vector, 2 by 1)

C     Purpose: Get displacement at point from *all* latent dislocations in material (from sources)
      
      implicit none
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      
C     output variables
      real(dp) :: disp(2)

C     local variables
      integer :: i
      integer :: mnum
      real(dp) :: time
      real(dp) :: sourcepos(2), dpos(2)
      integer :: isys
      real(dp) :: cost, sint
      real(dp) :: tnuc, lnuc
      real(dp) :: tauprev
      real(dp) :: bfudge
      real(dp) :: dispnew(2), dispnew2(2)

      disp = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, size(sources(mnumfe)%list)
      time = sources(mnumfe)%list(i)%time
      if (time > TOLCONST) then ! source is active
          sourcepos = sources(mnumfe)%list(i)%posn
          isys = sources(mnumfe)%list(i)%slipsys
          cost = slipsys(mnumfe)%trig(1,isys)
          sint = slipsys(mnumfe)%trig(2,isys)
          tnuc = sources(mnumfe)%list(i)%tnuc
          bfudge = time/tnuc
          lnuc = sources(mnumfe)%list(i)%lnuc
          tauprev = sources(mnumfe)%list(i)%tauprev
          dpos = (0.5_dp*lnuc*bfudge)*[cost,sint]
          if (tauprev < 0.0_dp) then ! flip dipole if tau is negative
              dpos = -dpos
          end if
          dispnew = getDispAtPoint(posn,sourcepos+dpos,
     &                             cost,sint,+1,0,mnum,bfudge)
          dispnew2 = getDispAtPoint(posn,sourcepos-dpos,
     &                              cost,sint,-1,0,mnum,bfudge)
          disp = disp + dispnew + dispnew2
      end if
      end do
      
      end function getLatentDispAtPointAll
************************************************************************
      function getLatentStressAtPointAll(posn,mnumfe) result(stress)

C     Inputs: posn - 2 by 1 vector of position of point at which stress is sought
C             mnumfe --- fe material in which point lies 

C     Outputs: stress --- stress at point (vector, 3 by 1)

C     Purpose: Get stress at point from *all* latent dislocations in materials (from sources)
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      integer :: i
      integer :: mnum
      real(dp) :: time
      real(dp) :: sourcepos(2), dpos(2)
      integer :: isys
      real(dp) :: cost, sint
      real(dp) :: tnuc, lnuc
      real(dp) :: tauprev
      real(dp) :: bfudge
      real(dp) :: stressnew(3), stressnew2(3)
      
      stress = 0.0_dp
      mnum = fematerials%list(mnumfe)
      do i = 1, size(sources(mnumfe)%list)
      time = sources(mnumfe)%list(i)%time
      if (time > TOLCONST) then ! source is active
          sourcepos = sources(mnumfe)%list(i)%posn
          isys = sources(mnumfe)%list(i)%slipsys
          cost = slipsys(mnumfe)%trig(1,isys)
          sint = slipsys(mnumfe)%trig(2,isys)
          tnuc = sources(mnumfe)%list(i)%tnuc
          bfudge = time/tnuc
          lnuc = sources(mnumfe)%list(i)%lnuc
          tauprev = sources(mnumfe)%list(i)%tauprev
          dpos = (0.5_dp*lnuc*bfudge)*[cost,sint]
          if (tauprev < 0.0_dp) then ! flip dipole if tau is negative
              dpos = -dpos
          end if
          stressnew = getStressAtPoint(posn,sourcepos+dpos,
     &                             cost,sint,+1,mnum,bfudge)
          stressnew2 = getStressAtPoint(posn,sourcepos-dpos,
     &                              cost,sint,-1,mnum,bfudge)
          stress = stress + stressnew + stressnew2
      end if
      end do
      
      end function getLatentStressAtPointAll
************************************************************************
      function getLatentStressOnSourceAll(sourcenum,mnumfe)
     &                                                   result(stress)

C     Inputs: sourcenum --- number of source
C             mnumfe --- fe material in which point lies 

C     Outputs: stress --- stress at point (vector, 3 by 1)

C     Purpose: Get stress on source at point from *all* latent dislocations in materials (from sources)
C     *except* latent dislocations from source itself (basically, no "self" stress)
      
      implicit none
      
C     input variables
      integer :: sourcenum
      integer :: mnumfe
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      integer :: i
      integer :: mnum
      real(dp) :: time
      real(dp) :: posn(2), sourcepos(2), dpos(2)
      integer :: isys
      real(dp) :: cost, sint
      real(dp) :: tnuc, lnuc
      real(dp) :: tauprev
      real(dp) :: bfudge
      real(dp) :: stressnew(3), stressnew2(3)
      
      stress = 0.0_dp
      posn = sources(mnumfe)%list(sourcenum)%posn
      mnum = fematerials%list(mnumfe)
      do i = 1, size(sources(mnumfe)%list)
          if (i /= sourcenum) then ! no "self"-stress
          time = sources(mnumfe)%list(i)%time
          if (time > TOLCONST) then ! source is active
              sourcepos = sources(mnumfe)%list(i)%posn
              isys = sources(mnumfe)%list(i)%slipsys
              cost = slipsys(mnumfe)%trig(1,isys)
              sint = slipsys(mnumfe)%trig(2,isys)
              tnuc = sources(mnumfe)%list(i)%tnuc
              bfudge = time/tnuc
              lnuc = sources(mnumfe)%list(i)%lnuc
              tauprev = sources(mnumfe)%list(i)%tauprev
              dpos = (0.5_dp*lnuc*bfudge)*[cost,sint]
              if (tauprev < 0.0_dp) then ! flip dipole if tau is negative
                  dpos = -dpos
              end if
              stressnew = getStressAtPoint(posn,sourcepos+dpos,
     &                                 cost,sint,+1,mnum,bfudge)
              stressnew2 = getStressAtPoint(posn,sourcepos-dpos,
     &                                  cost,sint,-1,mnum,bfudge)
              stress = stress + stressnew + stressnew2
          end if
          end if
      end do
      
      end function getLatentStressOnSourceAll
************************************************************************
      function getDispAtPoint(posn,dislpos,cost,sint,bsgn,bcut,
     &                       mnum,bfudge) result(disp)

C     Inputs: posn - 2 by 1 vector of position of point at which displacement is sought
C             dislpos - 2 by 1 vector of position of dislocation
C             cost, sint - trig for angle of dislocation
C             bsgn - sign of dislocation (+1 or -1)
C             bcut - branch cut of dislocation (0 if to the left, 1 if to the right)
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             bfudge - fudge of burgers for latent dislocations (equal to 1 for non-latent (real) dislocations)
                        
C     Outputs: disp - 2 by 1 vector of displacements of point (in original
C                     coordinate system, *not* dislocation coordinate system)
     
C     Purpose: Compute displacements at point due to a single 
C              (real/latent/ghost) edge dislocation.
C              Performs appropriate rotations so that
C              elasticity expressions can be used. Also adjusts for points
C              that are too close to the dislocation core.

C     Notes: 'n' in dxn, etc. refers to the 'new' (dislocation)
C            coordinate system; 'o' refers to the old (unrotated) system.
      
      implicit none
      
C     input variables
      real(dp) :: posn(2), dislpos(2)
      real(dp) :: cost, sint
      integer :: mnum
      integer :: bsgn, bcut
      real(dp) :: bfudge
      
C     output variables
      real(dp) :: disp(2)
      
C     local variables
      real(dp) :: dxo, dyo, dxn, dyn
      real(dp) :: uxn, uyn
      
      dxo = posn(1) - dislpos(1)
      dyo = posn(2) - dislpos(2)
      
C     adjust for branch cut (rotate everything by 180 degrees, reverse sign)
      if (bcut /= 0) then
          cost = -cost
          sint = -sint
          bsgn = -bsgn
      end if 
      
C     rotate coordinates into new coordinate system
      call rotateVec2d(cost,sint,dxo,dyo,dxn,dyn)
      
C     adjust positions if too close to core
      call adjustDxnDyn(dxn,dyn,materials(mnum)%rcoresq)
      
C     get displacements from analytical solution
      call getDispAtPointSub(dxn,dyn,bsgn,mnum,bfudge,uxn,uyn)  
      
C     rotate displacements back (by negative theta)
      call rotateVec2d(cost,-sint,uxn,uyn,disp(1),disp(2))
      
      end function getDispAtPoint
************************************************************************
      function getStressAtPoint(posn,dislpos,cost,sint,bsgn,mnum,bfudge)
     &                             result(stress)

C     Inputs: posn - 2 by 1 vector of position of point at which stress is sought
C             dislpos - 2 by 1 vector of position of dislocation
C             cost, sint - trig for angle of dislocation
C             bsgn - sign of dislocation (+1 or -1)
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             bfudge - fudge of burgers for latent dislocations (equal to 1 for non-latent (real) dislocations)
                        
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
      real(dp) :: cost, sint
      integer :: bsgn
      integer :: mnum
      real(dp) :: bfudge
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      real(dp) :: dxo, dyo, dxn, dyn
      real(dp) :: s11n, s22n, taun
      
      dxo = posn(1) - dislpos(1)
      dyo = posn(2) - dislpos(2)
      
C     cut is irrelevant for stresses
      
C     rotate coordinates into new coordinate system
      call rotateVec2d(cost,sint,dxo,dyo,dxn,dyn)
      
C     adjust positions if too close to core
      call adjustDxnDyn(dxn,dyn,materials(mnum)%rcoresq)
      
C     get stresses
      call getStressAtPointSub(dxn,dyn,bsgn,mnum,bfudge,s11n,s22n,taun)
      
C     rotate stresses back (by negative theta)
      call rotateStress2d(cost,-sint,s11n,s22n,taun,stress(1),
     &                                         stress(2),stress(3))
      
      end function getStressAtPoint
************************************************************************
      subroutine adjustDxnDyn(dxn,dyn,rcoresq)

C     Inputs: dxn, dyn - relative coordinates of point w.r.t. dislocation,
C                        in the dislocation coordinate system
C             rcoresq - radius of dislocation core squared
                        
C     Outputs: dxn, dyn - adjusted relative coordinates (moved outside
C                         of core)
     
C     Purpose: Adjust coordinates of point so it lies outside of disl.
C              core., to avoid spuriously large displacements/stresses

C     Notes/TODO: Should these fields simply be excluded instead?
     
      implicit none
      
C     input variables
      real(dp) :: rcoresq
      
C     in/out variables
      real(dp) :: dxn, dyn
      
C     local variables
      real(dp) :: rsq
      real(dp) :: rcorefac, rfac
      
C     (Points too close to core are treated somewhat inconsistently in
C     original DD code --- compare routines dissig.f and dislp.f)

      rsq = dxn*dxn + dyn*dyn      
C     
      if (rsq < rcoresq) then
          if (rsq < TOLCONST*rcoresq) then ! if point and dislocation are basically coincident
              rcorefac = sqrt(rcoresq/2.0_dp)
              dxn = sign(rcorefac,dxn)
              dyn = sign(rcorefac,dyn)
          else ! for other points inside core radius
              rfac = sqrt(rcoresq/rsq)
              dxn = dxn*rfac
              dyn = dyn*rfac
          end if    
      end if
      
      end subroutine adjustDxnDyn
************************************************************************
      subroutine getDispAtPointSub(dxn,dyn,bsgn,mnum,bfudge,uxn,uyn)

C     Inputs: dxn, dyn - (adjusted) relative coordinates of point w.r.t.
C                         dislocation, in the dislocation coordinate system
C             bsgn - sign of dislocation (+1 or -1)
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             bfudge - fudge of burgers for latent dislocations (equal to 1 for non-latent (real) dislocations)
                        
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
      real(dp) :: bfudge
      
C     output variables
      real(dp) :: uxn, uyn
      
C     local variables
      real(dp) :: nu
      real(dp) :: burgers, bsgndp
      real(dp) :: rsq, invrsq
      real(dp) :: prefac, phi
      
C     Initialize
      nu = materials(mnum)%nu
      burgers = materials(mnum)%burgers
      prefac = materials(mnum)%dispprefac
      
C     Auxiliary constants
      rsq = dxn*dxn + dyn*dyn
      invrsq = 1.0_dp/rsq
      bsgndp = real(bsgn,dp)
      prefac = sign(bfudge*prefac,bsgndp)
      phi = atan2(dyn,dxn)
      
C     Displacement fields

C     Long note: The ux displacement field in vdG and Needleman is not exactly
C     correct. Using the bare arctan(dx/dy) leads to two discontinuities,
C     one at x2 = 0, x1 > 0, the other at x2 = 0, x1 < 0. In fact, there
C     should only be one discontinuity; the former one is spurious.
C     This can also be seen in the polar coordinate solution of Bower (
C     Section 5.3.4), although his sign convention appears to be different.
C     Using atan2 with two arguments, and the Hirth and Lothe solution for ux
C     (which is the VDG and Needleman solution, modified using the identity
C     tan(x) + tan(1/x) = const.), gives only one discontinuity,
C     and in the correct location. This appears to be wrong in the old
C     CADD code and the old DD code.

      uxn = prefac*(dxn*dyn*invrsq + 2.0_dp*(1.0_dp-nu)*phi)
      uyn = prefac*(dyn*dyn*invrsq - (0.5_dp-nu)*log(rsq/burgers**2))
      
      end subroutine getDispAtPointSub
************************************************************************
      subroutine getStressAtPointSub(dxn,dyn,bsgn,mnum,bfudge,
     &                               sig11n,sig22n,taun)

C     Inputs: dxn, dyn - (adjusted) relative coordinates of point w.r.t.
C                         dislocation, in the dislocation coordinate system
C             bsgn - sign of dislocation (+1 or -1)
C             mnum - number of material that point and (edge) dislocation
C                    reside within
C             bfudge - fudge of burgers for latent dislocations (equal to 1 for non-latent (real) dislocations)
                        
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
      real(dp) :: bfudge
      
C     output variables
      real(dp) :: sig11n, sig22n, taun
      
C     local variables
      real(dp) :: bsgndp
      real(dp) :: dxnsq, dynsq, rsq, dxynsqdiff
      real(dp) :: prefac

C     Initialize
      prefac = materials(mnum)%stressprefac
      
C     Auxiliary constants
      dxnsq = dxn*dxn
      dynsq = dyn*dyn
      rsq = dxnsq + dynsq
      dxynsqdiff = dxnsq - dynsq
      bsgndp = real(bsgn,dp)
      prefac = sign(bfudge*prefac,bsgndp)/rsq**2
          
C     Stress fields
      sig11n = -prefac*dyn*(3.0_dp*dxnsq + dynsq)
      sig22n = prefac*dyn*dxynsqdiff
      taun = prefac*dxn*dxynsqdiff
      
      end subroutine getStressAtPointSub
************************************************************************

      end module