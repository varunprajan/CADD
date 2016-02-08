      module mod_pk_update
      
C     Purpose: Update dislocation positions/elements using Peach-Koehler
C     force, mobility law, and stresses from DD and FE fields (tilde and
C     hat stresses, respectively)

C     Possible extensions: Use Srinath's gradient correction to velocity (probably
C     requires multipole; otherwise, it would be horribly expensive)
      
      use mod_types, only: dp
      use mod_disl, only: disl, assignDislocation, addDislocation,
     &                    addImageDislocation, deleteDislocation
      use mod_disl_fields2, only: getPKTildeStressAll
      use mod_fe_elements, only: nfematerials, fematerials, feelements,
     &                           interfaceedges
      use mod_fe_main_2d, only: getFEStrainAtPoint
      use mod_materials, only: materials
      use mod_mesh_find, only: findInAllWithGuess
      use mod_disl_detect_pass, only: findInterfaceIntersection,
     &                                passContinuumtoAtomistic
      implicit none
      
      private
      public :: updatePK
      
      contains
************************************************************************
      subroutine updatePK(dt)

C     Subroutine: updatePK

C     Inputs: dt --- time increment for update

C     Outputs: None

C     Purpose: Update positions of dislocations within all materials,
C              and determine the new elements that these dislocations
C              belong to
      
      implicit none
      
C     input variables
      real(dp) :: dt
      
C     local variables
      integer :: i
      
      do i = 1, nfematerials
          call updatePosPKOneMat(i,dt)
      end do
      
      end subroutine updatePK
************************************************************************
      subroutine updatePosPKOneMat(mnumfe,dt)

C     Subroutine: updatePosPKOneMat

C     Inputs: mnumfe --- material in which to update dislocation positions
C             dt --- time increment for update

C     Outputs: None

C     Purpose: Update positions of dislocations within a single material
C              using Peach-Koehler force and mobility law. Uses hat and
C              tilde stresses to determine Peach-Koehler force
      
      implicit none
      
C     input variables
      integer :: mnumfe
      real(dp) :: dt
      
C     local variables
      integer :: i
      integer :: mnum
      integer :: element
      real(dp) :: stress(3), stresshat(3)
      real(dp) :: r, s
      real(dp) :: theta
      real(dp) :: burgers, disldrag, bfac, bfacabs
      real(dp) :: dislpos(2), dislposnew(2)
      real(dp) :: velocity(2)
      real(dp) :: C(3,3)
      integer :: nelnodes, neldof
      integer :: bsgn

      mnum = fematerials%list(mnumfe)
      burgers = materials(mnum)%burgers
      disldrag = materials(mnum)%disldrag
      bfacabs = burgers/disldrag
      C = materials(mnum)%elconst
      nelnodes = feelements(mnumfe)%nelnodes
      neldof = feelements(mnumfe)%neldof
      
C     loop over real dislocations
      do i = 1, disl(mnumfe)%ndisl
      element = disl(mnumfe)%element(i)
      if (element/=0) then
          r = disl(mnumfe)%localpos(1,i)
          s = disl(mnumfe)%localpos(2,i)
          dislpos = disl(mnumfe)%posn(1:2,i)
          theta = disl(mnumfe)%posn(3,i)
          bsgn = disl(mnumfe)%sgn(i)
          
          stress = getPKTildeStressAll(i,mnumfe)
          stresshat = matmul(C,
     &           getFEStrainAtPoint(mnumfe,element,nelnodes,neldof,r,s))
          stress = stress + stresshat
          bfac = bfacabs*bsgn
          velocity = getPKVelocity(stress,theta,bfac)
          dislposnew = dislpos + velocity*dt      

          call updateDislPos(dislpos,dislposnew,mnumfe,element,i,theta,
     &                       bsgn)
      end if    
      end do
      
      end subroutine updatePosPKOneMat
************************************************************************
      subroutine updateDislPos(dislposold,dislposnew,mnumfeold,
     &                         elementold,dislnum,theta,bsgn)
     
C     input variables
      real(dp) :: dislposold(2), dislposnew(2)
      integer :: mnumfeold, elementold
      integer :: dislnum
      real(dp) :: theta
      integer :: bsgn
     
C     local variables
      integer :: mnumfe, element
      real(dp) :: r, s
      logical :: badflip
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      real(dp) :: imagepos(2)

C     first, figure out if dislocation is still in mesh
      mnumfe = mnumfeold
      element = elementold
      call findInAllWithGuess(dislposnew(1),dislposnew(2),
     &                        mnumfe,element,r,s,badflip)
      if (badflip) then ! not in mesh
C         Two possibilities:
C         1) Crossed back to atomistic region
C         2) Left mesh, leaving a slip step (associated with "image" dislocation)
C         Check first possibility by seeing if it's crossed interface
C         between continuum and atomistic

          call findInterfaceIntersection(interfaceedges%array,.true., ! use undeformed positions
     &                        dislposold,dislposnew,pint,isint,edgenum)
          mnumfe = interfaceedges%array(3,edgenum)
          element = interfaceedges%array(4,edgenum)
          if (isint) then ! in atomistic region
              call passContinuumToAtomistic(dislposold,pint,bsgn,theta,
     &                                      dislnum,mnumfeold)
          else ! dislocation left mesh, add image dislocation at "infinity"
               ! (otherwise, displacements appear to be incorrect for angled bodies)
              imagepos = dislposold +
     &               (0.1_dp*huge(dislposold))*(dislposnew - dislposold)
              call addImageDislocation(imagepos(1),imagepos(2),
     &                                 theta,bsgn,mnumfeold)
          end if
          
      else ! still in mesh
          if (mnumfe/=mnumfeold) then ! moved to different material
              call addDislocation(mnumfe,element,dislposnew(1),
     &                            dislposnew(2),theta,bsgn,r,s)
              call deleteDislocation(mnumfeold,dislnum)
          else ! in the same material; simply update position
              call assignDislocation(mnumfe,element,dislposnew(1),
     &                            dislposnew(2),theta,bsgn,r,s,dislnum)
              end if
          end if
     
      end subroutine updateDislPos
************************************************************************
      function getPKVelocity(stress,theta,bfac) result(velocity)

C     Subroutine: getPKVelocity

C     Inputs: stress --- stress at point, vector of length 3 (sxx,syy,sxy)
C             theta --- orientation of dislocation (radians)
C             bfac --- equals bsgn*burgers/disldrag (speed = tau*bfac)

C     Outputs: velocity --- velocity of dislocation, vector of length 2 (vx, vy)

C     Purpose: Determine dislocation velocity using stress on dislocation,
C              using mobility law ("linear drag relation" --- see
C              vdG and Needleman, MSMSE, 1995)
      
      implicit none
      
C     input variables
      real(dp) :: stress(3)
      real(dp) :: theta
      real(dp) :: bfac
      
C     output variables
      real(dp) :: velocity(2)
      
C     local variables
      real(dp) :: mx, my
      real(dp) :: tx, ty
      real(dp) :: tau
      real(dp) :: speed
      
C     Equation 12 or 18, vdG and Needleman, MSMSE, 1995
      mx = cos(theta)
      my = sin(theta)
      tx = stress(1)*mx + stress(3)*my
      ty = stress(3)*mx + stress(2)*my
      tau = ty*mx - tx*my
      speed = tau*bfac
      velocity(1) = speed*mx
      velocity(2) = speed*my
      
      end function getPKVelocity
************************************************************************       
      end module