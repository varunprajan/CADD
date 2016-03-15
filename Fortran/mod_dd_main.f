      module mod_dd_main
      
C     Purpose: Collection of routines for running a single DD step:
C     updates dislocation positions, enforces obstacles, annihilates
C     opposite signed dislocations that cross, etc.

      use mod_types, only: dp
      use mod_disl_misc, only: dislmisc
      use mod_disl_try, only: disl, deleteDislocationSub, 
     & sortPlaneCheck, obstaclet, sortedplanedata,
     & deleteDislocationSub2, addDislocation, deleteDislocation,
     & zeroObstacles, obstacles, zeroDislDisp, sources,
     & addDislocation, sortDislPlanes, sourcet
      use mod_disl_escaped, only: addEscapedDislocation
      use mod_mesh_find, only: findInAllWithGuess
      use mod_disl_detect_pass, only: findInterfaceIntersectionDeformed,
     &   passContinuumtoAtomistic
      use mod_slip_sys, only: resolveStress, invResolveDisp, slipsys
      use mod_materials, only: materials
      use mod_fe_elements, only: fematerials, interfaceedges
      use mod_fe_main_2d, only: getFEStressAtPoint
      use mod_disl_fields2, only: getTildeStressOnSourceAll,
     &  getTildeStressAtPointAll
      use mod_math, only: findPointBetween, TOLCONST, sameSign
      use mod_dd_integrate, only: dispFromPK, assignDispFromPKOneMat
      implicit none
      
      private
      public :: runDDStep,
     & updateDislocations, enforceObstacles, isObstacleBetween,
     & getResolvedStressOnObstacle, getResolvedStressOnSource, 
     & insertionSortPlaneWithCrossing, annihilateDislocations,
     & annihilateDislocationsSub, updateDislPos, updateDislPosSub,
     & updateSources, updateDislRelPos, updateSource, createDipole,
     & updateDislPosSubNoAtoms, updateDislPosSub_ptr, assignDD

      procedure(Dummy), pointer :: updateDislPosSub_ptr => NULL()
      
C     HARD-CODED CONSTANTS
      real(dp), parameter :: OBSFUDGE = 0.9_dp ! see enforceObstacles
      
      contains
************************************************************************
      subroutine assignDD(simtype)

C     Subroutine: assignDD

C     Inputs: simtype --- simulation type: atomistic, fe, dd, cadd_nodisl, or cadd

C     Outputs: None

C     Purpose: Assign pointer for subroutines used in DD simulations
C     These may depend on simulation type
C     (e.g. updating dislocation position is different in dd and cadd simulations,
C     since in the former there is no atomistic domain to pass the dislocation into)
      
      implicit none
      
C     input variables
      character(len=*) :: simtype
      
      select case (trim(simtype))
          case ('cadd')
              updateDislPosSub_ptr => updateDislPosSub
          case ('dd')
              updateDislPosSub_ptr => updateDislPosSubNoAtoms
          case ('atomistic','cadd_nodisl','fe')
              return    
          case default
              write(*,*) 'Simulation type has not yet been defined'
              stop
      end select
      
      call assignDispFromPKOneMat(dislmisc%gradientcorrection)
      
      end subroutine assignDD
************************************************************************
      subroutine runDDStep(dt)

C     Subroutine: runDDStep

C     Inputs: dt --- time increment for DD update

C     Outputs: None

C     Purpose: Take a single step in main DD routine. Several substeps:
C     1) sort dislocations on planes (avoid unnecessary computation by only sorting
C     if resort is .true.)
C     2) compute displacements of dislocations using PK stress
C     3) update dislocation attributes, after enforcing obstacles, annihilating
C     dislocations, etc.
C     4) compute shear stresses on sources, update their timer, pop dislocation
C     dipole if abs(tau) > taucr for t > tnuc
C     5) zero out flags for obstacles, displacements for dislocations
      
      implicit none
      
C     input variables
      real(dp) :: dt
      
      call sortDislPlanes()
      call dispfromPK(dt)
      call updateDislocations()
      call updateSources(dt)
      call zeroObstacles()
      call zeroDislDisp()
      
      end subroutine runDDStep
************************************************************************
      subroutine updateDislocations()
      
C     Subroutine: updateDislocations

C     Inputs: None

C     Outputs: None

C     Purpose: Update dislocation structure after Peach-Koehler step, by enforcing constitutive
C     relations and keeping track of escaped dislocations. Several steps:
C     1) prevent dislocations from crossing obstacles, if obstacle is active
C     2) update relative position of dislocations along slip plane
C     3) using updated relative positions, check for disl. crossings
C        if dislocations of opposite sign cross, annihilate them
C     4) loop through slip plane again, checking for opposite sign dislocations
C        within Lannih of each other, and annihilate them if so
C     5) update remaining dislocation attributes (FE element, local position within element, x, y, etc.)

C     Notes: Have to be very careful with local vs. global copies of splane!
      
      implicit none
      
C     local variables
      integer :: i, j, k
      integer :: mnum
      real(dp) :: lannih
      
      do i = 1, size(disl)
          mnum = fematerials%list(i)
          lannih = materials(mnum)%lannih
          do j = 1, size(disl(i)%splanes)
              do k = 1, size(disl(i)%splanes(j)%splane)                  
                  call enforceObstacles(i,j,k,
     &                                  disl(i)%splanes(j)%splane(k))     
                  call updateDislRelpos(i,disl(i)%splanes(j)%splane(k))                  
                  call insertionSortPlaneWithCrossing(i,j,
     &                                     disl(i)%splanes(j)%splane(k))     
                  call annihilateDislocations(i,j,
     &                              disl(i)%splanes(j)%splane(k),lannih)
                  call updateDislPos(i,j,k)                  
              end do
          end do 
      end do
      
      end subroutine updateDislocations
************************************************************************
      subroutine enforceObstacles(mnumfe,isys,iplane,splane)
      
C     Subroutine: enforceObstacles

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system
C             iplane --- index of slip plane within slip system

C     In/out: splane --- sorted plane structure, disl(mnumfe)%splanes(isys)%splane(iplane)

C     Outputs: None

C     Purpose: Adjust displacements of dislocations on a particular slip plane
C              to avoid crossing active obstacles

C     Notes: We could avoid passing in splane, but then we'd have to write disl(mnumfe)%splanes(isys)%splane(iplane)
C     in front of everything (local copy is not possible as well, since we have to update global/module variable)
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: isys, iplane
      
C     in/out variables
      type(sortedplanedata) :: splane
      
C     local variables
      integer :: i
      integer :: dislnum
      real(dp) :: tau, taucr
      real(dp) :: disp, dispnew
      real(dp) :: relposold, relposnew, relposobs, relposdiff
      logical :: between
      integer :: obsnum
      logical :: active
      
      do i = 1, splane%nmax
          dislnum = splane%objnum(i)
          relposold = splane%relpos(i)
          disp = disl(mnumfe)%list(dislnum)%disp
          relposnew = relposold + disp
          call isObstacleBetween(relposold,relposnew,
     &                  obstacles(mnumfe)%splanes(isys)%splane(iplane),
     &                  between,obsnum,relposobs) ! check if obstacle is between new and old position
          if (between) then
              if (obstacles(mnumfe)%list(obsnum)%computed) then ! already computed
                  active = obstacles(mnumfe)%list(obsnum)%active
              else
                  tau = getResolvedStressOnObstacle(mnumfe,
     &                                 obstacles(mnumfe)%list(obsnum))
                  taucr = obstacles(mnumfe)%list(obsnum)%taucr
                  active = (abs(tau) < taucr) ! obstacle is active if tau < stress to break obstacle
                  obstacles(mnumfe)%list(obsnum)%computed = .true.
                  obstacles(mnumfe)%list(obsnum)%active = active
              end if
              if (active) then ! obstacle is active, adjust disp
                  relposdiff = relposobs - relposold
                  if (abs(relposdiff) < TOLCONST) then
                      dispnew = 0.0_dp ! if dislocation is very close to obstacle, don't move it
                  else
                      dispnew = OBSFUDGE*relposdiff ! otherwise, move it fractionally closer
                  end if
                  disl(mnumfe)%list(dislnum)%disp = dispnew
              end if
          end if
      end do
      
      end subroutine enforceObstacles
************************************************************************
      subroutine isObstacleBetween(pold,pnew,splane,
     &                             between,obsnum,pobs)

C     Subroutine: isObstacleBetween

C     Inputs: pold, pnew --- old and new relative positions along slip plane

C     In/out: splane --- sorted plane structure, disl(mnumfe)%splanes(isys)%splane(iplane)

C     Outputs: between --- logical indicating whether there's an obstacle between pold and pnew
C              obsnum --- if there is an obstacle, the obstacle's number
C              pobs --- if there is an obstacle, the obstacle's position

C     Purpose: Determine whether an obstacle lies between two points on a particular slip plane
C              Returns obstacle *closest* to pold
      
      implicit none
      
C     input variables
      real(dp) :: pold, pnew
      
C     in/out variables
      type(sortedplanedata) :: splane
      
C     output variables
      logical :: between
      integer :: obsnum
      real(dp) :: pobs
      
C     local variables
      integer :: n, iobs
      
      n = splane%ncount
      if (n > 0) then
          call findPointBetween(pold,pnew,splane%relpos(1:n),
     &                          n,between,iobs)
          pobs = splane%relpos(iobs)
          obsnum = splane%objnum(iobs)   
      else
          between = .false.
      end if
      
      end subroutine isObstacleBetween
************************************************************************
      subroutine updateDislRelpos(mnumfe,splane)

C     Subroutine: updateDislRelpos

C     Inputs: mnumfe --- fe material number

C     In/out: splane --- sorted plane structure containing dislocations

C     Purpose: Update relative positions of dislocations on slip plane
C     using the (adjusted) displacements

C     Notes/TODO: Is active check necessary, since splane has alreayd been sorted?

      implicit none
      
C     input variables
      integer :: mnumfe
      
C     in/out variables
      type(sortedplanedata) :: splane
      
C     local variables
      integer :: i
      integer :: dislnum
      real(dp) :: disp
      
      do i = 1, splane%nmax
          dislnum = splane%objnum(i)
          disp = disl(mnumfe)%list(dislnum)%disp
          if (disl(mnumfe)%list(dislnum)%active) then 
              splane%relpos(i) = splane%relpos(i) + disp
          end if    
      end do
      splane%resort = .true.
      
      end subroutine updateDislRelpos
************************************************************************
      subroutine insertionSortPlaneWithCrossing(mnumfe,isys,splane)

C     Subroutine: insertionSortPlaneWithCrossing

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system

C     In/out: splane --- sorted plane structure containing dislocations

C     Purpose: Sort dislocations on plane (using insertion sort),
C     keeping track of dislocation crossings and annihilating dislocations
C     of opposite sign that cross each other

C     Notes/TODO: Crossing of dislocations of same sign is ignored.
C     Is this correct?

C     Notes: Checking active is necessary since disl. can be deactivated via deletino
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: isys
      type(sortedplanedata) :: splane
      
C     local variables
      integer :: i, j
      real(dp) :: relpostemp
      integer :: dislnumtemp, dislnumother
      integer :: bsgntemp, bsgnother
      logical :: activetemp, activeother
      logical :: placedback
      
      splane%resort = .false. ! presumably, the list will be sorted at the end; however, if we annihilate, this is not so (taken care of in deleteDislocationSub2)
      do i = 2, splane%nmax
          j = i - 1
          relpostemp = splane%relpos(i)
          dislnumtemp = splane%objnum(i)
          activetemp = disl(mnumfe)%list(dislnumtemp)%active
          bsgntemp = disl(mnumfe)%list(dislnumtemp)%sgn
          placedback = .false.
          do while (j >= 1) 
              if (splane%relpos(j) <= relpostemp) then ! could combine with do while statement, but Fortran does not support short-circuiting
                  exit
              end if
C             need to swap; check for crossing of disl. of opposite signs
              dislnumother = splane%objnum(j)
              activeother = disl(mnumfe)%list(dislnumother)%active
              if (activetemp.and.activeother) then
                  bsgnother = disl(mnumfe)%list(dislnumother)%sgn
                  if (bsgntemp /= bsgnother) then ! opposite sign dislocations, need to annihilate
C                     place dislocation back before annihilating
                      splane%relpos(j+1) = relpostemp
                      splane%objnum(j+1) = dislnumtemp
                      placedback = .true.
                      call annihilateDislocationsSub(mnumfe,isys,
     &                                               splane,j,j+1)     
                      exit
                  end if
              end if
              splane%relpos(j+1) = splane%relpos(j)
              splane%objnum(j+1) = splane%objnum(j)
              j = j - 1
          end do
          if (.not.placedback) then ! if dislocation wasn't deleted, place it back
              splane%relpos(j+1) = relpostemp
              splane%objnum(j+1) = dislnumtemp
          end if
      end do
      
      call sortPlaneCheck(splane) ! resort plane, if resort = .true.
                                  ! (this is necessary before annihilateDislocations can be used)
      
      end subroutine insertionSortPlaneWithCrossing
************************************************************************
      subroutine annihilateDislocations(mnumfe,isys,splane,lannih)

C     Subroutine: annihilateDislocations

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system
C             lannih --- critical distance for annihilation of disl. of opposite signs

C     In/out: splane --- sorted plane structure containing dislocations

C     Purpose: Annihilate any remaining
C     opposite signed pairs within lannih of each other (plane must already be sorted)
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: isys
      real(dp) :: lannih

C     in/out variables
      type(sortedplanedata) :: splane
      
C     local variables
      integer :: counter
      logical :: delete
      integer :: dislnum1, dislnum2
      real(dp) :: relpos1, relpos2
      integer :: bsgn1, bsgn2
      integer :: nmax
      
      nmax = splane%nmax
      if (nmax > 1) then ! strictly unnecessary, but may save time
          counter = 1
          dislnum1 = splane%objnum(counter)
          relpos1 = splane%relpos(counter)
          bsgn1 = disl(mnumfe)%list(dislnum1)%sgn
          do while (counter < nmax)
              counter = counter + 1
              dislnum2 = splane%objnum(counter)
              relpos2 = splane%relpos(counter)
              bsgn2 = disl(mnumfe)%list(dislnum2)%sgn
              delete = ((bsgn1 /= bsgn2).and.           ! dislocations have to be of opposite sign
     &                  ((relpos2 - relpos1) < lannih)) ! and sufficiently close together
              if (delete) then
                  call annihilateDislocationsSub(mnumfe,isys,splane,
     &                                           counter-1,counter)
                  counter = counter + 1
                  if (counter < nmax) then
                      dislnum1 = splane%objnum(counter)
                      relpos1 = splane%relpos(counter)
                      bsgn1 = disl(mnumfe)%list(dislnum1)%sgn
                  end if
              else    
                  dislnum1 = dislnum2
                  relpos1 = relpos2
                  bsgn1 = bsgn2
              end if
          end do
      end if
      
      end subroutine annihilateDislocations
************************************************************************
      subroutine annihilateDislocationsSub(mnumfe,isys,splane,
     &                                     iobj1,iobj2)

C     Subroutine: annihilateDislocationsSub

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system
C             iobj1, iobj2 --- indices of dislocations along slip plane to be annihilated

C     In/out: splane --- sorted plane structure containing dislocations

C     Purpose: Annihilate pair of oppositely signed dislocations, paying
C     attention to whether the signs of the cuts are the same (if so, add escaped dislocations to account for displacement jump)

C     input variables
      integer :: mnumfe
      integer :: isys
      integer :: iobj1, iobj2
      
C     in/out variables
      type(sortedplanedata) :: splane
      
C     local variables
      integer :: dislnum1, dislnum2
      integer :: bcut1, bcut2
      integer :: bsgn1, bsgn2
      real(dp) :: p1(2), p2(2)

      dislnum1 = splane%objnum(iobj1)
      dislnum2 = splane%objnum(iobj2)      
      bcut1 = disl(mnumfe)%list(dislnum1)%cut
      bcut2 = disl(mnumfe)%list(dislnum2)%cut
      if (bcut1 /= bcut2) then ! create a slip step of magnitude b along entire plane
C         position is irrelevant, since cut extends along entire plane)
          p1 = disl(mnumfe)%list(dislnum1)%posn
          p2 = p1 + slipsys(mnumfe)%trig(:,isys)
          bsgn1 = disl(mnumfe)%list(dislnum1)%sgn
          bsgn2 = -bsgn1
          call addEscapedDislocation(p1,p2,isys,bsgn1,bcut1,mnumfe) ! positions are irrelevant (one branch cut will give nonzero answer, the other won't)
          call addEscapedDislocation(p1,p2,isys,bsgn2,bcut2,mnumfe) ! positions are irrelevant (one branch cut will give nonzero answer, the other won't)
      end if
      
C     might be dangerous to use deleteDislocation here because of side effects
      call deleteDislocationSub(mnumfe,dislnum1)
      call deleteDislocationSub2(mnumfe,splane,iobj1)
      call deleteDislocationSub(mnumfe,dislnum2)
      call deleteDislocationSub2(mnumfe,splane,iobj2)
      
      end subroutine annihilateDislocationsSub
************************************************************************
      subroutine updateDislPos(mnumfe,isys,iplane)

C     Subroutine: updateDislPos

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system
C             iplane --- index of slip plane within slip system

C     Purpose: Update disl. position (x, y) and associated attributes (localpos, element, etc.)
C     using (adjusted) disl. displacement along slip plane, for *all* dislocations

      implicit none

C     input variables
      integer :: mnumfe
      integer :: isys, iplane
      
C     local variables
      integer :: i
      integer :: dislnum

      do i = 1, disl(mnumfe)%splanes(isys)%splane(iplane)%nmax
          dislnum = disl(mnumfe)%splanes(isys)%splane(iplane)%objnum(i)
          if (disl(mnumfe)%list(dislnum)%active) then
              call updateDislPosSub(mnumfe,isys,iplane,i,dislnum)
          end if
      end do
      disl(mnumfe)%splanes(isys)%splane(iplane)%resort = .true.
      
      end subroutine updateDislPos
************************************************************************
      subroutine Dummy(mnumfe,isys,iplane,iobj,dislnum)
 
      implicit none
 
C     input variables 
      integer :: mnumfe
      integer :: isys, iplane, iobj, dislnum      
      
      end subroutine Dummy
************************************************************************
      subroutine updateDislPosSub(mnumfe,isys,iplane,iobj,dislnum)

C     Subroutine: updateDislPosSub

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system
C             iplane --- index of slip plane within slip system

C     Purpose: Update disl. position (x, y) and associated attributes (localpos, element, etc.)
C     using (adjusted) disl. displacement along slip plane, for *single* dislocation

C     Notes/TODO: multimaterial extension is not correct;
C     assumes slip systems have same numbers (isys) in different materials

      implicit none
     
C     input variables
      integer :: mnumfe
      integer :: isys, iplane, iobj, dislnum
     
C     local variables
      integer :: mnumfenew, element
      real(dp) :: r, s
      logical :: badflip
      real(dp) :: dislposold(2), dislposnew(2)
      integer :: bcut, bsgn
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      real(dp) :: disp(2)

      bcut = disl(mnumfe)%list(dislnum)%cut 
      dislposold = disl(mnumfe)%list(dislnum)%posn
      bsgn = disl(mnumfe)%list(dislnum)%sgn
      element = disl(mnumfe)%list(dislnum)%element
      
C     first, compute new position
      disp = invResolveDisp(mnumfe,isys,disl(mnumfe)%list(dislnum)%disp)
      dislposnew = dislposold + disp
      
C     then, figure out if dislocation is still in mesh
      mnumfenew = mnumfe
      call findInAllWithGuess(dislposnew(1),dislposnew(2),
     &                        mnumfenew,element,r,s,badflip)
      
      if (badflip) then ! not in mesh
          write(*,*) 'Badflip', badflip
C         Few possibilities:
C         1) Crossed back to atomistic region
C              a) Can be placed properly in atomistic region (past detection band)
C              b) Cannot be (reenters continuum or exits free surface)
C         2) Left mesh, leaving a slip step (associated with escaped dislocation)
C         Algorithm:
C         1) Check first possibility by seeing if it crossed atomistic-continuum interface
C         2) If so, passContinuumtoAtomistic takes care of logic in 1a, 1b
C         3) If not, then second possibility has occurred: so, delete the discrete dislocation and add an escaped one
 
          call findInterfaceIntersectionDeformed(interfaceedges%array,
     &                        dislposold,dislposnew,pint,isint,edgenum)
          if (isint) then ! in atomistic region
              call passContinuumToAtomistic(dislposold,pint,mnumfe,
     &                                      isys,iplane,iobj)
          else ! dislocation left mesh, add escaped dislocation
              call deleteDislocation(mnumfe,isys,iplane,iobj)
              call addEscapedDislocation(dislposold,dislposnew,
     &                                   isys,bsgn,bcut,mnumfe)
          end if
      else ! still in mesh
          if (mnumfenew /= mnumfe) then ! moved to different material
              call addDislocation(mnumfenew,element,dislposnew(1),
     &                            dislposnew(2),isys,bsgn,bcut) ! assumes slip systems have same numbers (isys) in different materials
              call deleteDislocation(mnumfe,isys,iplane,iobj)
          else ! in the same material; simply update position
              disl(mnumfe)%list(dislnum)%posn = dislposnew
              disl(mnumfe)%list(dislnum)%element = element
              disl(mnumfe)%list(dislnum)%localpos = [r,s]
          end if
      end if    
     
      end subroutine updateDislPosSub
************************************************************************
      subroutine updateDislPosSubNoAtoms(mnumfe,isys,iplane,
     &                                   iobj,dislnum)

C     Subroutine: updateDislPosSubNoAtoms.

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system
C             iplane --- index of slip plane within slip system

C     Purpose: Update disl. position (x, y) and associated attributes (localpos, element, etc.)
C     using (adjusted) disl. displacement along slip plane, for *single* dislocation.
C     Assumes no atoms, so no passing into atomistic region

C     Notes/TODO: multimaterial extension is not correct;
C     assumes slip systems have same numbers (isys) in different materials

      implicit none
     
C     input variables
      integer :: mnumfe
      integer :: isys, iplane, iobj, dislnum
     
C     local variables
      integer :: mnumfenew, element
      real(dp) :: r, s
      logical :: badflip
      real(dp) :: dislposold(2), dislposnew(2)
      integer :: bcut, bsgn
      real(dp) :: disp(2)

      bcut = disl(mnumfe)%list(dislnum)%cut 
      dislposold = disl(mnumfe)%list(dislnum)%posn
      bsgn = disl(mnumfe)%list(dislnum)%sgn
      element = disl(mnumfe)%list(dislnum)%element
      
C     first, compute new position
      disp = invResolveDisp(mnumfe,isys,disl(mnumfe)%list(dislnum)%disp)
      dislposnew = dislposold + disp
      
C     then, figure out if dislocation is still in mesh
      mnumfenew = mnumfe
      call findInAllWithGuess(dislposnew(1),dislposnew(2),
     &                        mnumfenew,element,r,s,badflip)
      
      if (badflip) then ! not in mesh
          call deleteDislocation(mnumfe,isys,iplane,iobj)
          call addEscapedDislocation(dislposold,dislposnew,
     &                               isys,bsgn,bcut,mnumfe)
          
      else ! still in mesh
          if (mnumfenew/=mnumfe) then ! moved to different material
              call addDislocation(mnumfenew,element,dislposnew(1),
     &                            dislposnew(2),isys,bsgn,bcut) ! assumes slip systems have same numbers (isys) in different materials
              call deleteDislocation(mnumfe,isys,iplane,iobj)
          else ! in the same material; simply update position
              disl(mnumfe)%list(dislnum)%posn = dislposnew
              disl(mnumfe)%list(dislnum)%element = element
              disl(mnumfe)%list(dislnum)%localpos = [r,s]
          end if
      end if    
     
      end subroutine updateDislPosSubNoAtoms
************************************************************************
      subroutine updateSources(dt)

C     Subroutine: updateSources

C     Inputs: dt --- time increment for DD update

C     Outputs: None

C     Purpose: Update timer for sources after computing resolved shear stress on each source
      
      implicit none
      
C     input variables
      real(dp) :: dt
      
C     local variables
      integer :: i, j
      real(dp) :: tau
      
      do i = 1, size(sources)
          do j = 1, size(sources(i)%list)
              tau = getResolvedStressOnSource(i,j)
              call updateSource(dt,tau,i,sources(i)%list(j))
          end do
      end do
      
      end subroutine updateSources
************************************************************************
      subroutine updateSource(dt,tau,mnumfe,source)

C     Subroutine: updateSources

C     Inputs: dt --- time increment for DD update
C             tau --- resolved shear stress on source
C             mnumfe --- fe material number

C     In/out: source --- dislocation source (type sourcet)

C     Outputs: None

C     Purpose: Update timer for individual source, popping dislocation dipole
C     if t > tnuc

      implicit none

C     input variables
      real(dp) :: dt
      real(dp) :: tau
      integer :: mnumfe
      
C     in/out variables
      type(sourcet) :: source

C     local variables
      real(dp) :: time
      real(dp) :: tauprev

      tauprev = source%tauprev
      time = source%time
      if (abs(tau) > source%taucr) then
          if (time < TOLCONST) then ! start timer
              time = time + dt
          else
              if (sameSign(tau,tauprev)) then ! same sign, so increment up
                  time = time + dt
                  if (time >= source%tnuc) then ! pop dislocation dipole
                      call createDipole(mnumfe,source,tau)
                      time = 0.0_dp ! reset timer
                  end if
              else ! restart timer
                  time = dt
              end if
          end if
      else
          if (time > TOLCONST) then
              time = max(time - dt,0.0_dp) ! decrement back to zero
              tau = tauprev ! keep record of previous tau > taucr stress, in case tau > taucr in future
          end if    
      end if
      source%tauprev = tau
      source%time = time
      
      end subroutine updateSource
************************************************************************
      subroutine createDipole(mnumfe,source,tau)

C     Subroutine: createDipole

C     Inputs: mnumfe --- fe material number
C             tau --- local shear stress on source

C     In/out: source --- dislocation source (type sourcet)

C     Outputs: None

C     Purpose: Create dislocation dipole corresponding to source, local shear stress
C     (dipole separation depends on material properties for material mnumfe)

      implicit none

C     input variables
      integer :: mnumfe
      type(sourcet) :: source
      real(dp) :: tau
      
C     local variables
      integer :: isys
      real(dp) :: posn(2), posnpos(2), posnneg(2), dpos(2)
      integer :: element

      posn = source%posn
      isys = source%slipsys
      element = source%element
      dpos = invResolveDisp(mnumfe,isys,0.5_dp*source%lnuc)
      if (tau < 0.0_dp) then ! flip dislocations if tau is negative
          dpos = -dpos
      end if
      posnpos = posn + dpos
      posnneg = posn - dpos  
      
      call addDislocation(mnumfe,element,
     &          posnpos(1),posnpos(2),isys,1,0)
      call addDislocation(mnumfe,element,
     &          posnneg(1),posnneg(2),isys,-1,0)
     
      end subroutine createDipole
************************************************************************
      function getResolvedStressOnSource(mnumfe,sourcenum) result(tau)

C     Function: getResolvedStressOnSource

C     Inputs: mnumfe --- fe material number
C             source --- number of dislocation source

C     Outputs: tau --- resolved shear stress on source

      implicit none

C     input variables
      integer :: mnumfe
      integer :: sourcenum
      
C     output variables
      real(dp) :: tau

C     local variables
      real(dp) :: r, s
      real(dp) :: stresstilde(3), stresshat(3), stress(3)
      type(sourcet) :: source
      
      source = sources(mnumfe)%list(sourcenum)
      r = source%localpos(1)
      s = source%localpos(2)
      stresstilde = getTildeStressOnSourceAll(sourcenum,mnumfe) ! exclude self-stress from latent dislocations associated with source itself
      stresshat = getFEStressAtPoint(mnumfe,source%element,r,s)
      stress = stresshat + stresstilde
      tau = resolveStress(mnumfe,source%slipsys,stress)
      
      end function getResolvedStressOnSource
************************************************************************
      function getResolvedStressOnObstacle(mnumfe,obstacle) result(tau)

C     Function: getResolvedStressOnSource

C     Inputs: mnumfe --- fe material number
C             obstacle --- dislocation obstacle (type obstaclet)

C     Outputs: tau --- resolved shear stress on obstacle

      implicit none      

C     input variables
      integer :: mnumfe
      type(obstaclet) :: obstacle
      
C     output variables
      real(dp) :: tau

C     local variables
      real(dp) :: r, s
      real(dp) :: stresstilde(3), stresshat(3), stress(3)
      
      r = obstacle%localpos(1)
      s = obstacle%localpos(2)
      stresstilde = getTildeStressAtPointAll(obstacle%posn,mnumfe)
      stresshat = getFEStressAtPoint(mnumfe,obstacle%element,r,s)
      stress = stresshat + stresstilde
      tau = resolveStress(mnumfe,obstacle%slipsys,stress)
      
      end function getResolvedStressOnObstacle
************************************************************************
      end module mod_dd_main