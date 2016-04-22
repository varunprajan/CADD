      module mod_moving_mesh_crack_cadd
      
C     Purpose: Implements moving mesh technique for CADD crack problem, similar
C     to that outlined in Chakravarthy and Curtin, 2011, MSMSE for the pure DD crack problem.
C     In the moving mesh method, the problem is shifted so that a crack at (xcrack, ycrackpos)
C     is moved to (xsmall, ycrackpos), where xsmall is close to zero
C     and ycrackpos probably is as well. To do this, displacements are found at nodes
C     in the original mesh by finding the corresponding displacement
C     at a point shifted by (xcrack, ycrackpos) using interpolation.
C     Also, dislocations are shifted by the same amount. In the end,
C     the result is that all of the field variables have been shifted by xcrack,
C     the crack is back at (xsmall, ycrackpos), where xsmall is very close to 0,
C     and the mesh is unchanged (stiffness matrix is the same, etc.).

C     Long note: There are many limitations/potential problems associated with the moving mesh technique for crack simulations.
C     1) No sources/obstacles are permitted to be directly ahead of the atomistic region,
C     since these objects are artificial and do not have any direct counterpart
C     in the atomistic region (unlike pure dislocations).
C     2) Dislocations that have escaped the body "ahead" of the crack tip
C     (in region REGIONMAINBODY, see mod_disl_escaped) can end up moving
C     "behind" the crack tip, leaving spurious steps on the crack plane.
C     3) Dislocations, sources, and obstacles far behind the crack *cannot* be moved out of the simulation cell;
C     thus, the zone for these objects should not extend too far behind the crack
C     4) Dislocations ahead of the crack, in the continuum region,
C     must be moved into the atomistic region, possible outside of the detection band.
C     The dislocations outside of the detection band are not permitted to be there.
C     So, the detection band is temporarily expanded *outward* to cover the entire atomistic region,
C     and dislocations are detected and passed out. Whether this is a good solution is unclear.
C     5) We also want to prevent atomistic dislocations from being
C     moved into the continuum (since the continuum obviously can't represent
C     the dislocation properly). Thus, before the moving mesh technique is applied,
C     the detection band is expanded *inward*, and dislocations are detected and passed out.
C     Again, whether this is a good solution is unclear.
C     6) The interpolation of displacements may fail in regions where the
C     deformation is non-affine. This may happen for large xshift = xcrack - xsmall
C     and when dislocation cuts are present in the left part of the domain
C     (which might happen if dislocations can be emitted backwards).

C     Issue #3 can probably be resolved with a better algorithm,
C     where dislocations are passed out of the simulation cell to become
C     escaped dislocations, and sources/obstacles are deactivated when
C     they leave the simulation cell (although this might cause unpinning
C     of dislocation pileups, etc.)

C     Notes/TODO: This module needs to be thoroughly tested using integration tests.
C     I'm particularly worried about the case where a discrete dislocation
C     leaves the simulation cell during the shifting process.

C     Other notes/TODO: Do dislocation sources and obstacles need to be
C     "replenished" as they are lost from the region ahead of the crack?

      use mod_utils, only: prettyPrintMat
      use mod_math, only: TOLCONST      
      use mod_types, only: dp
      use mod_nodes, only: nodes, getXYAtomBoundsDef
      use mod_misc, only: misc
      use mod_fe_el_2d, only: felib, getElTypeNum
      use mod_fe_elements, only: feelements
      use mod_fe_main_2d_assign, only: getTotalDispAtPoint_ptr,
     &                                 solveAll_ptr
      use mod_math, only: nearestMultiple, checkSameSide
      use mod_disl_try, only: disl, sources, obstacles,
     &                        addDislocation, deleteDislocation
      use mod_disl_escaped, only: escapeddisl
      use mod_disl_ghost, only: ghostdisl
      use mod_slip_sys, only: slipsys
      use mod_neighbors, only: genAtomBinsUndeformed, neighbors,
     &                isInAtomBox, getAtomBin, updateNeighborsNoCheck
      use mod_delaunay, only: delaunaydata, genDelaunay,
     & setDelaunayPosUndef
      use mod_mesh_find, only: findInAllWithGuessDef, getLocalCoords,
     & findInAllInitiallyUndef, findInAllWithGuessUndef
      use mod_disl_detect_pass, only: detectAndPassDislocations,
     & detection, processDetectionBand
      use mod_pad_atoms, only: updatePad
      use mod_kdispfield, only: applyKDispIsoSet
     
      implicit none
      
      type movingmeshdata
C     read-in
      real(dp) :: deltaxcrack
      real(dp) :: maxshift
C     processed
      type(delaunaydata) :: delaunay
      real(dp) :: xshift
      integer :: eltypenum
      end type
      
C     module variables (global)
      type(movingmeshdata) :: movingmesh    
      
      private
      public :: shiftPosnDisp, initCADDMovingMeshData,
     &  interpFromFEPoint, interpFromAtomPoint, readCADDMovingMeshData,
     &  writeCADDMovingMeshData, passDislocationsBefore, movingmesh,
     &  passDislocationsAfter, passDislocationsSub, checkTriangle,
     &  locateClosestAtomUndef, shiftDDPos, shiftSlipPlanes,
     &  updateDislPosMovingMesh, updateDislPosMovingMeshSub,
     &  updateSourcePosMovingMesh, updateObstaclePosMovingMesh,
     &  processCADDMovingMeshData, getNewDetectionBandAfter,
     &  getNewDetectionBandBefore, isInAtomBox, atomicDispInterpolation,
     &  assignCrackXShift, shiftAllPosn, updateAllShiftedDDObjects

      contains
************************************************************************
      subroutine initCADDMovingMeshData(movingmeshfile)
 
C     Inputs: movingmeshfile --- filename where moving mesh information is stored
C     (should be something like '[filepref]_caddmovingmesh')
 
C     Outputs: None
 
C     Purpose: Read, initialize data in "movingmesh" structure, which holds
C     information about the moving mesh method for CADD problems
      
      implicit none
      
C     input variables
      character(len=*) :: movingmeshfile
      
      call readCADDMovingMeshData(movingmeshfile)
      call processCADDMovingMeshData()      
      
      end subroutine initCADDMovingMeshData
************************************************************************
      subroutine readCADDMovingMeshData(movingmeshfile)
 
C     Inputs: movingmeshfile --- filename where moving mesh information is stored
C     (should be something like '[filepref]_caddmovingmesh')
 
C     Outputs: None
 
C     Purpose: Read moving mesh data from file
      
      implicit none
      
C     input variables
      character(len=*) :: movingmeshfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=movingmeshfile)
      read(iunit,*) movingmesh%deltaxcrack
      read(iunit,*) movingmesh%maxshift
      
      close(iunit)      
      
      end subroutine readCADDMovingMeshData
************************************************************************
      subroutine processCADDMovingMeshData()
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Set up node numbers for delaunay triangulation (used for interpolation),
C     etc.
      
      movingmesh%delaunay%nodenums = nodes%atomlist ! moving mesh triangulation uses all atoms, including pad
      movingmesh%eltypenum = getElTypeNum('CPE3') ! triangle
      
      end subroutine processCADDMovingMeshData
************************************************************************
      subroutine writeCADDMovingMeshData(movingmeshfile)
 
C     Inputs: movingmeshfile --- filename where moving mesh information is stored
C     (should be something like '[filepref]_caddmovingmesh')
 
C     Outputs: None
 
C     Purpose: Write moving mesh data to file (essentially
C     inverse of readCADDMovingMeshData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: movingmeshfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=movingmeshfile)
      write(iunit,*) movingmesh%deltaxcrack
      write(iunit,*) movingmesh%maxshift
      
      close(iunit)      
      
      end subroutine writeCADDMovingMeshData
************************************************************************
      subroutine assignCrackXShift(crackpos)
      
C     Inputs: crackpos --- current x, y coordinates of crack (real vector, length 2)
C     
C     Outputs: None
C     
C     Purpose: compute x-shift of the crack (= xsmall - xcrackpos).
C     Shift must be multiple of deltaxcrack, and no larger than maxshift
            
      implicit none
            
C     input variables
      real(dp) :: crackpos(2)
      
      movingmesh%xshift = -nearestMultiple(crackpos(1),
     &                                     movingmesh%deltaxcrack) ! shift must be integer multiple of deltaxcrack
                                                                    ! shift is backwards (hence negative sign)
                                                                   
      if (abs(movingmesh%xshift) > movingmesh%maxshift) then
          write(*,*) 'Crack shift is too large'
          stop
      end if
      
      end subroutine assignCrackXShift
************************************************************************
      subroutine shiftAllPosn()
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Assign new displacements to the nodes so that, in essence,
C     everything is shifted back by xshift, but we retain the old mesh (i.e., the stiffness matrix is the same).
C     Also, move all DD objects (dislocations of all types, sources, obstacles, and slip planes) back by xshift

      implicit none    

      call shiftPosnDisp()      
      call shiftDDPos()
      call shiftSlipPlanes()      
      
      end subroutine shiftAllPosn
************************************************************************
      subroutine updateAllShiftedDDObjects()
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Determine the new elements/local positions the dislocations, sources, obstacles
C     have moved to.

      implicit none    

      call updateDislPosMovingMesh()
      call updateSourcePosMovingMesh()
      call updateObstaclePosMovingMesh()      
      
      end subroutine updateAllShiftedDDObjects
************************************************************************
      subroutine passDislocationsBefore()
 
C     Inputs: None

C     Outputs: None
 
C     Purpose: Dislocations in the atomistic region must be passed out
C     before the moving mesh method is applied if they would otherwise
C     end up in the continuum region. The reason is that the continuum
C     cannot represent an atomistic dislocation properly with displacements alone
C     (we need an actual discrete dislocation object).

C     To do this, we expand the dislocation band inward, detect the dislocations
C     and pass them out. A problem may arise if the shift distance is too large, 
C     where we pass out dislocations very deep in the interior of the atomistic region.
C     So this possibility should be checked carefully
      
      implicit none
      
C     local variables
      real(dp) :: xmin, xmax, ymin, ymax      
      real(dp), allocatable :: paramsnew(:)
      real(dp) :: maxdisttointerface
      
C     get atom box
      call getXYAtomBoundsDef(xmin,xmax,ymin,ymax)
      maxdisttointerface = sqrt((xmax - xmin)**2 + (ymax - ymin)**2)
      
C     get new detection band
      paramsnew = getNewDetectionBandBefore(xmin,xmax,ymin,ymax)
      
C     pass
      call passDislocationsSub(maxdisttointerface,paramsnew)
      
      end subroutine passDislocationsBefore
************************************************************************
      function getNewDetectionBandBefore(xmin,xmax,ymin,ymax)
     &                                                result(paramsnew)
 
C     Inputs: xmin, xmax, ymin, ymax --- min, max x- and y-coordinates
C       of atomistic region in *deformed* configuration

C     Outputs: paramsnew --- new parameters for detection band (annulus),

C     Purpose: Construct a new detection band that captures all dislocations
C     that might otherwise be (spuriously) moved into the continuum. The strategy
C     is to shift the left inner boundary of the box to encompass more space.
C     This also causes the box's center to shift as well
     
      implicit none

C     input variables
      real(dp) :: xmin, xmax, ymin, ymax
      
C     output variables
      real(dp), allocatable :: paramsnew(:)

C     local variables
      real(dp) :: xc, xcnew
      real(dp) :: halfLxinner, halfLxinnernew
      real(dp) :: halfLxouter, halfLxouternew
      real(dp) :: xmininner1, xmininner2, xmininnernew
      real(dp) :: xminouter, xmaxinner
     
C     construct new detection band, with center and thickness shifted in x
C     (this is a bit complicated, and requires a picture to understand)
      paramsnew = detection%params
      if (detection%bandtype == 'rect_annulus') then
          xc = paramsnew(1)
          halfLxinner = paramsnew(3)
          halfLxouter = paramsnew(4)
          xmininner1 = xc - halfLxinner
          xmininner2 = xmin - movingmesh%xshift
          xmininnernew = max(xmininner1,xmininner2)
          xmaxinner = xc + halfLxinner
          xcnew = 0.5_dp*(xmininnernew + xmaxinner)
          halfLxinnernew = xcnew - xmininnernew
          xminouter = xc - halfLxouter                
          halfLxouternew = xcnew - xminouter
          paramsnew(1) = xcnew
          paramsnew(3) = halfLxinnernew
          paramsnew(4) = halfLxouternew
      else
          write(*,*) 'Method is undefined for non-rectangular bands'
          stop
      end if
      
      end function getNewDetectionBandBefore
************************************************************************
      subroutine passDislocationsAfter()
 
C     Inputs: None

C     Outputs: None
 
C     Purpose: After the moving mesh method is applied, dislocations
C     may reside outside of the detection band in the atomistic region.
C     (Due to continuum dislocations lying ahead of the atomistic region)
C     These cannot be detected, so they must be passed out before the 
C     normal CADD loop can resume.

C     To do this, we expand the dislocation band outward to cover the entire
C     atomistic region, and pass out any dislocations we find.

      implicit none

C     local variables
      real(dp) :: xmin, xmax, ymin, ymax      
      real(dp), allocatable :: paramsnew(:)
      real(dp) :: maxdisttointerface
      
C     get atom box
      call getXYAtomBoundsDef(xmin,xmax,ymin,ymax)
      maxdisttointerface = sqrt((xmax - xmin)**2 + (ymax - ymin)**2)

C     get new detection band
      paramsnew = getNewDetectionBandAfter(xmin,xmax,ymin,ymax) 
      
C     pass
      call passDislocationsSub(maxdisttointerface,paramsnew)
      
      end subroutine passDislocationsAfter
************************************************************************
      function getNewDetectionBandAfter(xmin,xmax,ymin,ymax)
     &                                                result(paramsnew)
 
C     Inputs: xmin, xmax, ymin, ymax --- min, max x- and y-coordinates
C       of atomistic region in *deformed* configuration

C     Outputs: paramsnew --- new parameters for detection band (annulus),

C     Purpose: Construct a new detection band that captures all dislocations
C     that might have been moved outside of the detection band as a result
C     of the moving mesh method. Basically, we simple expand the outer band
C     of the detection band to encompass the entire atomistic region.

      implicit none

C     input variables
      real(dp) :: xmin, xmax, ymin, ymax
      
C     output variables
      real(dp), allocatable :: paramsnew(:)
      real(dp) :: xc
      real(dp) :: xdist1, xdist2
      real(dp) :: halfLxouternew

C     construct new detection band, expanded outward to reach interface
      paramsnew = detection%params
      if (detection%bandtype == 'rect_annulus') then
          xc = paramsnew(1)
          xdist1 = xc - xmin
          xdist2 = xmax - xc
          halfLxouternew = max(xdist1,xdist2)
          paramsnew(4) = halfLxouternew
      else
          write(*,*) 'Method is undefined for non-rectangular bands'
          stop
      end if     
     
      end function getNewDetectionBandAfter
************************************************************************
      subroutine passDislocationsSub(maxdisttointerface,paramsnew)
 
C     Inputs: None

C     Outputs: None
 
C     Purpose: After the moving mesh method is applied, dislocations
C     may reside outside of the detection band in the atomistic region.
C     (Due to continuum dislocations lying ahead of the atomistic region)
C     These cannot be detected, so they must be passed out before the 
C     normal CADD loop can resume.

C     To do this, we expand the dislocation band outward to cover the entire
C     atomistic region, and pass out any dislocations we find.
      
      implicit none
      
C     input variables
      real(dp) :: maxdisttointerface
      real(dp) :: paramsnew(:)
      
C     local variables
      real(dp) :: tempmaxdisttointerface
      real(dp), allocatable :: tempparams(:), tempparamspadded(:)
      logical :: detected
      
C     regenerate neighbor lists
      call updateNeighborsNoCheck()
      
C     save old
      tempmaxdisttointerface = detection%maxdisttointerface
      tempparams = detection%params
      tempparamspadded = detection%paramspadded

C     assign new
      detection%maxdisttointerface = maxdisttointerface
      detection%params = paramsnew

C     regen detection band, pass dislocations
      call processDetectionBand()
      call detectAndPassDislocations(detected)
      
C     restore old parameters
      detection%maxdisttointerface = tempmaxdisttointerface
      detection%params = tempparams
      detection%paramspadded = tempparamspadded
      
      end subroutine passDislocationsSub
************************************************************************
      subroutine shiftPosnDisp()
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Re-interpolate the displacements so that, effectively,
C     the mesh, atoms, etc. are shifted back by xshift.
C     Operationally, we assign displacements to a node at point x
C     in the new mesh (with the crack back at the origin, or thereabouts)
C     by using the displacements of a node at point x + xshift
C     in the old mesh (with the crack at crackpos). x and x + xshift
C     correspond to the *undeformed* configuration
           
      implicit none
      
C     local variables
      real(dp), allocatable :: newposn(:,:)
      logical, allocatable :: alreadychecked(:)
      integer :: nodesshape(2)
      integer :: i, j, k
      integer :: node, nodetype
      integer :: m, n
      real(dp) :: posnundef(2), posnshift(2), shiftvec(2)
      real(dp) :: newdisp(2)
      
C     new posn
      nodesshape = shape(nodes%posn)
      m = nodesshape(1)
      n = nodesshape(2)
      allocate(newposn(m,n))
      newposn(3,:) = nodes%posn(3,:) ! z-position remains same
      newposn(6:7,:) = 0.0_dp ! velocities are zeroed
      
C     already checked nodes
      allocate(alreadychecked(nodes%nnodes))
      alreadychecked = .false.
      
C     shift vector (find corresponding point by shifting forwards)
      shiftvec = [-movingmesh%xshift,0.0_dp]
      
C     generate undeformed atomic bins, triangulation
      call genAtomBinsUndeformed()
      call setDelaunayPosUndef(movingmesh%delaunay)
      call genDelaunay(movingmesh%delaunay)
      
C     cycle over elements to generate good guesses
      do i = 1, size(feelements)
          do j = 1, size(feelements(i)%connect,2)
          do k = 1, size(feelements(i)%connect,1)
              node = feelements(i)%connect(k,j)
              nodetype = nodes%types(2,node)
              if (.not.alreadychecked(node).and.(nodetype==0)) then ! continuum node
                  posnundef = nodes%posn(1:2,node)-nodes%posn(4:5,node)
                  posnshift = posnundef + shiftvec
                  newdisp = interpFromFEPoint(posnshift,i,j)
                  newposn(1:2,node) = posnundef + newdisp
                  newposn(4:5,node) = newdisp
                  alreadychecked(node) = .true.
              end if
          end do
          end do
      end do
      
C     deal with rest of nodes (just atoms)
      do i = 1, nodes%nnodes
      if (.not.alreadychecked(i)) then
          posnundef = nodes%posn(1:2,i) - nodes%posn(4:5,i)
          posnshift = posnundef + shiftvec
          newdisp = interpFromAtomPoint(posnshift)
          newposn(1:2,i) = posnundef + newdisp
          newposn(4:5,i) = newdisp
      end if
      end do
      
      nodes%posn = newposn
      
      end subroutine shiftPosnDisp
************************************************************************
      function interpFromFEPoint(posnundef,mnumfeguess,elementguess)
     &                                                      result(disp)
 
C     Inputs: posnundef --- position of point in undeformed coordinates
C             mnumfe --- guess for fe material number that point is in
C             element --- guess for element that point is in
 
C     Outputs: disp --- displacement of point
 
C     Purpose: Find interpolated displacement of point, assuming that it lies within
C     the continuum domain (if that search fails, we search the atomistic domain instead)

      implicit none

C     input variables
      real(dp) :: posnundef(2)
      integer :: mnumfeguess
      integer :: elementguess
      
C     output variables
      real(dp) :: disp(2)
      
C     local variables
      integer :: eltypenum
      real(dp) :: r, s
      logical :: badflip, failure
      integer :: mnumfe, element

      mnumfe = mnumfeguess
      element = elementguess
      call findInAllWithGuessUndef(posnundef(1),posnundef(2), ! find in mesh (undeformed)
     &                             mnumfe,element,r,s,badflip)
      if (badflip) then ! if not successful, interpolate atomistic displacements
          call atomicDispInterpolation(posnundef,failure,disp)
          if (failure) then ! if still not successful, our point must be at the boundary,
                            ! so just use fake displacement (we'll correct boundary displacements later)
              disp = [0.0_dp,0.0_dp]
          end if
      else 
          eltypenum = feelements(mnumfe)%eltypenum
          disp = getTotalDispAtPoint_ptr(posnundef,
     &                                   mnumfe,eltypenum,element,r,s) ! compute displacement using fe + dd
      end if
              
      end function interpFromFEPoint
************************************************************************
      function interpFromAtomPoint(posnundef) result(disp)
 
C     Inputs: posnundef --- position of point in undeformed coordinates
 
C     Outputs: disp --- displacement of point
 
C     Purpose: Find interpolated displacement of point, assuming that it lies within
C     the atomistic domain (if that search fails, we search the continuum domain instead)
     
      implicit none

C     input variables
      real(dp) :: posnundef(2)
      
C     output variables
      real(dp) :: disp(2)

C     local variables
      integer :: mnumfe, element
      real(dp) :: r, s
      logical :: badflip, failure
    
      call atomicDispInterpolation(posnundef,failure,disp)
      if (failure) then ! point must be in mesh
          call findInAllInitiallyUndef(posnundef(1),posnundef(2), ! find in mesh (undeformed)
     &                                 mnumfe,element,r,s,badflip)
          if (badflip) then ! something went wrong
              write(*,*) 'Unable to find point'
              write(*,*) 'x', posnundef(1), 'y', posnundef(2)
              stop
          end if
      end if
      
      end function interpFromAtomPoint
************************************************************************
      subroutine atomicDispInterpolation(posnundef,failure,disp)
 
C     Inputs: posnundef --- position of point in undeformed coordinates

C     Outputs: failure --- flag indicating that search failed
C              disp --- displacement of point
 
C     Purpose: Find interpolated displacement of point, assuming that it lies within
C     atomistic region (if search fails, this is indicated in "failure").

C     Note: The position *must* lie very close to a vertex, since deformations
C     are not affine far within the atomistic region. Positions will lie
C     close to a vertex only if xshift is a multiple of the atomic spacing (deltaxcrack)
C     Some positions will not lie close to a vertex if they were originally
C     continuum positions (i.e., posnundef = posn_{continuum} + shiftvec),
C     but these are probably close to the interface where deformations should be affine.
      
      implicit none
      
C     input variables
      real(dp) :: posnundef(2)
      
C     output variables
      logical :: failure
      real(dp) :: disp(2)
      
C     local variables
      integer :: closestatom, closestnode
      real(dp) :: closestatompos(2), delpos(2)
      integer :: i, j
      
C     check if in atom box
      failure = (.not.isInAtomBox(posnundef(1),posnundef(2)))
      if (failure) then
          return
      end if
      
C     find closest atom in atom box
      closestatom = locateClosestAtomUndef(posnundef)
      
C     if atom is close enough, simply use its displacement
      closestnode = movingmesh%delaunay%nodenums(closestatom)
      closestatompos = movingmesh%delaunay%xy(:,closestatom)
      delpos = closestatompos - posnundef
      if (sqrt(sum(delpos**2)) < TOLCONST) then
          disp = nodes%posn(4:5,closestnode)
          failure = .false.
          return
      end if    
      
C     otherwise, find triangles that that atom belongs to, check if point lies in that triangle
      do i = 1, movingmesh%delaunay%numtri
          do j = 1, 3
              if (movingmesh%delaunay%trivert(j,i) == closestatom) then ! node belongs to a triangle, so check
                  call checkTriangle(posnundef,i,failure,disp)
                  if (.not.failure) then
                      return
                  end if        
              end if
          end do
      end do
     
      end subroutine atomicDispInterpolation
************************************************************************
      subroutine checkTriangle(posnundef,trinum,failure,disp)
 
C     Inputs: posnundef --- position of point in undeformed coordinates
C             trinum --- number of triangle in triangulation

C     Outputs: failure --- flag indicating that point was not in the triangle
C              disp --- displacement of point, linearly interpolating using vertices of triangle
C                       (only valid if failure == .false.)
 
C     Purpose: Check whether point lies within a triangle in the delaunay
C     triangulation. If it does, compute its displacement assuming deformation
C     is affine within the triangle.
      
      implicit none
      
C     input variables
      real(dp) :: posnundef(2)
      integer :: trinum
      
C     output variables
      logical :: failure
      real(dp) :: disp(2)
      
C     local variables
      real(dp) :: triposnundef(2,3)
      integer :: i
      integer :: atom, node
      integer :: idx1, idx2, idx3
      logical :: check
      real(dp), allocatable :: N(:)
      real(dp) :: udisp(3), vdisp(3)
      real(dp) :: r, s
      
C     triangle vertex positions, displacements
      do i = 1, 3
          atom = movingmesh%delaunay%trivert(i,trinum)
          node = movingmesh%delaunay%nodenums(atom)
          triposnundef(:,i) = movingmesh%delaunay%xy(:,atom)
          udisp(i) = nodes%posn(4,node)
          vdisp(i) = nodes%posn(5,node)
      end do
      
C     check all three sides
      do i = 1, 3
          idx1 = i
          idx2 = mod(idx1,3) + 1
          idx3 = mod(idx2,3) + 1
          check = checkSameSide(posnundef,triposnundef(:,idx3),
     &                       triposnundef(:,idx1),triposnundef(:,idx2))
          failure = (.not.check)
          if (failure) then ! we're not in triangle, so quit
              return
          end if    
      end do
      
C     if we've gotten this far, compute local position in the triangle
      call getLocalCoords(transpose(triposnundef),movingmesh%eltypenum,
     &                    posnundef(1),posnundef(2),r,s)
                          
C     displacement is given by shape function interpolation
      N = felib(movingmesh%eltypenum)%getN_2d_ptr(r,s)
      disp(1) = dot_product(N,udisp)
      disp(2) = dot_product(N,vdisp)               
      
      end subroutine checkTriangle
************************************************************************
      function locateClosestAtomUndef(posnundef) result(closestatom)
 
C     Inputs: posnundef --- position of point in undeformed coordinates

C     Outputs: closestatom --- index of atom (in nodes%atomlist) of closest atom to the point
 
C     Purpose: Finds closest atom to a point using atom bins
C              (already generated using getAtomBinsUndeformed)

      implicit none
      
C     input variables
      real(dp) :: posnundef(2)
      
C     output variables
      integer :: closestatom
      
C     local variables
      integer :: j, kx, ky
      integer :: nbinsx, nbinsy
      integer :: binx, biny, binxneigh, binyneigh
      integer :: atom, node
      real(dp) :: nodeposundef(2)
      real(dp) :: distmin, dist
      
C     get the atomic bin for our position
      call getAtomBin(posnundef(1),posnundef(2),binx,biny)
      
C     set distance to be very large value
      closestatom = 0
      distmin = huge(0.0_dp)
      
C     loop over atoms in neighboring/current bin
      nbinsx = size(neighbors%bincount,1)
      nbinsy = size(neighbors%bincount,2)
      do ky = -1, 1
      binyneigh = biny + ky
      do kx = -1, 1
      binxneigh = binx + kx
C     not outside box
      if ((binxneigh>=1).and.(binxneigh<=nbinsx).and.
     &    (binyneigh>=1).and.(binyneigh<=nbinsy)) then
C         loop through possible neighbors
          do j = 1, neighbors%bincount(binxneigh,binyneigh)
              atom = neighbors%binarray(j,binxneigh,binyneigh)
              node = nodes%atomlist(atom)
              nodeposundef = nodes%posn(1:2,node) - nodes%posn(4:5,node)
              dist = sum((posnundef - nodeposundef)**2)
              if (dist < distmin) then
                  closestatom = atom
                  distmin = dist
              end if    
          end do
      end if        
      end do
      end do
      
      end function
************************************************************************
      subroutine shiftDDPos()
 
C     Inputs: None

C     Outputs: None
 
C     Purpose: Shifts all DD objects by xshift: normal, ghost, and escaped
C     dislocations, sources, and obstacles

C     TODO: Does the shifting of ghost dislocations cause problems?

      implicit none
      
C     local variables
      integer :: i, j
      
C     normal dislocations
      do i = 1, size(disl)
          do j = 1, disl(i)%ndisl
          if (disl(i)%list(j)%active) then
              disl(i)%list(j)%posn(1) = disl(i)%list(j)%posn(1) +
     &                                  movingmesh%xshift
          end if
          end do
      end do
      
C     ghost dislocations
      do i = 1, size(ghostdisl)
          do j = 1, ghostdisl(i)%nghostdisl
              ghostdisl(i)%list(j)%posn(1) =
     &                  ghostdisl(i)%list(j)%posn(1) + movingmesh%xshift
          end do
      end do
      
C     escaped dislocations
      do i = 1, size(escapeddisl)
          do j = 1, escapeddisl(i)%nescapeddisl
              escapeddisl(i)%list(j)%posn(1) =
     &                escapeddisl(i)%list(j)%posn(1) + movingmesh%xshift
          end do
      end do
      
C     sources
      do i = 1, size(sources)
          do j = 1, size(sources(i)%list)
              sources(i)%list(j)%posn(1) =
     &                    sources(i)%list(j)%posn(1) + movingmesh%xshift
          end do
      end do    
      
C     obstacles
      do i = 1, size(obstacles)
          do j = 1, size(obstacles(i)%list)
              obstacles(i)%list(j)%posn(1) =
     &                  obstacles(i)%list(j)%posn(1) + movingmesh%xshift
          end do
      end do 
      
      end subroutine shiftDDPos
************************************************************************
      subroutine shiftSlipPlanes()
 
C     Inputs: None

C     Outputs: None
 
C     Purpose: Shifts all slip planes by shift
      
C     Long note: There are two ways to do this, each with their own advantages
C     and disadvantages.

C     1) The simpler/more elegant method is to simply shift 
C     the origin of the slip system by (xshift, 0).
C     The possible problem here is that, after successive shifts,
C     the origin may no longer lie outside of the body (see getSlipPlane).
C     In other words, slip planes may no longer cover the entire body if the origin is moved.
C     This method has the advantage that the slip plane spacing
C     can be arbitrary, and not necessarily a multiple of misc%deltaxshift

C     2) The other method is to move all of the objects (dislocations,
C     sources, obstacles) from plane isys to plane isys + disys, where disys*spacing
C     = xshift (with some trig factors). The relative positions may also need to be shifted.
C     The possible problem here is that isys + disys > nslipplanes
C     for the last few slip planes. So, either this array has to be reallocated,
C     or these dislocations must be eliminated at an earlier
C     stage. Or, dislocation activity should not be allowed on the last few planes...
C     This method requires that the slip plane spacing be a multiple of misc%deltaxshift

C     I have implemented method 1 below.

      implicit none
      
C     local variables
      integer :: i, j
      integer:: nslipsys

      do i = 1, size(slipsys)
          nslipsys = size(slipsys(i)%theta)
          do j = 1, nslipsys
              slipsys(i)%origin(1,j) = slipsys(i)%origin(1,j) +
     &                                   movingmesh%xshift
          end do
      end do      
      
      end subroutine shiftSlipPlanes
************************************************************************
      subroutine updateDislPosMovingMesh()

C     Inputs: None

C     Outputs: None

C     Purpose: Find new element and local position for all dislocations

      implicit none
      
C     local variables
      integer :: i, j, k, l
      integer :: dislnum
      
      do i = 1, size(disl)
      do j = 1, size(disl(i)%splanes)
      do k = 1, size(disl(i)%splanes(j)%splane)
      do l = 1, disl(i)%splanes(j)%splane(k)%nmax
          dislnum = disl(i)%splanes(j)%splane(k)%objnum(l)
          if (disl(i)%list(dislnum)%active) then
              call updateDislPosMovingMeshSub(i,j,k,l,dislnum)
          end if
      end do
      end do
      end do
      end do
      
      end subroutine updateDislPosMovingMesh
************************************************************************
      subroutine updateDislPosMovingMeshSub(mnumfe,isys,iplane,
     &                                      iobj,dislnum)

C     Inputs: mnumfe --- number of fe material
C             isys --- index of slip system
C             iplane --- index of slip plane within slip system
C             iobj --- index of dislocation within slip plane
C             dislnum --- number of dislocation

C     Outputs: None

C     Purpose: Find new element and local position for a single dislocation
C     If the dislocation moves into the atomistic region, it should
C     naturally be replaced with atomistic dislocation during interpolation
C     (hopefully).

C     Notes/TODO: I'm not sure how to deal with the case of the dislocation
C     leaving the simulation cell entirely

C     Notes/TODO: Doesn't work for multimaterial systems
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: isys, iplane, iobj
      integer :: dislnum
      
C     local variables
      integer :: element, elementnew, mnumfenew
      integer :: bcut, bsgn
      real(dp) :: dislposnew(2)
      real(dp) :: r, s
      logical :: badflip
      
      element = disl(mnumfe)%list(dislnum)%element
      bcut = disl(mnumfe)%list(dislnum)%cut 
      dislposnew = disl(mnumfe)%list(dislnum)%posn
      bsgn = disl(mnumfe)%list(dislnum)%sgn
      
C     figure out where the new dislocation is
      mnumfenew = mnumfe
      elementnew = element
      call findInAllWithGuessDef(dislposnew(1),dislposnew(2),
     &                           mnumfenew,elementnew,r,s,badflip)  
     
      if (badflip) then ! not in mesh. Either in atomistic region or left simulation cell. The second possibility is assumed not to occur
          write(*,*) 'Dislocation escaped while moving mesh'
          write(*,*) 'x', dislposnew(1), 'y', dislposnew(2)
          call deleteDislocation(mnumfe,isys,iplane,iobj)
      else ! still in mesh
          if (mnumfenew/=mnumfe) then ! moved to different material
              call addDislocation(mnumfenew,elementnew,dislposnew(1),
     &                            dislposnew(2),isys,bsgn,bcut) ! assumes slip systems have same numbers (isys) in different materials
              call deleteDislocation(mnumfe,isys,iplane,iobj)
          else ! in the same material; simply update position
              disl(mnumfe)%list(dislnum)%posn = dislposnew
              disl(mnumfe)%list(dislnum)%element = elementnew
              disl(mnumfe)%list(dislnum)%localpos = [r,s]
          end if
      end if
      
      end subroutine updateDislPosMovingMeshSub
************************************************************************
      subroutine updateSourcePosMovingMesh()

C     Inputs: None

C     Outputs: None

C     Purpose: Find new element and local position for all disl. sources

C     Notes/TODO: Does not work for multimaterial systems

      implicit none
      
C     local variables
      integer :: i, j
      integer :: element, mnumfenew
      real(dp) :: x, y
      real(dp) :: r, s
      logical :: badflip
      
      do i = 1, size(sources)
      do j = 1, size(sources(i)%list)
          element = sources(i)%list(j)%element
          x = sources(i)%list(j)%posn(1)
          y = sources(i)%list(j)%posn(2)
          mnumfenew = i
          call findInAllWithGuessDef(x,y,mnumfenew,element,r,s,badflip)
          if (badflip) then
              write(*,*) 'Could not place source within mesh'
              stop
          end if    
          if (mnumfenew /= i) then
              write(*,*) 'Attempting to move source between materials'
              write(*,*) 'Multimaterial extension is not supported'
              stop
          end if
          sources(i)%list(j)%element = element
          sources(i)%list(j)%localpos = [r,s]
      end do
      end do
      
      end subroutine updateSourcePosMovingMesh
************************************************************************
      subroutine updateObstaclePosMovingMesh()

C     Inputs: None

C     Outputs: None

C     Purpose: Find new element and local position for all disl. obstacles

C     Notes/TODO: Does not work for multimaterial systems

      implicit none
      
C     local variables
      integer :: i, j
      integer :: element, mnumfenew
      real(dp) :: x, y
      real(dp) :: r, s
      logical :: badflip
      
      do i = 1, size(obstacles)
      do j = 1, size(obstacles(i)%list)
          element = obstacles(i)%list(j)%element
          x = obstacles(i)%list(j)%posn(1)
          y = obstacles(i)%list(j)%posn(2)
          mnumfenew = i
          call findInAllWithGuessDef(x,y,mnumfenew,element,r,s,badflip)
          if (badflip) then
              write(*,*) 'Could not place obstacle within mesh'
              stop
          end if    
          if (mnumfenew /= i) then
              write(*,*) 'Attempting to move obstacle between materials'
              write(*,*) 'Multimaterial extension is not supported'
              stop
          end if
          obstacles(i)%list(j)%element = element
          obstacles(i)%list(j)%localpos = [r,s]
      end do
      end do
      
      end subroutine updateObstaclePosMovingMesh
************************************************************************
      end module mod_moving_mesh_crack_cadd