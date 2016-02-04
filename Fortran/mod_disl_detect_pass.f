      module mod_disl_detect_pass

C     Purpose: Reads/writes/stores information related to dislocation
C     detection and passing. Also contains procedures for detecting
C     dislocations (using procedure in Stukowski, 2014, JMPS) and passing dislocations
C     (both from the atomistic region to the continuum and vice versa).
      
C     TODO: check if use of deformed/undeformed positions is correct
C     Problem: Dislocation is found in mesh using mesh_find, which uses
C     undeformed positions, but this module uses deformed positions to track crossings, etc....
C     so, dislocation may be "lost"

C     TODO: Needs to be modified (heavily?) for 3D.      
      
      use mod_types, only: dp
      use mod_utils, only: readMatTransposeSize, writeMatTransposeSize,
     &                     readVecSize, writeVecSize
      use mod_delaunay, only: delaunay, getTriCenter, genDelaunay
      use mod_disl_ident_simple, only: computeCircuits
      use mod_utils, only: writeMat, prettyPrintMat
      use mod_math, only: getIntersectionTwoLines, normalizeVec,
     &                    piconst, tolconst
      use mod_nodes, only: nodes
      use mod_fe_elements, only: interfaceedges, fematerials
      use mod_groups, only: groups, getGroupNum, tempgroupname
      use mod_disl_fields2, only: getDispAtPoint
      use mod_disl_try, only: addDislocation, deleteDislocation, disl
      use mod_disl_escaped, only: addEscapedDislocation
      use mod_disl_ghost, only: addGhostDislocation
      use mod_slip_sys, only: slipsys
      use mod_mesh_find, only: findInAllWithGuess
      use mod_materials, only: materials, nmaterials
      use mod_integrate, only: loopVerlet
      use mod_neighbors, only: getAtomsInBoxGroupTemp
      use mod_misc, only: misc
      use mod_damping, only: dampingdata, readDampingDataSub,
     &                       writeDampingDataSub
      implicit none
 
      private
      public :: initDetectionData, readDetectionData, boxfudge,
     &  processDetectionData, writeDetectionData, assignDetectionPoints,
     &  assignInsideDetectionBand, assignBranchCut, insideAnnulus,
     &  getDislPropsFromBurgersVec, detectAndPassDislocations,
     &  passContinuumToAtomistic, passAtomistictoContinuum, detection,
     &  updateAtomsPassing, insideRectAnnulus, imposeDipoleDispOnAtoms,
     &  findInterfaceIntersectionDeformed, errorInterface,
     &  findInterfaceIntersectionUndeformed
     
      type compdata
      real(dp), allocatable :: burgersvec(:)
      integer :: isys
      integer :: bsgn
      end type

      type detectiondata
C     (read-in)
      character(len=20) :: bandtype
      type(dampingdata) :: damp
      integer, allocatable :: interfaceedges(:,:)
      integer :: mdnincrements
      real(dp) :: mdtimestep
      integer :: mnumfe
      real(dp), allocatable :: params(:)
      real(dp) :: passdistance
C     (processed)
      real(dp) :: burgers
      character(len=20) :: lattice
      type(compdata), allocatable :: comp(:)
      end type

C     module variables (private)
      type(detectiondata) :: detection
      procedure(Dummy), pointer :: insideDetectionBand_ptr

C     HARD-CODED CONSTANTS
      real(dp), parameter :: distfudge = 1.5_dp ! this parameter is not too important, see passAtomistictoContinuum
      real(dp), parameter :: distfudge2 = 3.0_dp ! see passAtomistictoContinuum. If dislocation path crosses interface at p,
                                                 !  dislocation is placed at p + distfudge2*b*m, where m is direction of travel
      real(dp), parameter :: boxfudge = 10.0_dp ! box border width; see updateAtomsPassing
      real(dp), parameter :: circumsqfachex = 2.0_dp ! see identifyLargeTri/processDetectionData
C                                                  ! ratio of circumradius**2 in largest dislocated triangle to circumradius**2 in equilibrium triangle
                                                   ! (must be less than 3, because otherwise edge triangles would be counted as "good")
      contains
************************************************************************
      subroutine initDetectionData(detectionfile)
 
C     Subroutine: initDetectionData
 
C     Inputs: detectionfile --- filename where dislocation detection
C     and passing data is stored
C     (should be something like '[filepref]_detection')
 
C     Outputs: None
 
C     Purpose: Read, initialize data in "detection" structure, which holds
C     information about dislocation detection and passing

      implicit none
      
C     input variables
      character(len=*) :: detectionfile
      
      call readDetectionData(detectionfile)
      call processDetectionData()
      
      end subroutine initDetectionData
************************************************************************
      subroutine readDetectionData(detectionfile)
 
C     Subroutine: readDetectionData
 
C     Inputs: detectionfile --- filename where dislocation detection
C     and passing data is stored
C     (should be something like '[filepref]_detection')
 
C     Outputs: None
 
C     Purpose: Read detection data from file
           
      implicit none
      
C     input variables
      character(len=*) :: detectionfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=detectionfile)
      
      read(iunit,*) detection%bandtype
      call readDampingDataSub(iunit,detection%damp)
      call readMatTransposeSize(iunit,detection%interfaceedges)
      read(iunit,*) detection%mdnincrements
      read(iunit,*) detection%mdtimestep
      read(iunit,*) detection%mnumfe
      call readVecSize(iunit,detection%params)
      read(iunit,*) detection%passdistance
      
      close(iunit)
      
      end subroutine readDetectionData
************************************************************************
      subroutine processDetectionData()
 
C     Subroutine: processDetectionData
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Assign lattice, burgers vector to detection structure,
C     compute all possible burgers vectors, initialize some attributes
C     of delaunay structure
      
C     Notes/TODO: Assumes detection band and adjacent FE continuum material all has same lattice, burgers vectors

C     Notes/TODO: Needs to be reworked for 3D...
      
C     local variables
      integer :: mnumfe, mnum
      integer :: isys, bsgn
      integer :: nslipsys, ncomp, counter
      real(dp) :: direction(2)
      
      mnumfe = detection%mnumfe
      mnum = fematerials%list(mnumfe)
      detection%lattice = materials(mnum)%lattice
      detection%burgers = materials(mnum)%burgers
      if (detection%lattice == 'hex') then
          nslipsys = size(slipsys(mnumfe)%theta)
          ncomp = 2*nslipsys
          allocate(detection%comp(ncomp))
          counter = 0
          do isys = 1, nslipsys
              do bsgn = -1, 1, 2 ! i.e. -1 and 1
                  counter = counter + 1
                  detection%comp(counter)%isys = isys
                  detection%comp(counter)%bsgn = bsgn
                  direction = slipsys(mnumfe)%trig(:,isys)
                  detection%comp(counter)%burgersvec =
     &                                 detection%burgers*bsgn*direction
              end do    
          end do
          delaunay%circumradiussqcutoff =
     &         circumsqfachex*(detection%burgers)**2/3.0_dp ! 1/sqrt(3) is factor for circumradius for equilateral triangle
      end if
      
      call assignInsideDetectionBand()
      delaunay%regen = .true.
      
      end subroutine processDetectionData
************************************************************************
      subroutine writeDetectionData(detectionfile)
 
C     Subroutine: writeDetectionData
 
C     Inputs: detectiontfile --- filename where where dislocation detection
C     and passing data is stored
C     (should be something like '[filepref]_detection')
 
C     Outputs: None
 
C     Purpose: Write dislocation detection and passing data to file (essentially
C     inverse of readDetectionData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: detectionfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=detectionfile)

      write(iunit,*) detection%bandtype  
      call writeDampingDataSub(iunit,detection%damp)
      call writeMatTransposeSize(iunit,detection%interfaceedges)
      write(iunit,*) detection%mdnincrements
      write(iunit,*) detection%mdtimestep
      write(iunit,*) detection%mnumfe
      call writeVecSize(iunit,detection%params) 
      write(iunit,*) detection%passdistance
      
      close(iunit)
      
      end subroutine writeDetectionData
************************************************************************
      subroutine assignDetectionPoints()
 
C     Subroutine: assignDetectionPoints
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Detection band is fixed in space, but points move. So,
C     whenever we regenerate the Delaunay triangulation, we need to reselect
C     the points.

C     Notes: Alternate approach, involving a fixed set of points in detection band,
C     may not work since band may shear into two and lose continuity.
      
      implicit none
      
C     local variables
      integer :: i
      integer :: node, counter
      real(dp) :: temp(2,nodes%nrealatoms)
      real(dp) :: posn(2)
      
      counter = 0
      do i = 1, nodes%nrealatoms
          node = nodes%realatomlist(i)
          posn = nodes%posn(1:2,node) ! current position
          if (insideDetectionBand_ptr(posn)) then
              counter = counter + 1
              temp(:,counter) = posn
          end if
      end do
      
C     final result, store in delaunay
      if (allocated(delaunay%xy)) then
          deallocate(delaunay%xy)
      end if    
      delaunay%xy = temp(:,1:counter)
      
      end subroutine assignDetectionPoints
************************************************************************
      subroutine assignInsideDetectionBand()
 
C     Subroutine: assignInsideDetectionBand
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Assigns insideDetectionBand_ptr to function, based on style
C     of detection band. Currently only implemented for circular and
C     rectangular annuli.
      
      select case (detection%bandtype)
          case ('annulus')
              insideDetectionBand_ptr => insideAnnulus
          case ('rect_annulus')
              insideDetectionBand_ptr => insideRectAnnulus
          case default
              write(*,*) 'Detection band style has not yet been defined'
              stop
      end select
      
      end subroutine assignInsideDetectionBand
************************************************************************
      function Dummy(posn) result(inside)
      
C     input variables
      real(dp) :: posn(2)
      
C     output variables
      logical :: inside
      
      end function Dummy
************************************************************************
      function insideAnnulus(posn) result(inside)
 
C     Subroutine: insideAnnulus
 
C     Inputs: posn --- vector (length 2) of coordinates of point of interest
 
C     Outputs: inside --- logical indicating whether point is inside annulus
 
C     Purpose: Figure out if point is inside annulus with parameters:
C     params(1) --- x-coord of center of (both) circles
C     params(2) --- y-coord of center of (both) circles
C     params(3) --- radius**2 of inner circle
C     params(4) --- radius**2 of outer circle      
          
C     input variables
      real(dp) :: posn(2)
      
C     output variables
      logical :: inside
      
C     local variables
      real(dp) :: dx, dy
      real(dp) :: rsq

      dx = posn(1) - detection%params(1)
      dy = posn(2) - detection%params(2)
      rsq = dx**2 + dy**2
      
      inside = .false.
      if (rsq >= detection%params(3)) then
          inside = (rsq <= detection%params(4))
      end if    
      
      end function insideAnnulus
************************************************************************
      function insideRectAnnulus(posn) result(inside)
 
C     Subroutine: insideRectAnnulus
 
C     Inputs: posn --- vector (length 2) of coordinates of point of interest
 
C     Outputs: inside --- logical indicating whether point is inside rectangular annulus
C     (space between two rectangles, assuming symmetry in both x- and y- directions)
C     (Like "rectangular hollow section" beam)
 
C     Purpose: Figure out if point is inside rectangular annulus with parameters:
C     params(1) --- x-coord of center of (both) rectangles
C     params(2) --- y-coord of center of (both) rectangles
C     params(3) --- Lx/2 for inner rectangle
C     params(4) --- Lx/2 for outer rectangle
C     params(5) --- Ly/2 for inner rectangle
C     params(6) --- Ly/2 for outer rectangle     
      
C     input variables
      real(dp) :: posn(2)
      
C     output variables
      logical :: inside
      
C     local variables
      real(dp) :: dx, dy
      
      ! params are: xcenter, ycenter, xinner, xouter, yinner, youter
      ! (assumes box is symmetric about center)
      dx = abs(posn(1) - detection%params(1))
      dy = abs(posn(2) - detection%params(2))
      
      inside = .false.
      if (dx > detection%params(3)) then
          if (dx < detection%params(4)) then
              inside = (dy < detection%params(6))
          end if    
      else
          if (dy > detection%params(5)) then
              inside = (dy < detection%params(6))
          end if
      end if
      
      end function insideRectAnnulus
************************************************************************
      function assignBranchCut(burgersvec,isys,bsgn,posn) result(bcut)
      
C     Notes/TODO: Need to fix for 3D...     
      
C     input variables
      real(dp) :: burgersvec(2)
      integer :: isys
      integer :: bsgn
      real(dp) :: posn(2)
      
C     output variables
      integer :: bcut
      
      if ((misc%iscrackproblem).and.
     &    (detection%lattice=='hex')) then
          if (posn(1) > 0.0_dp) then ! cut should be pointed towards crack tip, at (0,0)
              bcut = 0 
          else
              bcut = 1 ! reverse cut (equivalent to 180 degree rotation and sign flip)
          end if
      end if    
      
      end function assignBranchCut
************************************************************************
      subroutine getDislPropsFromBurgersVec(burgersvec,isys,bsgn)
 
C     Subroutine: getDislPropsFromBurgersVec
 
C     Inputs: burgersvec --- burgers vector for dislocation (i.e., found using burgers circuit)
 
C     Outputs: isys --- number of slip system
C              bsgn --- sign of dislocation
 
C     Purpose: Determines slip system, dislocation sign corresponding to a burgers vector
      
      implicit none
      
C     input variables
      real(dp) :: burgersvec(:)
      
C     output variables
      integer :: isys
      integer :: bsgn
      
C     local variables
      real(dp) :: bestnorm, norm
      integer :: i
      
      bestnorm = huge(dp)
      do i = 1, size(detection%comp)
          norm = maxval(
     &            abs((burgersvec - detection%comp(i)%burgersvec)))
          if (norm < bestnorm) then
              bestnorm = norm
              isys = detection%comp(i)%isys
              bsgn = detection%comp(i)%bsgn
          end if    
      end do
      
      end subroutine getDislPropsFromBurgersVec
************************************************************************
      subroutine detectAndPassDislocations()
      
C     Subroutine: detectAndPassDislocations
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Loop over all triangles in detection band triangulation,
C     find dislocations (if any), pass them to continuum
 
C     Notes/TODO: Assumes dislocation travels in direction of positive b...
 
      implicit none
     
C     local variables
      integer :: i
      integer :: mnumfe
      real(dp), allocatable :: circuit(:,:)
      real(dp) :: burgersinfnorm
      real(dp) :: burgersvec(2), disldisp(2), dislpos(2)
      integer :: bcut, bsgn, isys
      
C     regen delaunay, if necessary
      if (delaunay%regen) then
          call assignDetectionPoints()
          call genDelaunay()
      end if
      
C     compute burgers circuits
      mnumfe = detection%mnumfe
      circuit = computeCircuits()
      do i = 1, delaunay%numtri
          if (delaunay%trigood(i)) then
              burgersvec = circuit(:,i)
              burgersinfnorm = maxval(abs(burgersvec))
              if (burgersinfnorm > tolconst) then
                  call getDislPropsFromBurgersVec(burgersvec,isys,bsgn)
                  dislpos = getTriCenter(i)
                  bcut = assignBranchCut(burgersvec,isys,bsgn,dislpos)
                  disldisp = bsgn*slipsys(mnumfe)%trig(:,isys)
                  call passAtomisticToContinuum(dislpos,disldisp,
     &                                          bsgn,bcut,isys)
                  delaunay%regen = .true. ! since we passed dislocation, we'll need to regenerate delaunay next time
              end if
          end if
      end do
      
      end subroutine detectAndPassDislocations
************************************************************************
      subroutine passAtomisticToContinuum(dislpos,disldisp,bsgn,bcut,
     &                                    isys)
      
C     Subroutine: passAtomisticToContinuum
 
C     Inputs: dislpos --- position (length 2) of detected dislocation
C             disldisp --- vector (length 2) indicating direction of movement of disl.
C             bsgn --- sign of dislocation (+1 or -1)
C             bcut --- branch cut of dislocation (0 if to the left, 1 if to the right)
C             isys --- number of slip system
 
C     Outputs: None
 
C     Purpose: Pass a detected dislocation from the atomistic region to the
C              continuum region

C     Notes/TODO: Are all of these constants (distfudge) necessary?
      
      implicit none
      
C     input variables
      real(dp) :: dislpos(2)
      real(dp) :: disldisp(2)
      integer :: bsgn
      integer :: bcut
      integer :: isys
      
C     local variables
      real(dp), allocatable :: disldispnorm(:)
      real(dp) :: dislposnew(2)
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      integer :: mnumfe
      integer :: element
      real(dp) :: burgers
      real(dp) :: cost, sint
 
C     place end of dislocation path well within continuum
C     tacit assumption that continuum is sufficiently "thick"
C     (so there aren't multiple intersections)
C     if this is an issue, findInterfaceIntersection should be modified
C     to find *nearest* intersection
      disldispnorm = normalizeVec(disldisp)
      dislposnew = dislpos +
     &                    disldispnorm*detection%passdistance*distfudge
      call findInterfaceIntersectionUndeformed(interfaceedges%array,
     &                           dislpos,dislposnew,pint,isint,edgenum)
      mnumfe = interfaceedges%array(3,edgenum)
      element = interfaceedges%array(4,edgenum)
      if (.not.isint) then ! if no intersection, something went wrong
          call errorInterface(1,dislpos,dislposnew)
      end if
      
C     determine actual new position of dislocation
      burgers = detection%burgers
      dislposnew = pint + disldispnorm*distfudge2*burgers
      
C     add discrete dislocation to DD array
      call addDislocation(mnumfe,element,dislposnew(1),dislposnew(2),
     &                    isys,bsgn,bcut)
      
C     Add ghost dislocation with *opposite* sign at *old* location, inside
C     atomistic region (Must cancel out dislocation in continuum.)
      call addGhostDislocation(dislpos(1),dislpos(2),
     &                         isys,-bsgn,bcut,mnumfe)
      
C     update atom positions, including imposing dipole displacements
C     and minimizing atomic positions near core
      cost = slipsys(mnumfe)%trig(1,isys)
      sint = slipsys(mnumfe)%trig(2,isys)
      call updateAtomsPassing(dislpos,dislposnew,
     &                        burgers,bsgn,bcut,cost,sint)
      
      end subroutine passAtomisticToContinuum
************************************************************************
      subroutine passContinuumToAtomistic(dislpos,pint,mnumfe,isys,
     &                                    iplane,iobj)
      
C     Subroutine: passContinuumToAtomistic

C     Inputs: dislpos --- position (length 2) of discrete dislocation (before movement)
C             pint --- position (length 2) of intersection of dislocation path
C                      with interface between atomistic/continuum regions
C             mnumfe --- number of fe material
C             isys --- index of slip system
C             iplane --- index of slip plane within slip system
C             iobj --- index of dislocation along slip plane

C     Outputs: None

C     Purpose: Pass a discrete dislocation from the continuum region to the
C              atomistic region

C     Notes/TODO: Are all of these constants (distfudge) necessary?
      
      implicit none
      
C     input variables
      real(dp) :: dislpos(2)
      real(dp) :: pint(2)
      integer :: mnumfe
      integer :: isys, iplane, iobj
      
C     local variables
      integer :: dislnum
      real(dp) :: disldisp(2)
      real(dp), allocatable :: disldispnorm(:)
      real(dp) :: dislposnew(2)
      real(dp) :: pintnew(2)
      real(dp) :: burgers
      integer :: bsgn, bcut
      real(dp) :: cost, sint
      logical :: isint
      integer :: edgenum
    
      dislnum = disl(mnumfe)%splanes(isys)%splane(iplane)%objnum(iobj)
      bsgn = disl(mnumfe)%list(dislnum)%sgn
      bcut = disl(mnumfe)%list(dislnum)%cut
      cost = slipsys(mnumfe)%trig(1,isys)
      sint = slipsys(mnumfe)%trig(2,isys)

C     determine intersection of disl. path with detection band      
      disldisp = pint - dislpos
      disldispnorm = normalizeVec(disldisp)
      dislposnew = pint + disldispnorm*detection%passdistance*distfudge
      call findInterfaceIntersectionDeformed(detection%interfaceedges,
     &                           pint,dislposnew,pintnew,isint,edgenum)
      if (.not.isint) then ! if no intersection, something went wrong
          call errorInterface(2,pint,dislposnew)
      end if    
      
C     place dislocation just inside detection band
      burgers = detection%burgers
      dislposnew = pintnew + disldispnorm*distfudge2*burgers
      
C     delete discrete dislocation from DD array
      call deleteDislocation(mnumfe,isys,iplane,iobj)
      
C     Add ghost dislocation with *same* sign, at *new*
C     location (inside atomistic region)
      call addGhostDislocation(dislposnew(1),dislposnew(2),
     &                         isys,bsgn,bcut,mnumfe)
      
C     update atom positions, including imposing dipole displacements
C     and minimizing atomic positions near core
      call updateAtomsPassing(dislpos,dislposnew,
     &                        burgers,bsgn,bcut,cost,sint)
      
      end subroutine passContinuumToAtomistic
************************************************************************
      subroutine findInterfaceIntersectionDeformed(interfaceedgesarray,
     &                                         p0,p1,pint,isint,edgenum)

C     Subroutine: findInterfaceIntersectionDeformed
      
C     Inputs: interfaceedgesarray --- array containing interface edges (either of detection band, or of fe elements)
C             p0, p1 --- coordinates (vector, length 2) of points defining line
      
C     Outputs: pint --- coordinates of intersection point (if it exists) (vector, length 2)
C              isint --- boolean indicating whether intersection exists
C              edgenum --- number of intersected edge

C     Purpose: Find where a line between p0 and p1 (i.e. a path of a dislocation) intersects
C     an edge in the array edgesarray, using *deformed* coordinates. If search is successful,
C     return pint (intersection point) and edgenum
      
C     Notes: Assumes multiple intersections are not possible...

      implicit none
      
C     input variables
      integer :: interfaceedgesarray(:,:)
      real(dp) :: p0(2), p1(2)
      
C     output variables
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      
      call findInterfaceIntersectionSub(interfaceedgesarray,
     &                                p0,p1,pint,isint,edgenum,.false.)

      end subroutine findInterfaceIntersectionDeformed
************************************************************************
      subroutine findInterfaceIntersectionUndeformed(
     &                     interfaceedgesarray,p0,p1,pint,isint,edgenum)

C     Subroutine: findInterfaceIntersectionDeformed
      
C     Inputs: interfaceedgesarray --- array containing interface edges (either of detection band, or of fe elements)
C             p0, p1 --- coordinates (vector, length 2) of points defining line
      
C     Outputs: pint --- coordinates of intersection point (if it exists) (vector, length 2)
C              isint --- boolean indicating whether intersection exists
C              edgenum --- number of intersected edge

C     Purpose: Find where a line between p0 and p1 (i.e. a path of a dislocation) intersects
C     an edge in the array edgesarray, using *deformed* coordinates. If search is successful,
C     return pint (intersection point) and edgenum
      
C     Notes: Assumes multiple intersections are not possible...

      implicit none
      
C     input variables
      integer :: interfaceedgesarray(:,:)
      real(dp) :: p0(2), p1(2)
      
C     output variables
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      
      call findInterfaceIntersectionSub(interfaceedgesarray,
     &                                p0,p1,pint,isint,edgenum,.true.)

      end subroutine findInterfaceIntersectionUndeformed
************************************************************************
      subroutine findInterfaceIntersectionSub(interfaceedgesarray,
     &                             p0,p1,pint,isint,edgenum,undeformed)

C     Subroutine: findInterfaceIntersectionSub
      
C     Inputs: interfaceedgesarray --- array containing interface edges (either of detection band, or of fe elements)
C             p0, p1 --- coordinates (vector, length 2) of points defining line
C             undeformed --- flag indicating whether undeformed coordinates should be used
      
C     Outputs: pint --- coordinates of intersection point (if it exists) (vector, length 2)
C              isint --- boolean indicating whether intersection exists
C              edgenum --- number of intersected edge

C     Purpose: Find where a line between p0 and p1 (i.e. a path of a dislocation) intersects
C     an edge in the array edgesarray, using *undeformed coordinates*. If search is successful,
C     return pint (intersection point) and edgenum
      
C     Notes: Assumes multiple intersections are not possible...

      implicit none
      
C     input variables
      integer :: interfaceedgesarray(:,:)
      real(dp) :: p0(2), p1(2)
      logical :: undeformed
      
C     output variables
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      
C     local variables
      integer :: i
      integer :: node1, node2
      real(dp) :: posn1(2), posn2(2)

C     loop over interface edges
      do i = 1, size(interfaceedgesarray,2)
          node1 = interfaceedgesarray(1,i)
          node2 = interfaceedgesarray(2,i)
          posn1 = nodes%posn(1:2,node1)
          posn2 = nodes%posn(1:2,node2)
          if (undeformed) then
              posn1 = posn1 - nodes%posn(4:5,node1)
              posn2 = posn2 - nodes%posn(4:5,node2)
          end if    
          call getIntersectionTwoLines(p0,p1,posn1,posn2,pint,isint)
          if (isint) then
              edgenum = i
              return
          end if    
      end do
      
C     if we have gotten here, no intersection was found
      isint = .false.

      end subroutine findInterfaceIntersectionSub
************************************************************************
      subroutine errorInterface(option,dislposold,dislposnew)

C     Subroutine: errorInterface
      
C     Inputs: option --- 1 for actual interface, 2 for detection band interface
C             dislposold --- old position of dislocation
C             dislposnew --- new position of dislocation (having apparently crossed interface)
      
C     Outputs: None

C     Purpose: Print (helpful) error message indicating that
C     intersection of dislocation path with interface (either actual interface,
C     or inner boundary of detection band) was not found
      
      implicit none
      
C     input variables
      integer :: option
      real(dp) :: dislposold(2)
      real(dp) :: dislposnew(2)
      
      if (option == 1) then
          write(*,*) 'Dislocation intersection with interface'
          write(*,*) 'between atomistic and FE regions'
      else if (option == 2) then
          write(*,*) 'Dislocation intersection with detection band'
          write(*,*) 'interface'
      end if    
      write(*,*) 'could not be found.'
      write(*,*) 'Disl. path from', dislposold
      write(*,*) 'to', dislposnew
      stop
      
      end subroutine errorInterface
************************************************************************
      subroutine updateAtomsPassing(dislposold,dislposnew,burgers,
     &                              bsgn,bcut,cost,sint)

C     Subroutine: updateAtomsPassing

C     Inputs: dislposold --- old position (length 2) of dislocation
C             dislposnew --- new position (length 2) of dislocation
C             burgers --- magnitude of burgers vector
C             bsgn --- sign of dislocation (+1 or -1)
C             bcut --- branch cut of dislocation (0 if to the left, 1 if to the right)
C             cost, sint --- orientation of dislocation (radians)

C     Outputs: None

C     Purpose: Update atomic positions as a result of dislocation being passed
C     (in either direction). Helper for passAtomistictoContinuum and
C     passContinuumtoAtomistic

      implicit none
      
C     input variables
      real(dp) :: dislposold(2)
      real(dp) :: dislposnew(2)
      real(dp) :: burgers
      integer :: bsgn, bcut
      real(dp) :: cost, sint
      
C     local variables
      real(dp) :: xmin, xmax, ymin, ymax
      logical :: dampflag

C     create group of atoms for dislocation dipole displacements
C     we only need to get displacements correct for atoms sufficiently close to dipole...
C     (since dipole field decays rapidly in space)
C     so, we define a box around dipole with a "fudge" of boxfudge*burgers
C     (this is a bit hackish),
      xmin = min(dislposnew(1),dislposold(1)) - boxfudge*burgers
      xmax = max(dislposnew(1),dislposold(1)) + boxfudge*burgers
      ymin = min(dislposnew(2),dislposold(2)) - boxfudge*burgers
      ymax = max(dislposnew(2),dislposold(2)) + boxfudge*burgers          
      call getAtomsInBoxGroupTemp(xmin,xmax,ymin,ymax)
      
C     add dislocation dipole displacements (subtract disp. for dislposold; add for dislposnew)
      call imposeDipoleDispOnAtoms(dislposnew,dislposold,
     &                             bsgn,bcut,cost,sint,tempgroupname)

C     run damped dynamics to relax atoms (neighbor list automatically regenerated)
      dampflag = .true.
      call loopVerlet(detection%mdnincrements,detection%mdtimestep,
     &                tempgroupname,detection%damp)
      
      end subroutine updateAtomsPassing
************************************************************************
      subroutine imposeDipoleDispOnAtoms(dislpos1,dislpos2,
     &                                   bsgn1,bcut1,cost,sint,gname)
     
C     Subroutine: imposeDipoleDispOnAtoms

C     Inputs: dislpos1 --- position of first dislocation (vector, length 2)
C             dislpos2 --- position of second dislocation (vector, length 2)
C             bsgn1 --- sign of burgers vector of first dislocation (+1 or -1)
C                       (sign of second dislocation is -bsgn1, since we assume a dipole)
C             bcut1 --- branch cut of first dislocation
C             cost, sint --- orientation of dislocation
C             gname --- name of group for which to apply displacements to atoms (see mod_groups)
      
C     Outputs: None

C     Purpose: Apply displacements of a dislocation dipole to all atoms in a group (see mod_groups)
C              *including* pad atoms

C     Notes/TODO: Tacitly assumes no elastic mismatch in atomistic region
      
      implicit none
      
C     input variables
      real(dp) :: dislpos1(2), dislpos2(2)
      integer :: bsgn1, bcut1
      real(dp) :: cost, sint
      character(len=*) :: gname
      
C     local variables
      integer :: i
      integer :: node
      integer :: mnum
      integer :: gnum
      real(dp) :: atompos(2)
      real(dp) :: disp(2), disp2(2)
      
      gnum = getGroupNum(gname)
      do i = 1, nodes%natoms
      if (groups(gnum)%maskatoms(i)) then
          node = nodes%atomlist(i)
          atompos = nodes%posn(1:2,node)
          mnum = nodes%types(1,node) ! this tacitly assumes no elastic mismatch in atomistic region
          disp = getDispAtPoint(atompos,dislpos1,
     &                          cost,sint,bsgn1,bcut1,mnum,1.0_dp) ! no fudge
          disp2 = getDispAtPoint(atompos,dislpos2,
     &                           cost,sint,-bsgn1,bcut1,mnum,1.0_dp) ! no fudge
          disp = disp + disp2
          nodes%posn(1:2,node) = nodes%posn(1:2,node) + disp
          nodes%posn(4:5,node) = nodes%posn(4:5,node) + disp
      end if
      end do
      
      end subroutine imposeDipoleDispOnAtoms
************************************************************************     
      
      end module