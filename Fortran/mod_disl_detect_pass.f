      module mod_disl_detect_pass

C     Purpose: Reads/writes/stores information related to dislocation
C     detection and passing. Also contains procedures for detecting
C     dislocations (using procedure in Stukowski, 2014, JMPS) and passing dislocations
C     (both from the atomistic region to the continuum and vice versa).

C     Long note: There are a number of tricky algorithmic issues involved with dislocation detection and passing
C     1) Branch cut: Dislocation detection algorithms (including the original one in Shilkrot et al., 2004)
C     are unable to determine the branch cut of a dislocation. This is problem-specific information;
C     it is difficult (impossible?) to write a general algorithm.
C     For instance, for a crack problem, the cut of an emitted dislocation originates from the crack tip.

C     2) Direction of travel: Stukowski's algorithm is unable to determine the direction of travel of the dislocation.
C     (The algorithm in Shilkrot can do this, with some difficulty, but suffers from a number of other issues.)
C     What we need is an orientation to the detection band --- a way of specifying "outside" (closer to the continuum)
C     vs. "inside" (completely within the atoms). This is the purpose of the region function (insideDetectionBand_ptr).

C     3) There are two "corner cases" that need to be considered carefully.

C     First, when a dislocation apparently crosses the interface, placing it in the pure atomistic region, 
C     past the detection band, may not be possible. This can occur if the dislocation path cuts a corner
C     of the atomistic region and reenters the continuum region. In this case, I have chosen to simply pass
C     the dislocation through to the other side of the continuum.
C     (Placing it in the atomistic region is problematic --- how do we detect when it subsequently leaves?)

C     Second, the dislocation path during passing should not be allowed to cross over empty space
C     (such as the crack plane). This possibility is prevented by the user
C     defining an interface ("impermissibleedges") that cannot be crossed,
C     and checking if the proposed dislocation path crosses it. The interface
C     should be in the atomistic region.

C     The final logic (in passContinuumtoAtomistic and passAtomistictoContinuum)
C     ends up being fairly involved but (hopefully) robust.

C     TODO: Needs to be modified (heavily?) for 3D.      
      
      use mod_types, only: dp
      use mod_utils, only: readMatTransposeSize, writeMatTransposeSize,
     &                     readVecSize, writeVecSize
      use mod_delaunay, only: delaunay, getTriCenter, genDelaunay
      use mod_disl_ident_simple, only: computeCircuits, identsimple
      use mod_utils, only: writeMat, prettyPrintMat
      use mod_math, only: getIntersectionTwoLines, normalizeVec,
     &                    TOLCONST, sameSign
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
      use mod_neighbors, only: getAtomsInBoxGroupTemp, neighbors
      use mod_misc, only: misc
      use mod_damping, only: dampingdata, readDampingDataSub,
     &                       writeDampingDataSub
      implicit none
 
      private
      public :: initDetectionData, readDetectionData, placeDetectionSub,
     &  processDetectionData, writeDetectionData, assignDetectionPoints,
     &  assignDetectionBand, getDislBranchCut, insideAnnulus,
     &  getDislPropsFromBurgersVec, detectAndPassDislocations,
     &  passContinuumToAtomistic, passAtomistictoContinuum, detection,
     &  updateAtomsPassing, insideRectAnnulus, imposeDipoleDispOnAtoms,
     &  errorInterface, placeInsideDetection, placeOutsideDetection,
     &  findInterfaceIntersectionDeformed, BOXWIDTHNORM,
     &  getPaddedParamsAnnulus, getPaddedParamsRectAnnulus
     
      type compdata
      real(dp), allocatable :: burgersvec(:)
      integer :: isys
      integer :: bsgn
      end type

      type detectiondata
C     (read-in)
      character(len=20) :: bandtype
      type(dampingdata) :: damp
      integer, allocatable :: impermissibleedges(:,:)
      real(dp) :: maxdisttointerface
      integer :: mdnincrements
      real(dp) :: mdtimestep
      integer :: mnumfe
      real(dp), allocatable :: params(:)
      real(dp) :: passdistanceatoc
      real(dp) :: passdistancectoa
C     (processed)
      real(dp) :: burgers
      character(len=20) :: lattice
      type(compdata), allocatable :: comp(:)
      real(dp), allocatable :: paramspadded(:)
      end type

C     module variables (private)
      type(detectiondata) :: detection
      procedure(Dummy), pointer :: insideDetectionBand_ptr    

C     HARD-CODED CONSTANTS
      real(dp), parameter :: BOXWIDTHNORM = 7.5_dp ! see updateAtomsPassing
      real(dp), parameter :: DISTPAD = 1.5_dp ! see paramspadded/processDetectionData
      real(dp), parameter :: STEPFAC = 0.01_dp ! passing constants (not too important)
      real(dp), parameter :: FUDGEFAC = 0.02_dp ! passing constants (not too important)
      real(dp), parameter :: CIRCUMSQFACHEX = 2.0_dp ! see identifyLargeTri/processDetectionData
C                                                    ! ratio of circumradius**2 in largest dislocated triangle to circumradius**2 in equilibrium triangle
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
      call readMatTransposeSize(iunit,detection%impermissibleedges)
      read(iunit,*) detection%maxdisttointerface
      read(iunit,*) detection%mdnincrements
      read(iunit,*) detection%mdtimestep
      read(iunit,*) detection%mnumfe
      call readVecSize(iunit,detection%params)
      read(iunit,*) detection%passdistanceatoc
      read(iunit,*) detection%passdistancectoa
      
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
      
      implicit none
      
C     local variables
      integer :: mnumfe, mnum
      integer :: isys, bsgn
      integer :: nslipsys, ncomp, counter
      real(dp) :: direction(2)
      real(dp) :: padding
      
      mnumfe = detection%mnumfe
      mnum = fematerials%list(mnumfe)

C     identsimple stuff
      identsimple%mnum = mnum
      
C     detection stuff
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
     &         CIRCUMSQFACHEX*(detection%burgers)**2/3.0_dp ! 1/sqrt(3) is factor for circumradius for equilateral triangle
      end if
      
C     delaunay stuff
      padding = DISTPAD*detection%burgers
      call assignDetectionBand(padding)
      
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
      call writeMatTransposeSize(iunit,detection%impermissibleedges)
      write(iunit,*) detection%maxdisttointerface
      write(iunit,*) detection%mdnincrements
      write(iunit,*) detection%mdtimestep
      write(iunit,*) detection%mnumfe
      call writeVecSize(iunit,detection%params) 
      write(iunit,*) detection%passdistanceatoc
      write(iunit,*) detection%passdistancectoa
      
      close(iunit)
      
      end subroutine writeDetectionData
************************************************************************
      subroutine assignDetectionPoints()
 
C     Subroutine: assignDetectionPoints                                                                                                                                                                                                                                                                                                                                                                                                                                     ctionPoints
 
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
      real(dp) :: tempxy(2,nodes%nrealatoms)
      integer :: tempnodes(nodes%nrealatoms)
      real(dp) :: posn(2)
      
      counter = 0
      do i = 1, nodes%nrealatoms
          node = nodes%realatomlist(i)
          posn = nodes%posn(1:2,node) ! current position
          if (insideDetectionBand_ptr(posn,
     &                                detection%paramspadded) == 0) then ! use padded params to grab atoms just outside the band
              counter = counter + 1
              tempxy(:,counter) = posn
              tempnodes(counter) = node
          end if
      end do
      
C     final result, store in delaunay
      if (allocated(delaunay%xy)) then
          deallocate(delaunay%xy)
      end if
      if (allocated(delaunay%nodenums)) then
          deallocate(delaunay%nodenums)
      end if  
      delaunay%xy = tempxy(:,1:counter)
      delaunay%nodenums = tempnodes(1:counter)      
      
      end subroutine assignDetectionPoints
************************************************************************
      subroutine assignDetectionBand(padding)
 
C     Subroutine: assignDetectionBand
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Assigns arguments based on style
C     of detection band. Currently only implemented for circular and
C     rectangular annuli.

      implicit none
      
C     input variables
      real(dp) :: padding
      
      select case (trim(detection%bandtype))
          case ('annulus')
              insideDetectionBand_ptr => insideAnnulus
              detection%paramspadded =
     &              getPaddedParamsAnnulus(detection%params,padding)
          case ('rect_annulus')
              insideDetectionBand_ptr => insideRectAnnulus
              detection%paramspadded =
     &              getPaddedParamsRectAnnulus(detection%params,padding)
          case default
              write(*,*) 'Detection band style has not yet been defined'
              stop
      end select
      
      end subroutine assignDetectionBand
************************************************************************
      function Dummy(posn,params) result(inside)
      
C     input variables
      real(dp) :: posn(2)
      real(dp) :: params(:)
      
C     output variables
      integer :: inside
      
      end function Dummy
************************************************************************
      function insideAnnulus(posn,params) result(inside)
 
C     Subroutine: insideAnnulus
 
C     Inputs: posn --- vector (length 2) of coordinates of point of interest
C             params(1) --- x-coord of center of (both) circles
C             params(2) --- y-coord of center of (both) circles
C             params(3) --- radius**2 of inner circle
C             params(4) --- radius**2 of outer circle    
 
C     Outputs: inside --- integer indicating whether point is:
C                         0 --- within annulus (inside detection band)
C                         1 --- inside inner circle (closer to atomistic)
C                         -1 --- outside outer circle (closer to continuum)
 
C     Purpose: Figure out where point is with respect to annulus (detection band)
         
      implicit none
         
C     input variables
      real(dp) :: posn(2)
      real(dp) :: params(:)
      
C     output variables
      integer :: inside
      
C     local variables
      real(dp) :: dx, dy
      real(dp) :: rsq

      dx = posn(1) - params(1)
      dy = posn(2) - params(2)
      rsq = dx**2 + dy**2
      
      inside = -1
      if (rsq < params(3)) then
          inside = 1
      else if (rsq < params(4)) then
          inside = 0
      end if   
      
      end function insideAnnulus
************************************************************************
      function getPaddedParamsAnnulus(params,padding)
     &                                             result(paramspadded)
 
C     Function: getPaddedParamsAnnulus
 
C     Inputs: params --- parameters for old detection band (annulus). See insideAnnulus, above.
 
C     Outputs: paramspadded --- parameters for new detection band that is slightly larger
 
C     Purpose: Generate parameters for new, "padded" detection band (so that we grab atoms close to edge of band)
      
      implicit none
      
C     input variables
      real(dp) :: params(:)
      real(dp) :: padding
      
C     output variables
      real(dp), allocatable :: paramspadded(:)
      
C     local variables
      real(dp) :: rinner, router
      
      paramspadded = params
      rinner = sqrt(paramspadded(3)) - padding
      router = sqrt(paramspadded(4)) + padding
      paramspadded(3) = rinner**2
      paramspadded(4) = router**2
      
      end function getPaddedParamsAnnulus
************************************************************************
      function insideRectAnnulus(posn,params) result(inside)
 
C     Subroutine: insideRectAnnulus
 
C     Inputs: posn --- vector (length 2) of coordinates of point of interest
C             params(1) --- x-coord of center of (both) rectangles
C             params(2) --- y-coord of center of (both) rectangles
C             params(3) --- Lx/2 for inner rectangle
C             params(4) --- Lx/2 for outer rectangle
C             params(5) --- Ly/2 for inner rectangle
C             params(6) --- Ly/2 for outer rectangle
 
C     Outputs: inside --- integer indicating whether point is:
C                         0 --- within rectangular annulus (within detection band)
C                         1 --- inside inner box (closer to atomistic)
C                         -1 --- outside outer box (closer to continuum)
C                         
C     (Rectangular annulus is like space between two rectangles, assuming symmetry in both x- and y- directions)
C     (Like "rectangular hollow section" beam)
 
C     Purpose: Figure out where point is with respect to rectangular annulus (detection band)
      
      implicit none
      
C     input variables
      real(dp) :: posn(2)
      real(dp) :: params(:)
      
C     output variables
      integer :: inside
      
C     local variables
      real(dp) :: dx, dy
      
      dx = abs(posn(1) - params(1))
      dy = abs(posn(2) - params(2))
      
      inside = -1
      if (dx < params(3)) then
          if (dy < params(5)) then
              inside = 1
              return
          end if
      end if
      
      if (dx < params(4)) then
          if (dy < params(6)) then
              inside = 0
          end if
      end if
      
      end function insideRectAnnulus
************************************************************************
      function getPaddedParamsRectAnnulus(params,padding)
     &                                             result(paramspadded)
 
C     Function: getPaddedParamsRectAnnulus
 
C     Inputs: params --- parameters for old detection band (rect. annulus). See insideRectAnnulus, above.
 
C     Outputs: paramspadded --- parameters for new detection band that is slightly larger
 
C     Purpose: Generate parameters for new, "padded" detection band (so that we grab atoms close to edge of band)
      
      implicit none
      
C     input variables
      real(dp) :: params(:)
      real(dp) :: padding
      
C     output variables
      real(dp), allocatable :: paramspadded(:)
      
      paramspadded = params
      paramspadded(3) = paramspadded(3) - padding
      paramspadded(4) = paramspadded(4) + padding
      paramspadded(5) = paramspadded(5) - padding
      paramspadded(6) = paramspadded(6) + padding
      
      end function getPaddedParamsRectAnnulus
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
      
      bestnorm = huge(0.0_dp)
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
      function getDislBranchCut(mnumfe,isys,posn) result(bcut)
      
C     Function: getDislBranchCut

C     Inputs: mnumfe --- fe material number
C             isys --- number of slip system of dislocation
C             posn --- position of dislocation

C     Outputs: bcut --- branch cut of dislocation

C     Purpose: Determine branch cut of dislocation from attributes of the dislocation detected in the atomistic region
C     (the burgers circuit alone does not tell us this). For a crack problem,
C     we assume that the cut points towards the crack.

      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: isys
      real(dp) :: posn(2)
      
C     output variables
      integer :: bcut
      
      bcut = 0 ! anything else would require problem specific information
               ! this only fails (in 2D) for emitted dislocations that are traveling backwards, which seems unlikely
      
      end function getDislBranchCut
************************************************************************
      subroutine detectAndPassDislocations(detected)
      
C     Subroutine: detectAndPassDislocations
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Loop over all triangles in detection band triangulation,
C     find dislocations (if any), pass them to continuum
 
      implicit none
     
      logical :: detected ! FIX
     
C     local variables
      integer :: i
      integer :: node
      integer :: mnumfe
      real(dp), allocatable :: circuit(:,:)
      real(dp) :: burgersinfnorm
      real(dp) :: burgersvec(2), dislpos(2)
      integer :: bcut, bsgn, isys
    
C     regen delaunay, if necessary (note: neighbors%delaunayregen set to .true. whenever neighbor list is updated)
      if (neighbors%delaunayregen) then
          call assignDetectionPoints() ! reassign nodes/pts
          call genDelaunay() ! regenerate triangulation
      else ! otherwise, just update xy coordinates
          do i = 1, size(delaunay%nodenums)
              node = delaunay%nodenums(i)
              delaunay%xy(:,i) = nodes%posn(1:2,node) ! current position
          end do
      end if     
      
C     compute burgers circuits
      mnumfe = detection%mnumfe
      circuit = computeCircuits()
      detected = .false.
      do i = 1, delaunay%numtri
      if (delaunay%trigood(i)) then
          burgersvec = circuit(:,i)
          burgersinfnorm = maxval(abs(burgersvec))
          if (burgersinfnorm > TOLCONST) then
              dislpos = getTriCenter(i)
              if (insideDetectionBand_ptr(dislpos,
     &                                    detection%params) == 0) then ! dislocation must be inside actual detection band (not padded one)
                  write(*,*) 'Found dislocation'
                  detected = .true.
                  call getDislPropsFromBurgersVec(burgersvec,isys,bsgn)
                  write(*,*) 'Dislocation position', dislpos
                  bcut = getDislBranchCut(mnumfe,isys,dislpos)
                  call passAtomisticToContinuum(dislpos,bsgn,bcut,isys)
              end if
          end if
      end if
      end do
      
      end subroutine detectAndPassDislocations
************************************************************************
      subroutine passAtomisticToContinuum(detectpos,bsgn,bcut,isys)
      
C     Subroutine: passAtomisticToContinuum
 
C     Inputs: detectpos --- position (length 2) of detected dislocation
C             bsgn --- sign of dislocation (+1 or -1)
C             bcut --- branch cut of dislocation (0 if to the left, 1 if to the right)
C             isys --- number of slip system
 
C     Outputs: None
 
C     Purpose: Pass a detected dislocation from the atomistic region to the
C              continuum region
      
      implicit none
      
C     input variables
      real(dp) :: detectpos(2)
      integer :: bsgn
      integer :: bcut
      integer :: isys
      
C     local variables
      integer :: mnumfe
      real(dp) :: dispnorm(2)
      real(dp) :: cost, sint
      real(dp) :: pintdummy(2), pint(2), pintcrack(2)
      logical :: isint
      real(dp) :: trypos(2), posnew(2)
      integer :: edgenum, elguess
      real(dp) :: burgers
      integer, allocatable :: impedges(:,:)
      
C     figure out direction that dislocation is moving by determining where it exits the detection band
      mnumfe = detection%mnumfe
      dispnorm = slipsys(mnumfe)%trig(:,isys)
      cost = dispnorm(1)
      sint = dispnorm(2)
      call placeOutsideDetection(detectpos,dispnorm,pintdummy,isint)
 
C     place end of dislocation path past interface, find intersection with interface
      trypos = detectpos + dispnorm*detection%maxdisttointerface       
      call findInterfaceIntersectionDeformed(interfaceedges%array,
     &                             detectpos,trypos,pint,isint,edgenum)
      if (.not.isint) then ! if no intersection, something went wrong
          call errorInterface(1,detectpos,trypos)
      end if
      mnumfe = interfaceedges%array(3,edgenum)
      elguess = interfaceedges%array(4,edgenum)
      
C     is this path permissible? (e.g. does it cross a crack surface)
      impedges = detection%impermissibleedges
      call findInterfaceIntersectionDeformed(impedges,detectpos,
     &                                     pint,pintcrack,isint,edgenum)
      if (isint) then ! if path crosses crack face, leave dislocation where it is --- it will eventually exit by crossing crack surface
          return
      end if
      
C     otherwise, place dislocation away from interface, by an amount passdistanceatoc
      trypos = pint + dispnorm*detection%passdistanceatoc
      
C     is this path permissible? (e.g. does it cross a crack surface)
      call findInterfaceIntersectionDeformed(impedges,pint,
     &                                   trypos,pintcrack,isint,edgenum)
      if (isint) then
          posnew = 0.5_dp*(pint + pintcrack) ! place disl halfway between interface and crack surface
      else
          posnew = trypos ! otherwise place it normally
      end if
      
C     finally, place dislocation in continuum
      call addDislocation(mnumfe,elguess,posnew(1),posnew(2),
     &                    isys,bsgn,bcut)
      
C     Add ghost dislocation with *opposite* sign at *old* location, inside
C     atomistic region (Must cancel out dislocation in continuum.)
      call addGhostDislocation(detectpos(1),detectpos(2),
     &                         isys,-bsgn,bcut,mnumfe)
      
C     update atom positions, including imposing dipole displacements
C     and minimizing atomic positions near core
      burgers = detection%burgers
      call updateAtomsPassing(detectpos,posnew,
     &                        burgers,bsgn,bcut,cost,sint)
      
      end subroutine passAtomisticToContinuum
************************************************************************
      subroutine passContinuumToAtomistic(posold,pint,mnumfe,isys,
     &                                    iplane,iobj)
      
C     Subroutine: passContinuumToAtomistic

C     Inputs: posold --- position (length 2) of discrete dislocation (before movement)
C             pint --- position (length 2) of intersection of dislocation path
C                      with interface between atomistic/continuum regions
C             mnumfe --- number of fe material
C             isys --- index of slip system
C             iplane --- index of slip plane within slip system
C             iobj --- index of dislocation along slip plane

C     Outputs: None

C     Purpose: Pass a discrete dislocation from the continuum region to the
C              atomistic region
      
      implicit none
      
C     input variables
      real(dp) :: posold(2)
      real(dp) :: pint(2)
      integer :: mnumfe
      integer :: isys, iplane, iobj
      
C     local variables
      integer :: dislnum
      integer :: bsgn, bcut
      real(dp) :: burgers
      real(dp) :: cost, sint
      real(dp) :: disp(2), dispnorm(2)
      real(dp) :: pintdetect(2), pintfudge(2), pint2(2), pintcrack(2)
      real(dp) :: pintcrack2(2), pintdetect2(2), pintactual(2)
      real(dp) :: trypos(2), posnew(2)
      logical :: isint, isintdetect, isintcrack
      logical :: isintcrack2, isintdetect2, isintactual
      real(dp) :: distactual, disttry
      integer :: edgenum
      integer :: elguess
      integer, allocatable :: impedges(:,:)

      dislnum = disl(mnumfe)%splanes(isys)%splane(iplane)%objnum(iobj)
      bsgn = disl(mnumfe)%list(dislnum)%sgn
      bcut = disl(mnumfe)%list(dislnum)%cut
      burgers = detection%burgers
      cost = slipsys(mnumfe)%trig(1,isys)
      sint = slipsys(mnumfe)%trig(2,isys)
      impedges = detection%impermissibleedges
      
C     dislocation travel direction
      disp = pint - posold
      dispnorm = normalizeVec(disp)

C     try to place dislocation inside detection band
      call placeInsideDetection(pint,dispnorm,pintdetect,isintdetect)
      
      if (.not.isintdetect) then
      
C         dislocation is cutting the band obliquely
C         try to find second intersection with interface
          pintfudge = pint + dispnorm*burgers*0.1_dp ! fudge position slightly to avoid multiple intersections
          trypos = pintfudge + dispnorm*detection%maxdisttointerface
          call findInterfaceIntersectionDeformed(interfaceedges%array,
     &                       pintfudge,trypos,pint2,isint,edgenum)
          if (isint) then ! place dislocation in continuum. Need to figure out where, though
              trypos = pint2 + dispnorm*detection%passdistanceatoc
              call findInterfaceIntersectionDeformed(impedges,
     &                        pint2,trypos,pintcrack,isintcrack,edgenum)
              if (isintcrack) then
                  posnew = 0.5_dp*(pint2 + trypos) ! place disl halfway between interface and crack surface
              else
                  posnew = trypos
              end if
C             Finally, move dislocation to new location
C             Easiest to delete old one (see below) and add new one
              mnumfe = interfaceedges%array(3,edgenum)
              elguess = interfaceedges%array(4,edgenum)
              call addDislocation(mnumfe,elguess,posnew(1),
     &                            posnew(2),isys,bsgn,bcut)
          else
              call errorInterface(2,pintfudge,trypos)
          end if          
      else    
C         we need to figure out where in the atomistic region to place disl.
C         first, figure out if we've crossed crack
          call findInterfaceIntersectionDeformed(impedges,
     &                     pint,pintdetect,pintcrack,isintcrack,edgenum)
     
          if (isintcrack) then
C             place dislocation just beyond crack surface, in (hopefully) empty space
              posnew = pintcrack + dispnorm*burgers*FUDGEFAC
          else
C             try to place dislocation past detection band by some amount
              trypos = pintdetect + dispnorm*detection%passdistancectoa
C             does this new path intersect the crack?
              call findInterfaceIntersectionDeformed(impedges,
     &                 pintdetect,trypos,pintcrack2,isintcrack2,edgenum)
              if (isintcrack2) then
                  posnew = pintcrack2 + dispnorm*burgers*FUDGEFAC
              else
C                 does the path intersect the detection band again?
                  call recrossDetection(pintdetect,dispnorm,
     &                                  pintdetect2,isintdetect2)
                  if (isintdetect2) then
                      posnew = 0.5_dp*(pintdetect + pintdetect2) ! place halfway between detection band and other detection intersection
                                                                 ! I'm worried about possible oscillations back and forth, though
                  else
                      posnew = trypos
                  end if
              end if
          end if
          
C         Add ghost dislocation with *same* sign, at *new*
C         location (inside atomistic region)
          call addGhostDislocation(posnew(1),posnew(2),
     &                             isys,bsgn,bcut,mnumfe)
     
      end if
      
      write(*,*) 'Passed dislocation, c -> a'
      write(*,*) 'Old position', posold
      write(*,*) 'New position', posnew
      
C     delete old dislocation
      call deleteDislocation(mnumfe,isys,iplane,iobj)
      
C     update atom positions, including imposing dipole displacements
C     and minimizing atomic positions near core
      call updateAtomsPassing(posold,posnew,burgers,bsgn,bcut,cost,sint)
      
      end subroutine passContinuumToAtomistic
************************************************************************
      subroutine recrossDetection(pint,dispnorm,pintnew,isint)
      
C     Subroutine: recrossDetection

C     Inputs: pint --- old intersection point (between dislocation path and interface between atomistic and continuum region
C             dispnorm --- normalized vector indicating dislocation travel direction
      
C     Outputs: pintnew --- coordinates of new intersection point (between dislocation path and interface between detection band and "inside")
C              isint --- boolean indicating whether this intersection exists

C     Purpose: Find second intersection of dislocation path with detection band, if it exists. (Dislocation starts inside detection band)

      implicit none
      
C     input variables
      real(dp) :: pint(2)
      real(dp) :: dispnorm(2)
      
C     output variables
      real(dp) :: pintnew(2)
      logical :: isint
      
C     local variables
      real(dp) :: dispmax
      logical :: failed
      real(dp) :: step
      integer :: regdesired, regfailed
      
      step = STEPFAC*detection%burgers
      dispmax = detection%passdistancectoa
      regdesired = 0 ! within band
      regfailed = 2 ! fake region
           
      call placeDetectionSub(pint,dispnorm,regdesired,
     &                      regfailed,step,dispmax,pintnew,isint,failed)
      
      end subroutine recrossDetection
************************************************************************
      subroutine placeInsideDetection(pint,dispnorm,pintnew,isint)
      
C     Subroutine: placeInsideDetection

C     Inputs: pint --- old intersection point (between dislocation path and interface between atomistic and continuum region
C             dispnorm --- normalized vector indicating dislocation travel direction
      
C     Outputs: pintnew --- coordinates of new intersection point (between dislocation path and interface between detection band and "inside")
C              isint --- boolean indicating whether this intersection exists

C     Purpose: Find intersection of dislocation path with inner boundary of the detection band

      implicit none
      
C     input variables
      real(dp) :: pint(2)
      real(dp) :: dispnorm(2)
      
C     output variables
      real(dp) :: pintnew(2)
      logical :: isint
      
C     local variables
      real(dp) :: dispmax
      logical :: failed
      real(dp) :: step
      integer :: regdesired, regfailed
      
      step = STEPFAC*detection%burgers
      dispmax = detection%maxdisttointerface
      regdesired = 1 ! inside detection
      regfailed = 2 ! fake region
           
      call placeDetectionSub(pint,dispnorm,regdesired,
     &                      regfailed,step,dispmax,pintnew,isint,failed)
      
      end subroutine placeInsideDetection
************************************************************************
      subroutine placeOutsideDetection(posold,dispnorm,
     &                                 pintnew,isint)
     
C     Subroutine: placeOutsideDetection

C     Inputs: posold --- old position of dislocation
C     
C     In/out: dispnorm --- vector along which dislocation travels (in --- guessed direction; out --- actual direction)
      
C     Outputs: pintnew --- coordinates of new intersection point (between dislocation path and interface between detection band and "outside")
C              isint --- boolean indicating whether this intersection exists

C     Purpose: Find intersection of dislocation path with outer boundary of the detection band

      implicit none
      
C     input variables
      real(dp) :: posold(2)
      
C     in/out variables
      real(dp) :: dispnorm(2)
      
C     output variables
      real(dp) :: pintnew(2)
      logical :: isint
      
C     local variables
      real(dp) :: step, dispmax
      integer :: regdesired, regfailed
      logical :: failed

      step = STEPFAC*detection%burgers
      dispmax = detection%maxdisttointerface
      regdesired = -1 ! outside detection
      regfailed = 1 ! inside detection
      
C     first try moving dislocation in direction disldispnorm      
      call placeDetectionSub(posold,dispnorm,regdesired,regfailed,
     &                       step,dispmax,pintnew,isint,failed)
C     if that failed, move it in the opposite direction
      if (failed) then
          dispnorm = -dispnorm
          call placeDetectionSub(posold,dispnorm,regdesired,regfailed,
     &                      step,dispmax,pintnew,isint,failed)
      end if
      if (failed) then
          write(*,*) 'Could not place dislocation outside detection'
          stop
      end if
      
      end subroutine placeOutsideDetection
************************************************************************
      subroutine placeDetectionSub(posold,dispnorm,regdesired,regfailed,
     &                             step,dispmax,pint,isint,failed)
     
C     Subroutine: placeDetectionSub

C     Inputs: posold --- old position of dislocation
C             dispnorm --- unit vector along which dislocation travels
C             regdesired --- desired region that dislocation ends up in (0 - detection band, 1 - inside, -1 - outside)
C             regfailed --- if disl. path crosses this region, test fails ((0 - detection band, 1 - inside, -1 - outside)
C             step --- step along path in increments of this
C             dispmax --- maximum distance to step
      
C     Outputs: pint --- coordinates of intersection point (where search terminates)
C              isint --- boolean indicating whether this intersection exists
C              failed --- boolean indicating that path hit regfailed first

C     Purpose: Helper for place[Inside/Outside]Detection. Travel along a path,
C     try to hit desired region before hitting failed one; find intersection
C     with desired region boundary.
      
      implicit none
      
C     input variables
      real(dp) :: posold(2)
      real(dp) :: dispnorm(2)
      integer :: regdesired, regfailed
      real(dp) :: step
      real(dp) :: dispmax
      
C     output variables
      real(dp) :: mag
      real(dp) :: pint(2)
      logical :: isint, failed
      integer :: region
      
      isint = .false.
      failed = .false.
      mag = 0.0_dp
      do while (((.not.failed).and.(.not.isint)).and.(mag < dispmax))
          mag = mag + step
          pint = posold + mag*dispnorm
          region = insideDetectionBand_ptr(pint,detection%params)
          isint = (region == regdesired)
          failed = (region == regfailed)
      end do
      
      end subroutine placeDetectionSub
************************************************************************
      subroutine findInterfaceIntersectionDeformed(interfaceedgesarray,
     &                             p0,p1,pint,isint,edgenum)

C     Subroutine: findInterfaceIntersectionDeformed
      
C     Inputs: interfaceedgesarray --- array containing interface edges (either of detection band, or of fe elements)
C             p0, p1 --- coordinates (vector, length 2) of points defining line
      
C     Outputs: pint --- coordinates of intersection point (if it exists) (vector, length 2)
C              isint --- boolean indicating whether intersection exists
C              edgenum --- number of intersected edge

C     Purpose: Finds *closest* intersection between line from p0 and p1 (e.g. a path of a dislocation) and an interface
C     described by the array edgesarray, in the *deformed* coordinate system. If search is successful,
C     return pint (intersection point) and edgenum. "Closest" is defined as being the intersection point (pint)
C     closest to p0.

      implicit none
      
C     input variables
      integer :: interfaceedgesarray(:,:)
      real(dp) :: p0(2), p1(2)
      
C     output variables
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      
C     local variables
      integer :: i
      integer :: node1, node2
      real(dp) :: posn1(2), posn2(2)
      real(dp) :: distsq, distsqnew
      real(dp) :: pintnew(2)
      logical :: isintnew

      distsq = huge(0.0_dp)
      isint = .false.
      
C     loop over interface edges
      do i = 1, size(interfaceedgesarray,2)
          node1 = interfaceedgesarray(1,i)
          node2 = interfaceedgesarray(2,i)
          posn1 = nodes%posn(1:2,node1)
          posn2 = nodes%posn(1:2,node2) 
          call getIntersectionTwoLines(p0,p1,posn1,posn2,
     &                                 pintnew,isintnew)
          if (isintnew) then
              distsqnew = sum((pintnew - p0)**2)
              if (distsqnew < distsq) then
                  edgenum = i
                  distsq = distsqnew
                  pint = pintnew
                  isint = .true.
              end if
          end if    
      end do

      end subroutine findInterfaceIntersectionDeformed
************************************************************************
      subroutine errorInterface(option,dislposold,dislposnew)

C     Subroutine: errorInterface
      
C     Inputs: option --- 1 for atomistic to continuum, 2 for continuum to atomistic
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
          write(*,*) 'Cannot place dislocation in the atomistic region'
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
      real(dp) :: boxwidth

C     create group of atoms for dislocation dipole displacements
C     we only need to get displacements correct for atoms sufficiently close to dipole...
C     (since dipole field decays rapidly in space)
C     so, we define a box around dipole with a "fudge" of boxfudge*burgers
C     (this is a bit hackish)
      boxwidth = BOXWIDTHNORM*burgers
      xmin = min(dislposnew(1),dislposold(1)) - boxwidth
      xmax = max(dislposnew(1),dislposold(1)) + boxwidth
      ymin = min(dislposnew(2),dislposold(2)) - boxwidth
      ymax = max(dislposnew(2),dislposold(2)) + boxwidth      
      call getAtomsInBoxGroupTemp(xmin,xmax,ymin,ymax)
      
C     add dislocation dipole displacements (subtract disp. for dislposold; add for dislposnew)
      call imposeDipoleDispOnAtoms(dislposnew,dislposold,
     &                             bsgn,bcut,cost,sint,tempgroupname)

C     run damped dynamics to relax atoms (neighbor list automatically regenerated)
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