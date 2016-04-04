      module mod_moving_mesh_crack_cadd
      
C     Purpose: Implements moving mesh technique for CADD crack problem, similar
C     to that outlined in Chakravarthy and Curtin, 2011, MSMSE for the DD crack problem.
C     In the moving mesh method, the problem is shifted so that a crack at (xcrack, 0)
C     is moved to (0, 0). To do this, displacements are found at nodes
C     in the original mesh by finding the corresponding displacement
C     at a point shifted by (xcrack, 0) using interpolation.
C     Also, dislocations are shifted by the same amount. In the end,
C     the result is that all of the field variables have been shifted by xcrack,
C     the crack is back at (xsmall, 0), where xsmall is very close to 0,
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

C     Issue #3 can probably be resolved with a better algorithm,
C     where dislocations are passed out of the simulation cell to become
C     escaped dislocations, and sources/obstacles are deactivated when
C     they leave the simulation cell (although this might cause unpinning
C     of dislocation pileups, etc.)

      use mod_types, only: dp
      use mod_nodes, only: nodes, getXYAtomBoundsDef
      use mod_misc, only: misc
      use mod_fe_el_2d, only: felib, getElTypeNum
      use mod_fe_elements, only: feelements
      use mod_fe_main_2d_assign, only: getTotalDispAtPoint_ptr
      use mod_math, only: nearestMultiple, checkSameSide
      use mod_disl_try, only: disl, sources, obstacles,
     &                        addDislocation, deleteDislocation
      use mod_disl_escaped, only: escapeddisl
      use mod_disl_ghost, only: ghostdisl
      use mod_slip_sys, only: slipsys
      use mod_neighbors, only: genAtomBinsUndeformed, neighbors,
     &                         getAtomBin, updateNeighborsNoCheck
      use mod_delaunay, only: delaunaydata, genDelaunay
      use mod_mesh_find, only: findInAllWithGuessDef, getLocalCoords,
     & findInAllInitiallyUndef, findInAllWithGuessUndef
      use mod_disl_detect_pass, only: detectAndPassDislocations,
     & detection, processDetectionBand
      implicit none
      
      type movingmeshdata
C     read-in
      logical :: active
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
      public :: moveMeshCrack, shiftPosnDisp, initCADDMovingMeshData,
     &  interpFromFEPoint, interpFromAtomPoint, readCADDMovingMeshData,
     &  genAtomTriangulationUndeformed, writeCADDMovingMeshData

      contains
************************************************************************
      subroutine initCADDMovingMeshData(movingmeshfile)
      
      implicit none
      
C     input variables
      character(len=*) :: movingmeshfile
      
      call readCADDMovingMeshData(movingmeshfile)
      
      end subroutine initCADDMovingMeshData
************************************************************************
      subroutine readCADDMovingMeshData(movingmeshfile)
      
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
      
      movingmesh%delaunay%nodenums = nodes%atomlist
      allocate(movingmesh%delaunay%xy(2,nodes%natoms))
      
      end subroutine processCADDMovingMeshData
************************************************************************
      subroutine writeCADDMovingMeshData(movingmeshfile)
      
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
      subroutine moveMeshCrack(crackpos)
      
C     input variables
      real(dp) :: crackpos(2)
      
      movingmesh%xshift = -nearestMultiple(crackpos(1),
     &                                     movingmesh%deltaxcrack) ! shift must be integer multiple of deltaxcrack
                                                                   ! shift is backwards (hence negative sign)
      movingmesh%eltypenum = getElTypeNum('CPE3') ! triangle
      call passDislocationsBefore()    
      call shiftPosnDisp()
      call shiftDDPos()
      call shiftSlipPlanes()
      call updateDislPosMovingMesh()
      call updateSourcePosMovingMesh()
      call updateObstaclePosMovingMesh()
      call passDislocationsAfter()
      call updateNeighborsNoCheck() ! force update of neighbors
      
      end subroutine moveMeshCrack
************************************************************************
      subroutine passDislocationsBefore()
      
      call passDislocationsSub('inner')
      
      end subroutine passDislocationsBefore
************************************************************************
      subroutine passDislocationsAfter()
      
      call passDislocationsSub('outer')
      
      end subroutine passDislocationsAfter
************************************************************************
      subroutine passDislocationsSub(option)
      
C     input variables
      character(len=*) :: option
      
C     local variables
      real(dp) :: xmin, xmax, ymin, ymax
      real(dp) :: tempmaxdisttointerface
      real(dp), allocatable :: tempparams(:), tempparamspadded(:)
      real(dp) :: xc, xmin1, xmin2, xminnew, xdist1, xdist2
      real(dp) :: halfLxinner, halfLxouter
      logical :: detected
      
C     regenerate neighbor lists
      call updateNeighborsNoCheck()
      
C     save old
      tempmaxdisttointerface = detection%maxdisttointerface
      tempparams = detection%params
      tempparamspadded = detection%paramspadded
      
C     get atom box
      call getXYAtomBoundsDef(xmin,xmax,ymin,ymax)
      detection%maxdisttointerface =
     &                      sqrt((xmax - xmin)**2 + (ymax - ymin)**2)
      
C     resize detection band
      xc = detection%params(1)
      if (detection%bandtype == 'rect_annulus') then
          if (trim(option) == 'inner') then ! expand dislocation band inwards
              halfLxinner = detection%params(3)
              xmin1 = xc - halfLxinner
              xmin2 = xmin + movingmesh%maxshift
              xminnew = min(xmin1,xmin2)
              halfLxinner = xc - xminnew
              detection%params(3) = halfLxinner
          else if (trim(option) == 'outer') then ! expand dislocation band outwards
              xdist1 = xc - xmin
              xdist2 = xmax - xc
              halfLxouter = max(xdist1,xdist2)
              detection%params(4) = halfLxouter
          end if
      else
          write(*,*) 'Method is undefined for non-rectangular bands'
          stop
      end if 

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
 
C     Subroutine: shiftPosnDisp
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Shift the positions and displacements of all nodes by delta_x = movingmesh%xshift
           
      implicit none
      
C     local variables
      real(dp), allocatable :: newposn(:,:)
      logical, allocatable :: alreadychecked(:)
      integer :: nodesshape(2)
      integer :: i, j, k
      integer :: node
      integer :: m, n
      real(dp) :: posnundef(2), newdisp(2)
      real(dp) :: shiftvec(2)
      
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
      call genAtomTriangulationUndeformed()
      
C     cycle over elements to generate good guesses
      do i = 1, size(feelements)
          do j = 1, size(feelements(i)%connect,2)
          do k = 1, size(feelements(i)%connect,1)
              node = feelements(i)%connect(k,j)
              if (.not.alreadychecked(node)) then
                  posnundef = nodes%posn(1:2,node)-nodes%posn(4:5,node)
                  posnundef = posnundef + shiftvec
                  newdisp = interpFromFEPoint(node,posnundef,i,j)
                  newposn(1:2,i) = posnundef + newdisp
                  newposn(4:5,i) = newdisp
                  alreadychecked(node) = .true.
              end if
          end do
          end do
      end do
      
C     deal with rest of nodes (just atoms)
      do i = 1, nodes%nnodes
      if (.not.alreadychecked(i)) then
          posnundef = nodes%posn(1:2,i) - nodes%posn(4:5,i)
          newdisp = interpFromAtomPoint(posnundef)
          newposn(1:2,i) = posnundef + newdisp
          newposn(4:5,i) = newdisp
      end if
      end do
      
      nodes%posn = newposn
      
      end subroutine shiftPosnDisp
************************************************************************
      subroutine genAtomTriangulationUndeformed()
      
C     local variables
      integer :: i
      integer :: node
      
      do i = 1, nodes%natoms
          node = nodes%atomlist(i)
          movingmesh%delaunay%xy(:,i) = nodes%posn(1:2,node) -
     &                                  nodes%posn(4:5,node)
      end do
      call genDelaunay(movingmesh%delaunay)
      
      end subroutine genAtomTriangulationUndeformed
************************************************************************
      function interpFromFEPoint(node,posnundef,mnumfe,element)
     &                                                     result(disp)

C     input variables
      integer :: node
      real(dp) :: posnundef(2)
      integer :: mnumfe
      integer :: element
      
C     output variables
      real(dp) :: disp(2)
      
C     local variables
      integer :: eltypenum
      real(dp) :: r, s
      logical :: badflip, failure

      call findInAllWithGuessUndef(posnundef(1),posnundef(2), ! find in mesh (undeformed)
     &                             mnumfe,element,r,s,badflip)
      if (badflip) then ! if not successful, interpolate atomistic displacements
          call atomicDispInterpolation(posnundef,failure,disp)
          if (failure) then ! if still not successful, our point must be at the boundary,
                            ! so just use the original node displacements (we'll correct boundary displacements later)
              disp = nodes%posn(node,4:5)
          end if
      else 
          eltypenum = feelements(mnumfe)%eltypenum
          disp = getTotalDispAtPoint_ptr(posnundef,
     &                                   mnumfe,eltypenum,element,r,s) ! compute displacement using fe + dd
      end if
              
      end function interpFromFEPoint
************************************************************************
      function interpFromAtomPoint(posnundef) result(disp)
     
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
     
C     input variables
      real(dp) :: posnundef(2)
      
C     output variables
      logical :: failure
      real(dp) :: disp(2)
      
C     local variables
      integer :: closestatom
      integer :: i, j
      
      closestatom = locateClosestAtomUndef(posnundef)
      failure = .true.
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
      
C     input variables
      real(dp) :: posnundef(2)
      integer :: trinum
      
C     output variables
      logical :: failure
      real(dp) :: disp(2)
      
C     local variables
      real(dp) :: posn(2,3)
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
          node = nodes%atomlist(atom)
          posn(:,i) = nodes%posn(node,1:2)
          udisp(i) = nodes%posn(node,4)
          vdisp(i) = nodes%posn(node,5)
      end do
      
C     check all three sides
      do i = 1, 3
          idx1 = i
          idx2 = mod(idx1,3) + 1
          idx3 = mod(idx2,3) + 1
          check = checkSameSide(posnundef,posn(:,idx3),
     &                          posn(:,idx1),posn(:,idx2))
          failure = (.not.check)
          if (failure) then ! we're not in triangle, so quit
              return
          end if    
      end do
      
C     if we've gotten this far, compute local position in the triangle
      call getLocalCoords(transpose(posn),movingmesh%eltypenum,
     &                    posnundef(1),posnundef(2),r,s)
                          
C     displacement is given by shape function interpolation
      N = felib(movingmesh%eltypenum)%getN_2d_ptr(r,s)
      disp(1) = dot_product(N,udisp)
      disp(2) = dot_product(N,vdisp)               
      
      end subroutine checkTriangle
************************************************************************
      function locateClosestAtomUndef(posn) result(closestatom)
      
C     input variables
      real(dp) :: posn(2)
      
C     output variables
      integer :: closestatom
      
C     local variables
      real(dp) :: x, y
      integer :: j, kx, ky
      integer :: nbinsx, nbinsy
      integer :: binx, biny, binxneigh, binyneigh
      integer :: atom, node
      real(dp) :: nodepos(2)
      real(dp) :: distmin, dist
      
      ! get the atomic bin for our position
      x = posn(1)
      y = posn(2)
      call getAtomBin(x,y,binx,biny)
      
      ! set distance to be very large value
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
              nodepos = nodes%posn(1:2,node) - nodes%posn(4:5,node)
              dist = sum((posn - nodepos)**2)
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
      
C     Notes: There are two ways to do this, each with their own advantages
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

C     input variables
      integer :: mnumfe
      integer :: isys, iplane
      
C     local variables
      integer :: i, j, k, l
      integer :: dislnum
      
      do i = 1, size(disl)
      do j = 1, size(disl(i)%splanes)
      do k = 1, size(disl(i)%splanes(j)%splane)
      do l = 1, disl(i)%splanes(j)%splane(k)%nmax
          dislnum = disl(mnumfe)%splanes(isys)%splane(iplane)%objnum(i)
          if (disl(mnumfe)%list(dislnum)%active) then
              call updateDislPosMovingMeshSub(mnumfe,isys,
     &                                        iplane,i,dislnum)
          end if
      end do
      end do
      end do
      end do
      
      end subroutine updateDislPosMovingMesh
************************************************************************
      subroutine updateDislPosMovingMeshSub(mnumfe,isys,iplane,
     &                                      iobj,dislnum)
      
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
     
      if (badflip) then ! not in mesh. Not sure how to deal with case of leaving simulation cell entirely.
                        ! If it moves into atomistic region, it should be replaced with atomistic dislocation naturally during interpolation (we hope),
                        ! so simply delete it
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

      implicit none
      
C     local variables
      integer :: i, j
      integer :: element, mnumfenew
      real(dp) :: x, y
      real(dp) :: r, s
      logical :: badflip
      
      do i = 1, size(obstacles)
      do j = 1, size(obstacles(i)%list)
          element = sources(i)%list(j)%element
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