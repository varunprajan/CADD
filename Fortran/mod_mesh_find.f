      module mod_mesh_find
      
C     Purpose: Module for finding in which element in the mesh and at which local position
C     a particular point lies.

C     Possible extensions: See below.

C     TODO: Need to time algorithm, see if it's fast enough.
C     TODO: Also, change brute force algorithm to generate guess in a smarter way
C     TODO: Also, need to check if there are problems with crossing over "internal boundaries" (e.g. crack face)
C     TODO: Need to fix findInOneMatSub to deal with points lying on corners/edges
C     TODO: Check pad atom local positions after deformation.

C     Long note: Merits of finding deformed vs. undeformed positions
      
C     1) Deformed positions are expensive to calculate --- we have to compute
C     dislocation displacement fields at every single FE node. In other words, we only want to calculate
C     the nodal displacements when needed (i.e. before a dump or restart).
C     2) I think there might be an inconsistency in using local positions (e.g, for dislocations)
C     in the deformed coordinate system if the stiffness matrix is assembled
C     using positions in the undeformed coordinate system. But perhaps
C     this is fine under the small strain assumption.
C     3) Dislocations can be "lost" or present in the atomistic region 
C     if undeformed positions are used.
C     4) Passing becomes extremely tricky if
C     undeformed positions for the interface are used. The dislocation
C     can be within the (original) mesh in the undeformed coordinates,
C     but actually in empty space (having crossed the crack plane)
C     in the deformed coordinates. This makes passing the dislocation
C     from the continuum to the atomistic region basically impossible.
      
C     So, I have adopted the following solution:
C     1) Pad atom positions are found once and for all at the beginning
C     of the program (during initialization). Presumably, displacements are zero
C     at this point, so there is no difference between deformed and undeformed positions.
C     Even assuming the elements deform, the pad atom local position
C     in the element should remain constant (since the pad atom deformation
C     tracks that of the surrounding nodes), so this approach should work
C     even for restart files (for instance). (This needs to be checked.)
C     2) Dislocation positions are found at every DD timestep with respect
C     to the deformed elements. This should be fine even if these elements
C     become highly distorted.
C     3) A dislocation is passed or escaped if it leaves the deformed mesh,
C     e.g., by crossing the deformed interface. Since all of the passing
C     is done with respect to the deformed coordinates (e.g., of the atoms),
C     everything should be self-consistent, and no dislocations should be "lost"

      use mod_types, only: dp
      use mod_fe_elements, only: feelements, nfematerials
      use mod_nodes, only: nodes
      use mod_math, only: checkSameSide
      use mod_fe_el_2d, only: felib, getElTypeNum
      implicit none
      
      private
      public :: findInAllWithGuess, findInAllSub,
     &   getLocalCoords, getElGuessBrute, findInOneMat,
     &   findinOneMatSub, getNodeGuessBrute, findInOneMatInitially,
     &   findInAllWithGuessDef, findInAllWithGuessUndef,
     &   findInAllInitiallyDef, findInAllInitiallyUndef,
     &   findInOneMatInitiallyDef

C     module variables (private)
C     HARD-CODED CONSTANTS
      integer, parameter :: COUNTERMAX = 1000
      integer, parameter :: COUNTERMAX2 = 20
      integer, parameter :: MAXEL = 10
      real(dp), parameter :: NORMCONST = 1.0e-5_dp
      
      contains
************************************************************************
      subroutine findInAllWithGuessUndef(xp,yp,
     &                                  mnumfeguess,elguess,r,s,badflip)

C     Inputs: xp, yp --- coordinates of point to search

C     In/Out: mnumfeguess --- material to search in/material point was found in
C             elguess --- guess for starting element/element point was found in

C     Outputs: r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was unsuccessful

C     Purpose: findInAllWithGuess, using *undeformed* coordinates
      
      implicit none
      
C     input variables
      real(dp) :: xp, yp
      
C     in/out variables
      integer :: mnumfeguess, elguess
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      
      call findInAllWithGuess(.true.,xp,yp,
     &                        mnumfeguess,elguess,r,s,badflip)
      
      end subroutine findInAllWithGuessUndef
************************************************************************
      subroutine findInAllWithGuessDef(xp,yp,
     &                                  mnumfeguess,elguess,r,s,badflip)

C     Inputs: xp, yp --- coordinates of point to search

C     In/Out: mnumfeguess --- material to search in/material point was found in
C             elguess --- guess for starting element/element point was found in

C     Outputs: r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was unsuccessful

C     Purpose: findInAllWithGuess, using *deformed* coordinates
      
      implicit none
      
C     input variables
      real(dp) :: xp, yp
      
C     in/out variables
      integer :: mnumfeguess, elguess
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      
      call findInAllWithGuess(.false.,xp,yp,
     &                        mnumfeguess,elguess,r,s,badflip)
      
      end subroutine findInAllWithGuessDef
************************************************************************       
      subroutine findInAllWithGuess(undeformed,xp,yp,
     &                              mnumfeguess,elguess,r,s,badflip)

C     Inputs: undeformed --- flag indicating whether undeformed coordinates should be used
C             xp, yp --- coordinates of point to search

C     In/Out: mnumfeguess --- material to search in/material point was found in
C             elguess --- guess for starting element/element point was found in

C     Outputs: r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was unsuccessful

C     Purpose: Search over FE mesh for all materials to
C              find which material and element a point belongs to. Also
C              return local coordinates of point in element. Resorts to brute
C              force search if point has left material mnumfeguess.
      
      implicit none
      
C     input variables
      logical :: undeformed
      real(dp) :: xp, yp
      
C     in/out variables
      integer :: mnumfeguess, elguess
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      
      call findInOneMat(mnumfeguess,undeformed,
     &                  elguess,xp,yp,r,s,badflip)
C     if search was unsuccessful
      if (badflip) then
          call findInAllSub([mnumfeguess],undeformed,
     &                      xp,yp,mnumfeguess,elguess,r,s,badflip)
      end if
      
      end subroutine findInAllWithGuess
************************************************************************
      subroutine findInAllInitiallyUndef(
     &                                 xp,yp,mnumfe,element,r,s,badflip)

C     Inputs: xp, yp --- coordinates of point to search

C     Outputs: mnumfe --- material that point was found in
C              element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was unsuccessful

C     Purpose: Used for initialization, where no information is available about where
C              a point might be. Search over FE mesh using brute force search. Return
C              error if point is not found. Uses *undeformed* coordinates.
      
      implicit none
      
C     input variables
      real(dp) :: xp, yp
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      integer :: mnumfe, element
      
      call findInAllSub([0],.true.,xp,yp,mnumfe,element,r,s,badflip)   
      
      end subroutine findInAllInitiallyUndef
************************************************************************
      subroutine findInAllInitiallyDef(xp,yp,mnumfe,element,r,s,badflip)

C     Inputs: xp, yp --- coordinates of point to search

C     Outputs: mnumfe --- material that point was found in
C              element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was unsuccessful

C     Purpose: Used for initialization, where no information is available about where
C              a point might be. Search over FE mesh using brute force search. Return
C              error if point is not found. Uses *deformed* coordinates
      
      implicit none
      
C     input variables
      real(dp) :: xp, yp
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      integer :: mnumfe, element
      
      call findInAllSub([0],.false.,xp,yp,mnumfe,element,r,s,badflip)   
      
      end subroutine findInAllInitiallyDef
************************************************************************
      subroutine findInAllSub(alreadysearched,undeformed,
     &                        xp,yp,mnumfe,element,r,s,badflip)

C     Inputs: alreadysearched --- list of materials *not* to search in
C             undeformed --- flag indicating whether undeformed coordinates should be used
C             xp, yp --- coordinates of point to search

C     Outputs: mnumfe --- material that point was found in
C              element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating that search was unsuccessful

C     Purpose: Helper routine for findInAllInitially and findInAllWithGuess.
C     Loops over materials. For each, generates guess using brute force,
C     and then tries to find point using findInOneMat
      
      implicit none
      
C     input variables
      integer :: alreadysearched(:)
      logical :: undeformed
      real(dp) :: xp, yp
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      integer :: mnumfe, element
      
C     local variables
      logical :: search
      integer :: i, j
      
C     loop over all materials, except the ones that have been searched already
      do i = 1, nfematerials
          search = .true.
          do j = 1, size(alreadysearched)
              if (alreadysearched(j) == i) then
                  search = .false.
              end if
          end do    
          if (search) then
              call findInOneMatInitially(i,undeformed,
     &                                   xp,yp,element,r,s,badflip)
              if (.not.badflip) then
                  mnumfe = i
                  return
              end if    
          end if
      end do
      
      end subroutine findInAllSub
************************************************************************
      subroutine findInOneMatInitiallyDef(mnumfe,
     &                                 xp,yp,element,r,s,badflip)

C     Inputs: mnumfe --- material to search in
C             xp, yp --- coordinates of point to search

C     Outputs: element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating that search was unsuccessful

C     Purpose: Helper routine for findInAllInitially and findInAllWithGuess.
C     Generates guess for element using brute force,
C     and then tries to find point using findInOneMat.
      
      implicit none
      
C     input variables
      integer :: mnumfe
      real(dp) :: xp, yp
      
C     output variables
      integer :: element
      real(dp) :: r, s
      logical :: badflip
      
      call findInOneMatInitially(mnumfe,.false.,
     &                           xp,yp,element,r,s,badflip)
     
      end subroutine findInOneMatInitiallyDef
************************************************************************
      subroutine findInOneMatInitially(mnumfe,undeformed,
     &                                 xp,yp,element,r,s,badflip)

C     Inputs: mnumfe --- material to search in
C             undeformed --- flag indicating whether undeformed coordinates should be used
C             xp, yp --- coordinates of point to search

C     Outputs: element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating that search was unsuccessful

C     Purpose: Helper routine for findInAllInitially and findInAllWithGuess.
C     Generates guess for element using brute force,
C     and then tries to find point using findInOneMat.
      
      implicit none
      
C     input variables
      integer :: mnumfe
      logical :: undeformed
      real(dp) :: xp, yp
      
C     output variables
      integer :: element
      real(dp) :: r, s
      logical :: badflip
      
C     local variables
      integer :: i
      integer, allocatable :: elements(:)
      
C     get very good guesses for elements
      elements = getElGuessBrute(mnumfe,undeformed,xp,yp)
C     search in material
      do i = 1, size(elements)
          element = elements(i)
          call findInOneMat(mnumfe,undeformed,element,xp,yp,r,s,badflip)
          if (.not.badflip) then
              return
          end if
      end do    
     
      end subroutine findInOneMatInitially
************************************************************************
      subroutine findInOneMat(mnumfe,undeformed,
     &                        elguess,xp,yp,r,s,badflip)

C     Inputs: mnumfe --- material to search in
C             undeformed --- flag indicating whether undeformed coordinates should be used         
C             xp, yp --- coordinates of point to search

C     In/out: elguess --- (in) guess for starting element; (out) element that point was found in 

C     Outputs: r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating that search was unsuccessful

C     Purpose: Tries to find point in mesh using elguess as the starting element.
C     Calls jump-and-walk algorithm.

      implicit none
      
C     input variables
      integer:: mnumfe
      logical :: undeformed    
      integer :: elguess
      real(dp) :: xp, yp
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      
C     local variables
      integer:: eltypenum
      
      eltypenum = feelements(mnumfe)%eltypenum      
      call findInOneMatSub(mnumfe,eltypenum,undeformed,
     &                     elguess,xp,yp,r,s,badflip)
      
      end subroutine findInOneMat
************************************************************************
      subroutine findInOneMatSub(mnumfe,eltypenum,undeformed,
     &                              elguess,xp,yp,r,s,badflip)

C     Inputs: mnumfe --- material to search in
C             eltypenum --- number of element type in fe element library
C             undeformed --- flag indicating whether undeformed coordinates should be used
C             xp, yp --- coordinates of point to search

C     In/out: elguess --- (in) guess for starting element; (out) element that point was found in

C     Outputs: r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was successful

C     Purpose: Search over FE mesh for one material (mnumfe) to
C              determine if point is located in the mesh of that material.
C              If so, return element in which it is located and the
C              local coordinates of point in the element. If not, 
C              set badflip = .true. (flag indicating failure)

C     Algorithm: True jump and walk algorithm:
C       1) start with an element,
C       2) pick an edge
C       3) determine if point of interest lies on same side of edge as another node in element
C       4) if not, flip to adjacent element (from neighbor list)
C       5) if so, loop over other edges
C       6) if this loop is successful, we've found the correct element,
C        and can iterate to find local coordinates using Newton-Raphson routine
      
C     Notes: Supplying a good guess is very important, especially
C            for non-convex domains (where getting stuck is a possibility)
C            Might be a good idea to find the object (dislocation, pad atom)
C            initially using a brute force search (mod_mesh_find_brute).
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: eltypenum
      logical :: undeformed
      real(dp) :: xp, yp

C     output variables
      real(dp) :: r, s
      logical :: badflip
      
C     in/out variables
      integer :: elguess
      
C     local variables
      integer :: counter, i
      real(dp) :: posn(2,felib(eltypenum)%nelnodes)
      real(dp) :: pt(2)
      integer :: node
      integer :: nelnodes
      integer :: idx1, idx2, idx3
      logical :: proceed, check
      integer :: eltry

      pt = [xp,yp]
      nelnodes = felib(eltypenum)%nelnodes
      do counter = 1, COUNTERMAX
          do i = 1, nelnodes
              node = feelements(mnumfe)%connect(i,elguess)
              posn(:,i) = nodes%posn(1:2,node)
              if (undeformed) then
                  posn(:,i) = posn(:,i) - nodes%posn(4:5,node)
              end if    
          end do
          proceed = .true.
          badflip = .false.
          do i = 1, felib(eltypenum)%nedges
              if (proceed) then
                  idx1 = i
                  idx2 = mod(idx1,nelnodes) + 1
                  idx3 = mod(idx2,nelnodes) + 1
                  check = checkSameSide(pt,posn(:,idx3),
     &                                  posn(:,idx1),posn(:,idx2))
                  if (.not.check) then ! try to flip to adjacent element
                      eltry = feelements(mnumfe)%neighbors(i,elguess)
                      proceed = (eltry == 0)  ! if element doesn't exist, proceed
                      if (proceed) then
                          badflip = .true.
                      else
                          elguess = eltry
                      end if
                  end if
              end if
          end do
          if (proceed) then  ! we're either inside the element or a bad flip happened
              if (.not.badflip) then ! inside the element
                  call getLocalCoords(transpose(posn),
     &                                eltypenum,xp,yp,r,s)
              end if
              return 
          end if
      end do
      
C     if we've gotten to this point, search was unsuccessful
      write(*,*) 'WARNING: Search in findInOneMatSub
     &            took way too long'
      badflip = .true. 
      
      end subroutine findInOneMatSub
************************************************************************
      subroutine getLocalCoords(posn,eltypenum,xp,yp,r,s)

C     Inputs: posn --- array, nelnodes by 2, of nodal positions (each row
C                      is a coordinate pair)
C             eltypenum --- number code of fe element in felib
C             xp, yp --- coordinates of target point

C     Outputs: r, s --- local coordinates of point in element
                      
C     Purpose: Determine local coordinates of point w.r.t. element

      implicit none

C     input variables
      real(dp) :: posn(:,:)
      integer :: eltypenum
      real(dp) :: xp, yp
      
C     output variables
      real(dp) :: r, s
      
C     local variables
      integer :: i
      real(dp) :: norm
      real(dp) :: rold, sold

      norm = huge(0.0_dp)
      r = 0.0_dp
      s = 0.0_dp
      do i = 1, COUNTERMAX2
          if (norm < NORMCONST) then
              return
          else    
              rold = r
              sold = s             
              call felib(eltypenum)%findinElement_ptr(posn,xp,yp,r,s)
              norm = sqrt((r - rold)**2 + (s - sold)**2)
          end if    
      end do
      write(*,*) 'Getting local coordinates took way too long'
      
      end subroutine getLocalCoords
************************************************************************
      function getElGuessBrute(mnumfe,undeformed,xp,yp) result(elguess)

C     Inputs: xp, yp --- coordinates of point to search
C             undeformed --- flag indicating whether undeformed coordinates should be used
C             mnumfe --- continuum material to search in

C     Outputs: elguess --- list of elements that are reasonably close to point

C     Purpose: Find many elements that are reasonably close the point of interest
C              within a certain material, using brute force search
C              (i.e. to supply initial guess for more sophisticated
C              point finding schemes)
      
      implicit none
      
C     input variables
      integer :: mnumfe
      logical :: undeformed    
      real(dp) :: xp, yp
      
C     output variables
      integer, allocatable :: elguess(:)
      
C     local variables
      integer :: j, k
      integer :: closestnode
      integer :: counter
      integer :: elguesstemp(MAXEL)

      closestnode = getNodeGuessBrute(mnumfe,undeformed,xp,yp)
      counter = 0
      do j = 1, size(feelements(mnumfe)%connect,2)
          do k = 1, size(feelements(mnumfe)%connect,1)
              if (feelements(mnumfe)%connect(k,j) == closestnode) then
                  counter = counter + 1
                  elguesstemp(counter) = j
              end if    
          end do
      end do
      
      elguess = elguesstemp(1:counter)
      
      end function getElGuessBrute
************************************************************************
      function getNodeGuessBrute(mnumfe,undeformed,xp,yp)
     &                                              result(closestnode)

C     Inputs: mnumfe --- continuum material to search in
C             undeformed --- flag indicating whether undeformed coordinates should be used
C             xp, yp --- coordinates of point to search 

C     Outputs: closestnode --- number of node closest to point of interest

C     Purpose: Find closest fe node within material to point of interest
C     (in the deformed configuration)

      implicit none

C     input variables
      integer :: mnumfe
      logical :: undeformed  
      real(dp) :: xp, yp  
      
C     output variables
      integer :: closestnode
      
C     local variables
      integer :: j
      integer :: node
      real(dp) :: distsq, distsqtry
      real(dp) :: pt(2), posn(2)
      
      pt = [xp,yp]
      distsq = huge(0.0_dp)
      do j = 1, size(feelements(mnumfe)%nodelist)
          node = feelements(mnumfe)%nodelist(j)
          posn = nodes%posn(1:2,node) ! deformed positions
          if (undeformed) then
              posn = posn - nodes%posn(4:5,node)
          end if    
          distsqtry = sum((posn-pt)**2)
          if (distsqtry < distsq) then
              distsq = distsqtry
              closestnode = node
          end if
      end do
      
      end function getNodeGuessBrute
************************************************************************ 
      end module