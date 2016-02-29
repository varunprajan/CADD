      module mod_mesh_find
      
C     Purpose: Module for finding in which element in the mesh a
C     particular point lies.

C     Possible extensions: See below.

C     TODO: Need to time algorithm, see if it's fast enough.
C     TODO: Also, change brute force algorithm to generate guess in a smarter way
C     TODO: Also, need to check if there are problems with crossing over "internal boundaries" (e.g. crack face)
C     TODO: Need to fix findInOneMatSub to deal with points lying on corners/edges

      use mod_types, only: dp
      use mod_fe_elements, only: feelements, nfematerials
      use mod_nodes, only: nodes
      use mod_math, only: checkSameSide
      use mod_fe_el_2d, only: felib, getElTypeNum
      implicit none
      
      private
      public :: findInAllWithGuess, findInAllSub, findInAllInitially,
     &   getLocalCoords, getElGuessBrute, findInOneMat, findInOneMatSub,
     &   findinOneMatSubAlt, getNodeGuessBrute, findInOneMatInitially

C     module variables (private)
C     HARD-CODED CONSTANTS
      integer, parameter :: countermax = 1000
      real(dp), parameter :: normconst = 1.0e-5_dp
      
      contains
************************************************************************
      subroutine findInAllWithGuess(xp,yp,mnumfeguess,elguess,
     &                               r,s,badflip)
      
C     Subroutine: findInAllWithGuess

C     Inputs: xp, yp --- coordinates of point to search

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
      real(dp) :: xp, yp
      
C     in/out variables
      integer :: mnumfeguess, elguess
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      
      call findInOneMat(mnumfeguess,elguess,xp,yp,r,s,badflip)
C     if search was unsuccessful
      if (badflip) then
          call findInAllSub([mnumfeguess],xp,yp,mnumfeguess,elguess,
     &                                    r,s,badflip)
      end if
      
      end subroutine findInAllWithGuess
************************************************************************
      subroutine findInAllInitially(xp,yp,mnumfe,element,r,s,badflip)

C     Subroutine: findInAllInitially

C     Inputs: xp, yp --- coordinates of point to search

C     Outputs: mnumfe --- material that point was found in
C              element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was unsuccessful

C     Purpose: Used for initialization, where no information is available about where
C              a point might be. Search over FE mesh using brute force search. Return
C              error if point is not found.
      
C     input variables
      real(dp) :: xp, yp
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      integer :: mnumfe, element
      
      call findInAllSub([0],xp,yp,mnumfe,element,r,s,badflip)   
      
      end subroutine findInAllInitially
************************************************************************
      subroutine findInAllSub(alreadysearched,xp,yp,mnumfe,element,
     &                                        r,s,badflip)

C     Subroutine: findInAllSub

C     Inputs: alreadysearched --- list of materials *not* to search in
C             xp, yp --- coordinates of point to search

C     Outputs: mnumfe --- material that point was found in
C              element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating that search was unsuccessful

C     Purpose: Helper routine for findInAllInitially and findInAllWithGuess.
C     Loops over materials. For each, generates guess using brute force,
C     and then tries to find point using findInOneMat
      
C     input variables
      integer :: alreadysearched(:)
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
              call findInOneMatInitially(i,xp,yp,element,r,s,badflip)
              if (.not.badflip) then
                  mnumfe = i
                  return
              end if    
          end if
      end do
      
      end subroutine findInAllSub
************************************************************************
      subroutine findInOneMatInitially(mnumfe,xp,yp,element,r,s,badflip)

C     Subroutine: findInOneMatInitially

C     Inputs: mnumfe --- material to search in
C             xp, yp --- coordinates of point to search

C     Outputs: element --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating that search was unsuccessful

C     Purpose: Helper routine for findInAllInitially and findInAllWithGuess.
C     Loops over materials. For each, generates guess using brute force,
C     and then tries to find point using findInOneMat
      
C     input variables
      integer :: mnumfe
      real(dp) :: xp, yp
      
C     output variables
      integer :: element
      real(dp) :: r, s
      logical :: badflip
      
C     get a (very good) guess
      element = getElGuessBrute(mnumfe,xp,yp)
C     search in material
      call findInOneMat(mnumfe,element,xp,yp,r,s,badflip)
     
      end subroutine findInOneMatInitially
************************************************************************
      subroutine findInOneMat(mnumfe,elguess,xp,yp,r,s,badflip)
      
C     input variables
      integer:: mnumfe
      integer :: elguess
      real(dp) :: xp, yp
      
C     output variables
      real(dp) :: r, s
      logical :: badflip
      
C     local variables
      integer:: eltypenum
      
      eltypenum = feelements(mnumfe)%eltypenum      
      call findInOneMatSubAlt(mnumfe,eltypenum,elguess,
     &                        xp,yp,r,s,badflip)
      
      end subroutine findInOneMat
************************************************************************      
      subroutine findInOneMatSub(mnumfe,eltypenum,
     &                           elguess,xp,yp,r,s,badflip)

C     Subroutine: findInOneMatSub

C     Inputs: mnumfe --- material to search in
C             eltypenum --- number of element type in fe element library
C             elguess --- guess for starting element
C             xp, yp --- coordinates of point to search

C     Outputs: elguess --- element that point was found in
C              r, s --- local (element-level) coordinates of point
C              badflip --- flag indicating whether search was unsuccessful

C     Purpose: Search over FE mesh for one material (mnumfe) to
C              determine if point is located in the mesh of that material.
C              If so, return element in which it is located and the
C              local coordinates of point in the element. If not, 
C              set badflip = .true. (flag indicating failure)

C     Algorithm: Based on the algorithm described in Allievi and Bermejo, JCP, 1997
C     Similar to a jump and walk algorithm:
C       1) start with an element,
C       2) determine approximate local coordinates of point w.r.t. element
C       3) if points lie within bounds, then iterate to find correct local coordinates
C       4) if not, flip to adjacent element (from neighbor list), depending on local coordinates
      
C     Notes: Supplying a good guess is very important, especially
C            for non-convex domains (where getting stuck is a possibility)

C     IMPORTANT NOTE: This algorithm appears to be flawed for points lying on element
C                     vertices...not sure how to fix. See bad result in findInOneMatSub_test.txt
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: eltypenum
      real(dp) :: xp, yp

C     output variables
      real(dp) :: r, s
      logical :: badflip
      
C     in/out variables
      integer :: elguess
      
C     local variables
      integer :: counter, i
      real(dp) :: posn(felib(eltypenum)%nelnodes,2)
      integer :: node
      logical :: proceed, check
      integer :: eltry

      do counter = 1, countermax
          do i = 1, felib(eltypenum)%nelnodes
              node = feelements(mnumfe)%connect(i,elguess)
              posn(i,:) = nodes%posn(1:2,node) - nodes%posn(4:5,node) ! undeformed positions
          end do
          r = 0.0_dp
          s = 0.0_dp
          call felib(eltypenum)%findinElement_ptr(posn,xp,yp,r,s)
          proceed = .true.
          badflip = .false.
          do i = 1, felib(eltypenum)%nedges
              if (proceed) then
                  check = felib(eltypenum)%checkElement_ptr(i,r,s)
                  if (.not.check) then ! try to flip to adjacent element
                      eltry = feelements(mnumfe)%neighbors(i,elguess)
                      proceed = (eltry == 0) ! if element doesn't exist, proceed
                      if (proceed) then
                          badflip = .true.
                      else
                          elguess = eltry
                      end if
                  end if
              end if
          end do
          if (proceed) then ! we're either inside the element or a bad flip happened
              if (.not.badflip) then ! inside the element
                  call getLocalCoords(posn,eltypenum,xp,yp,r,s)
              end if
              return
          end if
      end do
      
C     if we've gotten to this point, search was unsuccessful
      write(*,*) 'Search in findInOneMatSub took way too long'
      badflip = .true.
      
      end subroutine findInOneMatSub
************************************************************************
      subroutine findInOneMatSubAlt(mnumfe,eltypenum,
     &                              elguess,xp,yp,r,s,badflip)

C     Subroutine: findInOneMatSubAlt

C     Inputs: mnumfe --- material to search in
C             eltypenum --- number of element type in fe element library
C             elguess --- guess for starting element
C             xp, yp --- coordinates of point to search

C     Outputs: elguess --- element that point was found in
C              r, s --- local (element-level) coordinates of point
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
      do counter = 1, countermax
          do i = 1, nelnodes
              node = feelements(mnumfe)%connect(i,elguess)
              posn(:,i) = nodes%posn(1:2,node) - nodes%posn(4:5,node) ! undeformed positions
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
                  r = 0.0_dp
                  s = 0.0_dp
                  call getLocalCoords(transpose(posn),
     &                                eltypenum,xp,yp,r,s)
              end if
              return 
          end if
      end do
      
C     if we've gotten to this point, search was unsuccessful
      write(*,*) 'Search in findInOneMatSubAlt took way too long'
      badflip = .true. 
      
      end subroutine findInOneMatSubAlt
************************************************************************
      subroutine getLocalCoords(posn,eltypenum,xp,yp,r,s)

C     Subroutine: getLocalCoords

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
      
C     in/out variables
      real(dp) :: r, s
      
C     local variables
      integer :: i
      real(dp) :: norm
      real(dp) :: rold, sold

      norm = huge(0.0_dp)
      do i = 1, 10
          if (norm < normconst) then
              return
          else    
              rold = r
              sold = s
              call felib(eltypenum)%findinElement_ptr(posn,xp,yp,r,s)
              norm = sqrt((r - rold)**2 + (s - sold)**2)
          end if    
      end do
      
      end subroutine getLocalCoords
************************************************************************
      function getElGuessBrute(mnumfe,xp,yp) result(elguess)
      
C     Subroutine: getElGuessBrute

C     Inputs: xp, yp --- coordinates of point to search
C             mnumfe --- continuum material to search in

C     Outputs: elguess --- element that is reasonably close to point

C     Purpose: Find an element reasonably close the point of interest
C              within a certain material, using brute force search
C              (i.e. to supply initial guess for more sophisticated
C              point finding schemes)
      
      implicit none
      
C     input variables
      integer :: mnumfe
      real(dp) :: xp, yp
      
C     output variables
      integer :: elguess
      
C     local variables
      integer :: j, k
      integer :: closestnode

      closestnode = getNodeGuessBrute(mnumfe,xp,yp)
      do j = 1, size(feelements(mnumfe)%connect,2)
          do k = 1, size(feelements(mnumfe)%connect,1)
              if (feelements(mnumfe)%connect(k,j) == closestnode) then
                  elguess = j
                  return
              end if    
          end do
      end do
      
      end function getElGuessBrute
************************************************************************
      function getNodeGuessBrute(mnumfe,xp,yp) result(closestnode)

C     Subroutine: getNodeGuessBrute

C     Inputs: xp, yp --- coordinates of point to search
C             mnumfe --- continuum material to search in

C     Outputs: closestnode --- number of node closest to point of interest

C     Purpose: Find closest fe node within material to point of interest
C     (in the deformed configuration)

      implicit none

C     input variables
      real(dp) :: xp, yp
      integer :: mnumfe
      
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
          posn = nodes%posn(1:2,node) - nodes%posn(4:5,node)
          distsqtry = sum((posn-pt)**2)
          if (distsqtry < distsq) then
              distsq = distsqtry
              closestnode = node
          end if
      end do
      
      end function getNodeGuessBrute
************************************************************************ 
      end module