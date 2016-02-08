      module mod_mesh_find
      
C     Purpose: Module for finding in which element in the mesh a
C     particular point lies.

C     Possible extensions: See below.

C     TODO: Need to time algorithm, see if it's fast enough.
C     TODO: Also, change brute force algorithm to generate guess in a smarter way

      use mod_types, only: dp
      use mod_fe_elements, only: feelements, nfematerials
      use mod_nodes, only: nodes
      use mod_math, only: invertMat2
      use mod_fe_el_2d, only: getN_2d
      implicit none
      
      private
      public :: findInAllWithGuess, findInAllSub, findInAllInitially,
     &   getLocalCoords, findinCPE4, getElGuessBrute,
     &   checkSameSide, checkCPE3, checkCPE4, findinCPE3, findInOneMat,
     &   findinOneMatAlt, getNodeGuessBrute, findInElement,
     &   checkElement, findInOneMatInitially

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
       
      call findInOneMat(mnumfeguess,feelements(mnumfeguess)%elname,
     &                              feelements(mnumfeguess)%nelnodes,
     &                              elguess,xp,yp,r,s,badflip)
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
      real(dp) :: r, s
      logical :: badflip
      
C     local variables
      integer :: element
      
C     get a (very good) guess
      element = getElGuessBrute(mnumfe,xp,yp)
C     search in material
      call findInOneMat(mnumfe,feelements(mnumfe)%elname,
     &           feelements(mnumfe)%nelnodes,element,xp,yp,r,s,badflip)
     
      end subroutine findInOneMatInitially
************************************************************************      
      subroutine findInOneMat(mnumfe,elname,nelnodes,elguess,
     &                        xp,yp,r,s,badflip)

C     Subroutine: findInOneMat

C     Inputs: mnumfe --- material to search in
C             elname --- name of element, character string, follows ABAQUS convention
C             nelnodes --- number of nodes per element
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
C            I'm a bit concerned that this algorithm might yield bad
C            results for search points lying on edges or vertices (nodes)...
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: nelnodes
      character(len=*) :: elname
      real(dp) :: xp, yp

C     output variables
      real(dp) :: r, s
      logical :: badflip
      
C     in/out variables
      integer :: elguess
      
C     local variables
      integer :: counter, i
      real(dp) :: posn(nelnodes,2)
      integer :: node
      logical :: proceed, check
      integer :: nedges
      integer :: eltry

      nedges = feelements(mnumfe)%nedges
      do counter = 1, countermax
          do i = 1, nelnodes
              node = feelements(mnumfe)%connect(i,elguess)
              posn(i,:) = nodes%posn(1:2,node) - nodes%posn(4:5,node) ! undeformed pos.
          end do
          r = 0.0_dp
          s = 0.0_dp
          call findInElement(posn,nelnodes,elname,xp,yp,r,s)
          proceed = .true.
          badflip = .false.
          do i = 1, nedges
              if (proceed) then
                  check = checkElement(i,elname,r,s)
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
          if (proceed) then
              if (badflip) then  ! not inside the element, since we tried to flip and failed
                  return
              else ! inside the element
                  call getLocalCoords(posn,nelnodes,elname,xp,yp,r,s)
                  return 
              end if 
          end if
      end do
      
C     if we've gotten to this point, search was unsuccessful
      badflip = .true.
      
      end subroutine findInOneMat
************************************************************************
      subroutine getLocalCoords(posn,nelnodes,elname,xp,yp,r,s)

C     Subroutine: getLocalCoords

C     Inputs: posn --- array, nelnodes by 2, of nodal positions (each row
C                      is a coordinate pair)
C             nelnodes --- number of nodes per element
C             elname --- name of element, character string, follows ABAQUS convention
C             xp, yp --- coordinates of target point


C     Outputs: r, s --- local coordinates of point in element
                      
C     Purpose: Determine local coordinates of point w.r.t. element

C     input variables
      real(dp) :: posn(nelnodes,2)
      integer :: nelnodes
      character(len=*) :: elname
      real(dp) :: xp, yp
      
C     in/out variables
      real(dp) :: r, s
      
C     local variables
      integer :: i
      real(dp) :: norm
      real(dp) :: rold, sold

      norm = huge(dp)
      do i = 1, 10
          if (norm < normconst) then
              return
          else    
              rold = r
              sold = s
              call findinElement(posn,nelnodes,elname,xp,yp,r,s)
              norm = sqrt((r - rold)**2 + (s - sold)**2)
          end if    
      end do
      
      end subroutine getLocalCoords
************************************************************************
      subroutine findInElement(posn,nelnodes,elname,xp,yp,r,s)
      
C     input variables
      real(dp) :: posn(nelnodes,2)
      integer :: nelnodes
      character(len=*) :: elname
      real(dp) :: xp, yp      

C     in/out variables
      real(dp) :: r, s
      
      if (elname == 'CPE3') then
          call findinCPE3(posn,nelnodes,xp,yp,r,s)
      else if (elname == 'CPE4') then
          call findinCPE4(posn,nelnodes,xp,yp,r,s)
      end if
      
      end subroutine findInElement
************************************************************************     
      subroutine findInCPE3(posn,nelnodes,xp,yp,r,s)

C     Subroutine: findInCPE3

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             nelnodes --- number of element nodes (in this case, 3)
C             xp, yp --- coordinates of point to search for
C             r, s --- local coordinates of point in element

C     Outputs: r, s --- updated local coordinates of point in element

C     Purpose: Determine local coordinates of point in CPE3 element, using
C              procedure described in Allievi and Bermejo, JCP, 1997
C              (if point does not lie in element, r, s will be outside bounds)
      
      implicit none
      
C     input variables
      real(dp) :: posn(nelnodes,2)
      integer :: nelnodes
      real(dp) :: xp, yp
      
C     in/out variables
      real(dp) :: r, s
      
C     local variables
      real(dp) :: diffx2, diffx3
      real(dp) :: diffy2, diffy3
      real(dp) :: diffxp, diffyp
      real(dp) :: N(nelnodes)
      real(dp) :: xold, yold
      real(dp) :: invdet
      
      diffx2 = posn(2,1) - posn(1,1)
      diffx3 = posn(3,1) - posn(1,1)
      diffy2 = posn(2,2) - posn(1,2)
      diffy3 = posn(3,2) - posn(1,2)
      N = getN_2d(r,s,'CPE3',nelnodes)
      xold = dot_product(posn(:,1),N)
      yold = dot_product(posn(:,2),N)
      diffxp = xp - xold
      diffyp = yp - yold
      invdet = 1.0_dp/(diffx2*diffy3 - diffx3*diffy2)
      r = r + invdet*(diffy3*diffxp - diffx3*diffyp)
      s = s + invdet*(diffx2*diffyp - diffy2*diffxp)
      
      end subroutine findInCPE3
************************************************************************
      subroutine findInCPE4(posn,nelnodes,xp,yp,r,s)

C     Subroutine: findInCPE4

C     Inputs: posn --- global coordinates of element nodes, 2 by nelnodes
C                      (i.e. each column is a coordinate pair)
C             nelnodes --- number of element nodes (in this case, 4)
C             xp, yp --- coordinates of point to search for
C             r, s --- local coordinates of point in element

C     Outputs: r, s --- updated local coordinates of point in element

C     Purpose: Determine local coordinates of point in CPE4 element, using
C              procedure described in Allievi and Bermejo, JCP, 1997
C              (if point does not lie in element, r, s will be outside bounds)

      implicit none
      
C     input variables
      real(dp) :: posn(nelnodes,2)
      integer :: nelnodes
      real(dp) :: xp, yp
      
C     in/out variables
      real(dp) :: r, s
      
C     local variables
      real(dp) :: diffx12, diffx34
      real(dp) :: diffy12, diffy34
      real(dp) :: diffxp, diffyp
      real(dp) :: N(nelnodes)
      real(dp) :: xold, yold
      real(dp) :: invdet
      real(dp) :: a1, a2, a3, b1, b2, b3
      real(dp) :: mat11, mat12, mat21, mat22
      
      N = getN_2d(r,s,'CPE4',nelnodes)
      xold = dot_product(posn(:,1),N)
      yold = dot_product(posn(:,2),N)
      diffxp = xp - xold
      diffyp = yp - yold
      diffx12 = posn(1,1) - posn(2,1)
      diffx34 = posn(3,1) - posn(4,1)
      diffy12 = posn(1,2) - posn(2,2)
      diffy34 = posn(3,2) - posn(4,2)
      a1 = 0.25_dp*(diffx34 - diffx12)
      a2 = 0.25_dp*(-posn(1,1) - posn(2,1) + posn(3,1)  + posn(4,1))
      a3 = 0.25_dp*(diffx34 + diffx12)
      b1 = 0.25_dp*(diffy34 - diffy12)
      b2 = 0.25_dp*(-posn(1,2) - posn(2,2) + posn(3,2)  + posn(4,2))
      b3 = 0.25_dp*(diffy34 + diffy12)
      invdet = 1.0_dp/(a1*b2-a2*b1 + (a1*b3-a3*b1)*r + (a3*b2-a2*b3)*s)
      mat11 = b2 + b3*r
      mat12 = a2 + a3*r
      mat21 = b1 + b3*s
      mat22 = a1 + a3*s
      r = r + invdet*(mat11*diffxp - mat12*diffyp)
      s = s + invdet*(mat22*diffyp - mat21*diffxp)
      
      end subroutine findInCPE4
************************************************************************
      function checkElement(edge,elname,r,s) result(check)

C     Function: checkElement

C     Inputs: edge --- number of edge of element
C             elname --- name of element, character string, follows ABAQUS convention
C             r, s --- local coordinates of point in element

C     Outputs: check --- logical value, giving whether the local coordinates
C                        lie on the correct side of the edge
                      
C     Purpose: Determine whether local coordinates from findinElement
C              lie on the correct side of edge (if checkElement(edge,r,s)
C              is true for all edges, then point lies inside element)
      
C     input variables
      integer :: edge
      character(len=*) :: elname
      real(dp) :: r, s
      
C     output variables
      logical :: check
      
      if (elname == 'CPE3') then
          check = checkCPE3(edge,r,s)
      else if (elname == 'CPE4') then
          check = checkCPE4(edge,r,s)
      end if
      
      end function checkElement
************************************************************************ 
      function checkCPE3(edge,r,s) result(check)

C     Function: checkCPE3

C     Inputs: edge --- number of edge of element
C             r, s --- local coordinates of point in element

C     Outputs: check --- logical value, giving whether the local coordinates
C                        lie on the correct side of the edge
                      
C     Purpose: Determine whether local coordinates from findinCPE3
C              lie on the correct side of edge (if checkCPE3(edge,r,s)
C              is true for all edges, then point lies inside element)
      
C     input variables
      integer :: edge
      real(dp) :: r, s
      
C     output variables
      logical :: check
      
      if (edge == 1) then
          check = (s >= 0)
      else if (edge == 2) then
          check = (r + s <= 1)
      else
          check = (r >= 0)
      end if
      
      end function checkCPE3
************************************************************************
      function checkCPE4(edge,r,s) result(check)

C     Function: checkCPE4

C     Inputs: edge --- number of edge of element
C             r, s --- local coordinates of point in element

C     Outputs: check --- logical value, giving whether the local coordinates
C                        lie on the correct side of the edge
                      
C     Purpose: Determine whether local coordinates from findinCPE4
C              lie on the correct side of edge (if checkCPE4(edge,r,s)
C              is true for all edges, then point lies inside element)
      
C     input variables
      integer :: edge
      real(dp) :: r, s
      
C     output variables
      logical :: check
      
      if (edge == 1) then
          check = (s >= -1)
      else if (edge == 2) then
          check = (r <= 1)
      else if (edge == 3) then
          check = (s <= 1)
      else
          check = (r >= -1)
      end if
      
      end function checkCPE4
************************************************************************  
      subroutine findInOneMatAlt(mnumfe,elname,nelnodes,elguess,
     &                        xp,yp,r,s,badflip)

C     Subroutine: findInOneMatAlt

C     Inputs: mnumfe --- material to search in
C             elname --- name of element, character string, follows ABAQUS convention
C             nelnodes --- number of nodes per element
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
      character(len=*) :: elname
      integer :: nelnodes
      real(dp) :: xp, yp

C     output variables
      real(dp) :: r, s
      logical :: badflip
      
C     in/out variables
      integer :: elguess
      
C     local variables
      integer :: counter, i
      real(dp) :: posn(2,nelnodes)
      real(dp) :: tposn(nelnodes,2)
      real(dp) :: pt(2)
      integer :: node
      integer :: nedges
      integer :: idx1, idx2, idx3
      logical :: proceed, check
      integer :: eltry

      pt = [xp,yp]
      nedges = feelements(mnumfe)%nedges
      do counter = 1, countermax
          do i = 1, nelnodes
              node = feelements(mnumfe)%connect(i,elguess)
C             undeformed positions
              posn(:,i) = nodes%posn(1:2,node) - nodes%posn(4:5,node)
          end do
          proceed = .true.
          badflip = .false.
          do i = 1, nedges
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
          if (proceed) then
              if (badflip) then ! not inside the element, since we tried to flip and failed
                  return
              else ! inside the element
                  r = 0.0_dp
                  s = 0.0_dp
                  tposn = transpose(posn)
                  call getLocalCoords(tposn,
     &                                nelnodes,elname,xp,yp,r,s)
                  return  
              end if 
          end if
      end do
      
C     if we've gotten to this point, search was unsuccessful
      badflip = .true. 
      
      end subroutine findInOneMatAlt
************************************************************************
      function checkSameSide(pt1,pt2,vert1,vert2) result(check)

C     Subroutine: checkSameSide

C     Inputs: pt1, pt2 --- coordinates of two points of interest (length 2)
C             vert1, vert2 --- coordinates of two points on a line (length 2)

C     Outputs: check --- true or false, depending on whether pt1 and pt2
C                        lie on the same or opposite sides of the line formed
C                        by vert1 and vert2

C     Purpose: Determine whether pt1 and pt2 lie on same side of line formed
C              by vert1 and vert2. Used as helper for findInOneMatAlt
      
C     input variables
      real(dp) :: pt1(2), pt2(2)
      real(dp) :: vert1(2), vert2(2)
      
C     output variables
      logical :: check
      
C     local variables
      real(dp) :: vec(2)
      real(dp) :: res1, res2
      
      vec = vert2 - vert1
      res1 = vec(1)*(pt1(2) - vert1(2)) - vec(2)*(pt1(1) - vert1(1))
      res2 = vec(1)*(pt2(2) - vert1(2)) - vec(2)*(pt2(1) - vert1(1)) 
C     check if res1*res2 >= 0
      if (res1 >= 0) then
          check = (res2 >= 0)
      else
          check = (res2 <= 0)
      end if    
      
      end function checkSameSide
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
      real(dp) :: xp, yp
      integer :: mnumfe
      
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

C     input variables
      real(dp) :: xp, yp
      integer :: mnumfe
      
C     output variables
      integer :: closestnode
      
C     local variables
      integer :: j
      integer :: node
      real(dp) :: distsq, distsqtry
      real(dp) :: pt(2)
      
      pt = [xp,yp]
      distsq = huge(dp)
      do j = 1, size(feelements(mnumfe)%nodelist)
          node = feelements(mnumfe)%nodelist(j)
C         could use undeformed positions,
C         but this is probably close enough
          distsqtry = norm2(nodes%posn(1:2,node)-pt)
          if (distsqtry < distsq) then
              distsq = distsqtry
              closestnode = node
          end if
      end do
      
      end function getNodeGuessBrute
************************************************************************ 
      end module