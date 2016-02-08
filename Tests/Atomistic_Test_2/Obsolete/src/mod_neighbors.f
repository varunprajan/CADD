      module mod_neighbors

C     Purpose: Has routines for generating neighbor lists, and checking
C     displacements to see if neighbor lists should be regenerated.
C     Also contains information related to neighbors (cutoff distance,
C     etc.)

C     Possible extensions: ?

      use mod_types, only: dp
      use mod_utils, only: prettyPrintIntMat, prettyPrintRealMat
      use mod_math, only: tolconst, isMultiple
      use mod_nodes, only: nodes
      use mod_materials, only: materials, nmaterials
      use mod_potentials, only: potentials, npotentials
      use mod_groups, only: groups, genGroupMaskAtoms
      use mod_misc, only: misc
      implicit none
      
      private
      public :: neighbors, updateNeighbors, initNeighborData,
     &          writeNeighborData, genNeighList, checkDisp,
     &          updatePosSinceLastCheck, readNeighborData,
     &          processNeighborData, getXYAtomBounds, getAtomsInBoxGroup
      
      type neighbordata
C     read-in
      logical :: checkdisp
      integer :: delay
      integer :: every
      integer :: images
      real(dp) :: Lz
      real(dp) :: skin
C     processed
      integer, allocatable :: bincount(:,:)
      integer, allocatable :: binarray(:,:,:)
      integer, allocatable :: neighlist(:,:)
      integer, allocatable :: neighcount(:)
      integer :: nmaxneigh
      integer :: incrsincelastupdate
      real(dp), allocatable :: possincelastcheck(:,:)
      real(dp) :: rhomax
      real(dp) :: rneigh
      end type
      
      type atomboxdata
      real(dp) :: xmin
      real(dp) :: xmax
      real(dp) :: ymin
      real(dp) :: ymax
      end type
      
C     module variables (global)
      type(neighbordata) :: neighbors
      type(atomboxdata) :: atombox
C     HARD-CODED CONSTANTS
C     factors of safety for computing max. neighbors, etc.
C     (see initNeighbors and genNeighList)
      real(dp) :: nmaxneighfac = 2.0_dp
      real(dp) :: nmaxbinfac = 2.5_dp
      
      contains
************************************************************************      
      subroutine initNeighborData(neighborfile)

C     Subroutine: initNeighborData

C     Inputs: neighborfile - name of file containing neighbors information
C     (e.g. example_fortran_neighbors)

C     Outputs: None

C     Purpose: Read, initialize data in "neighbors" structure, which
C     contains information about building neighbor lists
      
      implicit none
      
C     input variables
      character(len=*) :: neighborfile
            
      call readNeighborData(neighborfile)
      call processNeighborData()
      
      end subroutine initNeighborData
************************************************************************      
      subroutine readNeighborData(neighborfile)
      
C     Subroutine: readNeighborData

C     Inputs: neighborfile - name of file containing neighbors information
C     (e.g. example_fortran_neighbors)

C     Outputs: None

C     Purpose: Read data into "neighbors" from file
      
      implicit none
      
C     input variables
      character(len=*) :: neighborfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=neighborfile)
      
      read(iunit,*) temp
C     use explicit integer -> logical conversion
      neighbors%checkdisp = (temp /= 0)
      read(iunit,*) neighbors%delay
      read(iunit,*) neighbors%every
      read(iunit,*) neighbors%images
      read(iunit,*) neighbors%Lz
      read(iunit,*) neighbors%skin
      
      close(iunit)
      
      end subroutine readNeighborData
************************************************************************      
      subroutine processNeighborData()
      
C     Subroutine: processNeighborData

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize various quantities (neighboring radius, etc.)
C     used for building neighbor lists, and allocates memory for arrays.
C     Also calls genNeighList (which generates neighbor lists).
C     To be used only once, at the beginning of a simulation.

C     Notes: Requires materials, potentials, nodes structures to have been
C     read/initialized already
      
C     local variables
      integer :: i
      
C     get rhomax (used to calculate nmaxbin for building bins)
      neighbors%rhomax = materials(1)%rho
      do i = 2,nmaterials
          if (materials(i)%rho > neighbors%rhomax) then
              neighbors%rhomax = materials(i)%rho
          end if
      end do

C     get rneigh (for calculating neighbor lists)
C     first calculate largest force cutoff, then add skin distance
      neighbors%rneigh = potentials(1)%forcecutoff
      do i = 2,npotentials
          if (potentials(i)%forcecutoff > neighbors%rneigh) then
              neighbors%rneigh = potentials(i)%forcecutoff
          end if
      end do
      neighbors%rneigh = neighbors%rneigh + neighbors%skin

C     calculate nmaxneigh using area of circle, density,
C     factor of safety of nmaxneighfac, then allocate neighlist
      neighbors%nmaxneigh = ceiling(nmaxneighfac*neighbors%rhomax*
     &               3.14_dp*(neighbors%rneigh*neighbors%rneigh))
      allocate(neighbors%neighlist(neighbors%nmaxneigh,nodes%natoms))
      
C     allocate neighcount
      allocate(neighbors%neighcount(nodes%natoms))
      
C     initialize neighbor list
      call genNeighList()
      
C     allocate possincelastcheck/incrsincelastupdate, initialize
      allocate(neighbors%possincelastcheck(2,nodes%natoms))
      call updatePosSinceLastCheck()
      neighbors%incrsincelastupdate = misc%incrementcurr
      
      end subroutine processNeighborData
************************************************************************      
      subroutine writeNeighborData(neighborfile)
      
C     Subroutine: writeNeighborData

C     Inputs: neighborfile - name of file containing neighbors information
C     (e.g. example_fortran_neighbors)

C     Outputs: None

C     Purpose: Write data from "neighbors" from file (essentially, inverse
C     of readNeighborData). Useful for writing "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: neighborfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=neighborfile)
      
C     use explicit logical -> integer conversion
      if (neighbors%checkdisp) then
          temp = 1
      else
          temp = 0
      end if    
      write(iunit,*) temp
      write(iunit,*) neighbors%delay
      write(iunit,*) neighbors%every
      write(iunit,*) neighbors%images
      write(iunit,*) neighbors%Lz
      write(iunit,*) neighbors%skin
      
      close(iunit)
      
      end subroutine writeNeighborData
************************************************************************
      subroutine updateNeighbors()

C     Subroutine: updateNeighbors

C     Outputs: None

C     Purpose: Recreate neighbor list. If check, then first check
C     if displacement of atoms is too large (see checkDisp).
C     Also update possincelastcheck, incrsincelastupdate
C     if neighbor list is recreated.
      
      implicit none
      
C     local variables
      logical :: update
      integer :: delayincr
      
      delayincr = misc%incrementcurr - neighbors%incrsincelastupdate
C     check delay and every conditions
      update = ((delayincr >= neighbors%delay).and.
     &          isMultiple(delayincr,neighbors%every,.true.))  
      
C     if check is on, need to check displacements
      if (update) then
          if (neighbors%checkdisp) then
              update = checkDisp()
          end if
      end if    
      
C     if we need to update, then regenerate neighbor lists,
C     and update displacements since last check array
      if (update) then
          write(*,*) 'Updated neighbors'
          call genNeighList()
          call updatePosSinceLastCheck()
          neighbors%incrsincelastupdate = misc%incrementcurr
      end if
      
      end subroutine updateNeighbors
************************************************************************
      subroutine updatePosSinceLastCheck()
      
C     Subroutine: updatePosSinceLastCheck

C     Inputs: None

C     Outputs: None

C     Purpose: Update nodes%possincelastcheck, which contains nodal
C     positions used in dispcheck. This subroutine should be invoked
C     only when neighbor lists are regenerated (using genNeighList).
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1,nodes%natoms 
          neighbors%possincelastcheck(:,i) =
     &                            nodes%posn(1:2,nodes%atomlist(i))
      end do
      
      end subroutine updatePosSinceLastCheck
************************************************************************
      function checkDisp() result(dispcheck)
      
C     Function: checkDisp

C     Inputs: None

C     Outputs: dispcheck - logical, indicates whether neighbor lists
C     should be recreated

C     Purpose: Check whether neighbor lists should be recreated.
      
C     Algorithm: Get two largest atomic displacements, see if their sum
C     exceeds skin distance. Relies on atom positions at the last check
C     being stored in the array nodes%possincelastcheck. Includes pad
C     atoms, since they are important for computing forces.
      
      implicit none
      
C     output variables
      logical :: dispcheck
      
C     local variables
      integer :: i, atom
      real(dp) :: dispmax1, dispmax2, disp
      real(dp) :: posold(2), posnew(2)
      
C     goal: get two largest displacements since last check
      dispmax1 = 0.0_dp
      dispmax2 = 0.0_dp
      do i = 1,nodes%natoms
          atom = nodes%atomlist(i)
          posold = neighbors%possincelastcheck(:,i)
          posnew = nodes%posn(1:2,atom)
          disp = norm2(posold-posnew)
          if (disp > dispmax2) then
              if (disp > dispmax1) then
                  dispmax2 = dispmax1
                  dispmax1 = disp
              else
                  dispmax2 = disp
              end if
          end if
      end do
      dispcheck = (dispmax1 + dispmax2 > neighbors%skin)
      
      end function checkDisp
************************************************************************
      subroutine genNeighList()
      
C     Subroutine: genNeighList

C     Inputs: None

C     Outputs: None

C     Purpose: Generate neighbor list, given atom positions
      
C     Algorithm: Goal is to find all atoms within radius rneigh of given
C     atom. Use binning algorithm. Bin atoms in square bins with side
C     length rneigh. After binning, go through binlist to get neighbors.
C     Exploit symmetry to avoid duplicate work (i.e. use fact that if i
C     is a neighbor of j, j is a neighbor of i)
      
C     Notes: This algorithm is efficient only when atoms are roughly
C     contiguous (otherwise, bins are too large).
C     Also, note that the result (neighbors%neighlist) contains entries
C     that are indices of nodes%atomlist. To get the actual node number,
C     use nodes%atomlist(x).
      
      implicit none
      
C     local variables
      integer :: binlist(2,nodes%natoms)
      real(dp) :: x, y
      integer :: i, j, kx, ky
      integer :: atom, binx, biny, bcount
      integer :: nmaxbin, nbinsx, nbinsy
      integer :: neigh, neighatom, binxneigh, binyneigh
      integer :: neighcount1, neighcount2
      real(dp) :: atompos(2), neighpos(2)
      real(dp) :: dist
      
C     wipe previous lists
      neighbors%neighlist = 0
      neighbors%neighcount = 0
      binlist = 0
      
C     get x, y bounds, store in atombox structure
      call getXYAtomBounds()
      
C     figure out # of bins, fudge end points slightly b/c of finite prec.
      atombox%xmin = atombox%xmin - tolconst
      atombox%ymin = atombox%ymin - tolconst
      nbinsx = ceiling((atombox%xmax - atombox%xmin)/neighbors%rneigh)
      nbinsy = ceiling((atombox%ymax - atombox%ymin)/neighbors%rneigh)
      
C     allocate/initialize arrays
      nmaxbin = ceiling(nmaxbinfac*neighbors%rhomax*
     &                        (neighbors%rneigh*neighbors%rneigh))
      if (allocated(neighbors%binarray)) then
          deallocate(neighbors%binarray)
      end if 
      allocate(neighbors%binarray(nmaxbin,nbinsx,nbinsy))
      neighbors%binarray = 0
      if (allocated(neighbors%bincount)) then
          deallocate(neighbors%bincount)
      end if    
      allocate(neighbors%bincount(nbinsx,nbinsy))
      neighbors%bincount = 0

C     assign to bins
      do i = 1, nodes%natoms
          atom = nodes%atomlist(i)
          x = nodes%posn(1,atom)
          y = nodes%posn(2,atom)
          binx = ceiling((x - atombox%xmin)/neighbors%rneigh)
          biny = ceiling((y - atombox%ymin)/neighbors%rneigh)
          binlist(1,i) = binx
          binlist(2,i) = biny
          bcount = neighbors%bincount(binx,biny) + 1
          if (bcount > nmaxbin) then
              write(*,*) 'nmaxbin is too low'
              write(*,*) 'increase nmaxbinfac'
              stop
          end if    
          neighbors%bincount(binx,biny) = bcount
          neighbors%binarray(bcount,binx,biny) = i
      end do
      
C     use bins to generate neighborlist
      do i = 1,nodes%natoms
          atom = nodes%atomlist(i)
          binx = binlist(1,i)
          biny = binlist(2,i)   
          atompos = nodes%posn(1:2,atom) 
C         loop over neighboring/current bin
          do ky = -1, 1
          binyneigh = biny + ky
          do kx = -1, 1
          binxneigh = binx + kx
C         not on edge
          if ((binxneigh>=1).and.(binxneigh<=nbinsx).and.
     &        (binyneigh>=1).and.(binyneigh<=nbinsy)) then
C             loop through possible neighbors
              do j = 1, neighbors%bincount(binxneigh,binyneigh)
                  neigh = neighbors%binarray(j,binxneigh,binyneigh)
C                 use symmetry (if i is neigh. of j, j is neigh. of i)
                  if (i < neigh) then
                      neighatom = nodes%atomlist(neigh)
                      neighpos = nodes%posn(1:2,neighatom)
                      dist = norm2(atompos-neighpos) 
                      if (dist < neighbors%rneigh) then
                          neighcount1 = neighbors%neighcount(i) + 1
                          neighcount2 = neighbors%neighcount(neigh) + 1
                          if ((neighcount1 > neighbors%nmaxneigh).or.
     &                        (neighcount2 > neighbors%nmaxneigh)) then
                              write(*,*) 'nmaxneigh is too low'
                              write(*,*) 'increase nmaxneighfac'
                              stop
                          end if
                          neighbors%neighlist(neighcount1,i) = neigh
                          neighbors%neighlist(neighcount2,neigh) = i
                          neighbors%neighcount(i) = neighcount1
                          neighbors%neighcount(neigh) = neighcount2
                      end if
                  end if
              end do
          end if
          end do
          end do
      end do
      
      end subroutine genNeighList
************************************************************************
      subroutine getXYAtomBounds()
      
C     Subroutine: getXYAtomBounds

C     Inputs: None

C     Outputs: None

C     Purpose: Get max, min x- and y-positions for all atoms
C     (including pad atoms), store in atombox list

C     output variables
      real(dp) :: xmin, xmax
      real(dp) :: ymin, ymax
      
C     local variables
      integer :: i, atom
      real(dp) :: xcurr, ycurr

C     get x, y bounds by looping over all atoms, including pad atoms
C     could use min, max, but this would result in two loops
      xmin = huge(xmin)
      xmax = -huge(xmax)
      ymin = huge(ymin)
      ymax = -huge(ymax)
      do i = 1, nodes%natoms
          atom = nodes%atomlist(i)
          xcurr = nodes%posn(1,atom)
          ycurr = nodes%posn(2,atom)
          xmin = min(xmin,xcurr)
          xmax = max(xmax,xcurr)
          ymin = min(ymin,ycurr)
          ymax = max(ymax,ycurr) 
      end do
      atombox%xmin = xmin
      atombox%xmax = xmax
      atombox%ymin = ymin
      atombox%ymax = ymax
      
      end subroutine getXYAtomBounds
************************************************************************
      subroutine getAtomsInBoxGroup(xmin,xmax,ymin,ymax,gnum)
      
C     Subroutine: getAtomsInBoxGroup

C     Inputs: xmin, xmax --- x-bounds of box
C             ymin, ymax --- y-bounds of box
C             gnum --- group to assign atoms to

C     Outputs: None

C     Purpose: Get atoms within box, assign them to gnum (used by damped minimization
C     routine in dislocation passing)
      
C     input variables
      real(dp) :: xmin, xmax
      real(dp) :: ymin, ymax
      integer :: gnum
      
C     local variables
      integer :: binxmin, binxmax
      integer :: binymin, binymax
      integer :: j, kx, ky
      integer :: listtemp(nodes%natoms)
      integer :: counter
      integer :: atom
      integer :: node
      real(dp) :: posn(2)
      
      binxmin = ceiling((xmin - atombox%xmin)/neighbors%rneigh)
      binymin = ceiling((ymin - atombox%ymin)/neighbors%rneigh)
      binxmax = ceiling((xmax - atombox%xmin)/neighbors%rneigh)
      binymax = ceiling((ymax - atombox%ymin)/neighbors%rneigh)
      
C     get number of atoms
      counter = 0
      do ky = binymin, binymax
      do kx = binxmin, binxmax
          do j = 1, neighbors%bincount(kx,ky)
              atom = neighbors%binarray(j,kx,ky)
              node = nodes%atomlist(atom)
              posn = nodes%posn(1:2,node)
              if ((posn(1) > xmin).and.(posn(1) < xmax).and.
     &            (posn(2) > ymin).and.(posn(2) < ymax)) then
                  counter = counter + 1
                  listtemp(counter) = node
              end if    
          end do    
      end do
      end do
      
C     allocate groups%members
      if (allocated(groups(gnum)%members)) then
          deallocate(groups(gnum)%members)
      end if
      allocate(groups(gnum)%members(counter))
      
C     assign atoms to group
      groups(gnum)%members = listtemp(1:counter)
      
C     get masks
      call genGroupMaskAtoms(gnum)
      
      end subroutine getAtomsInBoxGroup
************************************************************************
      end module mod_neighbors