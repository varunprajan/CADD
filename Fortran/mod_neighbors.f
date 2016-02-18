      module mod_neighbors

C     Purpose: Has routines for generating neighbor lists, and checking
C     displacements to see if neighbor lists should be regenerated.
C     Also contains information related to neighbors (cutoff distance,
C     etc.)

C     Possible extensions: ?

      use mod_types, only: dp
      use mod_math, only: tolconst, piconst, isMultiple,
     &                    logicalToInt, intToLogical
      use mod_nodes, only: nodes, getXYAtomBounds
      use mod_materials, only: materials, nmaterials
      use mod_potentials, only: potentials, npotentials
      use mod_groups, only: groups, genGroupMaskAtoms, tempgroupname,
     &                      getGroupNum
      use mod_sort, only: mergeSort
      implicit none
      
      private
      public :: neighbors, updateNeighborsCheck, initNeighborData,
     &          writeNeighborData, genNeighList, checkDisp,
     &          updatePosSinceLastCheck, readNeighborData, genAtomBins,
     &          processNeighborData, updateNeighIncrementCurr,
     &          getAtomsInBoxGroupTemp, updateNeighborsNoCheck
      
      type neighbordata
C     read-in
      logical :: checkdisp
      integer :: delay
      integer :: every
      integer :: images
      real(dp) :: Lz
      real(dp) :: skin
C     processed
      integer :: incrementcurr
      integer ,allocatable :: binlist(:,:)
      integer, allocatable :: bincount(:,:)
      integer, allocatable :: binarray(:,:,:)
      integer, allocatable :: neighlist(:,:)
      integer, allocatable :: neighcount(:)
      integer :: nmaxneigh
      integer :: nmaxbin
      integer :: incrsincelastupdate
      real(dp), allocatable :: possincelastcheck(:,:)
      real(dp) :: rhomax
      real(dp) :: rneigh
      real(dp) :: rneighsq
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
      neighbors%checkdisp = intToLogical(temp)
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
      
      implicit none
      
C     local variables
      integer :: i
      
C     get rhomax (used to calculate nmaxbin for building bins)
      neighbors%rhomax = 0.0_dp
      do i = 1, nmaterials
          if (materials(i)%rho > neighbors%rhomax) then
              neighbors%rhomax = materials(i)%rho
          end if
      end do

C     get rneigh (for calculating neighbor lists)
C     first calculate largest force cutoff, then add skin distance
      neighbors%rneigh = 0.0_dp
      do i = 1, npotentials
          if (potentials(i)%forcecutoff > neighbors%rneigh) then
              neighbors%rneigh = potentials(i)%forcecutoff
          end if
      end do
      neighbors%rneigh = neighbors%rneigh + neighbors%skin
      neighbors%rneighsq = (neighbors%rneigh)**2

C     calculate nmaxneigh using area of circle, density,
C     factor of safety of nmaxneighfac, then allocate neighlist
      neighbors%nmaxneigh = ceiling(nmaxneighfac*neighbors%rhomax*
     &                              piconst*neighbors%rneighsq)
      allocate(neighbors%neighlist(neighbors%nmaxneigh,nodes%natoms))
      
C     similar for nmaxbin
      neighbors%nmaxbin = ceiling(nmaxbinfac*neighbors%rhomax*
     &                            neighbors%rneighsq)
      
C     allocate other arrays
      allocate(neighbors%neighcount(nodes%natoms))
      allocate(neighbors%possincelastcheck(2,nodes%natoms))
      allocate(neighbors%binlist(2,nodes%natoms)) 

C     initialize neighbor list (this takes care of timestep initialization, etc.)
      call updateNeighborsNoCheck()
      
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
      temp = logicalToInt(neighbors%checkdisp)  
      write(iunit,*) temp
      write(iunit,*) neighbors%delay
      write(iunit,*) neighbors%every
      write(iunit,*) neighbors%images
      write(iunit,*) neighbors%Lz
      write(iunit,*) neighbors%skin
      
      close(iunit)
      
      end subroutine writeNeighborData
************************************************************************
      subroutine updateNeighIncrementCurr(deltaincrement)
      
C     Subroutine: updateNeighIncrementCurr

C     Inputs: deltaincrement --- change in incrementcurr (incrementcurrnew - incrementcurrold)

C     Outputs: None

C     Purpose: Update current increment...usually, usage is updateIncrementCurr(1), to advance by 1 step
      
      implicit none
      
C     input variables
      integer :: deltaincrement
      
      neighbors%incrementcurr = neighbors%incrementcurr + deltaincrement
      
      end subroutine updateNeighIncrementCurr
************************************************************************
      subroutine updateNeighborsCheck()

C     Subroutine: updateNeighborsCheck

C     Outputs: None

C     Purpose: Regenerates neighbor list
C     *after checking* "check", "delay", and "every" conditions.
      
      implicit none
      
C     local variables
      logical :: update
      integer :: delayincr
      
C     check delay
      delayincr = neighbors%incrementcurr -
     &             neighbors%incrsincelastupdate
      update = (delayincr >= neighbors%delay)
      
C     check every
      if (update) then
          if (neighbors%every /= 0) then ! disregard check if every = 0
              update = isMultiple(delayincr,neighbors%every)
          end if    
      end if    
      
C     if check is on, need to check displacements
      if (update) then
          if (neighbors%checkdisp) then
              update = checkDisp()
          end if
      end if    
      
C     if we need to update, then regenerate neighbor lists,
C     and update displacements since last check array
      if (update) then
          call updateNeighborsNoCheck()
      end if
      
      end subroutine updateNeighborsCheck
************************************************************************
      subroutine updateNeighborsNoCheck()

C     Subroutine: updateNeighborsNoCheck

C     Outputs: None

C     Purpose: Regenerates neighbor list, pos since last check (for displacement check)
C     *without checking* "check", "delay", and "every" conditions.

      implicit none
      
      write(*,*) 'Updated neighbors'
      call genAtomBins()
      call genNeighList()
      call updatePosSinceLastCheck()
      neighbors%incrementcurr = 0
      neighbors%incrsincelastupdate = neighbors%incrementcurr
      
      end subroutine
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
      integer :: atom
      
      do i = 1,nodes%natoms
          atom = nodes%atomlist(i)
          neighbors%possincelastcheck(:,i) = nodes%posn(1:2,atom)
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
      real(dp) :: dispsqmax1, dispsqmax2, dispsq
      real(dp) :: posold(2), posnew(2)
      
C     goal: get two largest displacements since last check
      dispsqmax1 = 0.0_dp
      dispsqmax2 = 0.0_dp
      do i = 1,nodes%natoms
          atom = nodes%atomlist(i)
          posold = neighbors%possincelastcheck(:,i)
          posnew = nodes%posn(1:2,atom)
          dispsq = sum((posold-posnew)**2)
          if (dispsq > dispsqmax2) then
              if (dispsq > dispsqmax1) then
                  dispsqmax2 = dispsqmax1
                  dispsqmax1 = dispsq
              else
                  dispsqmax2 = dispsq
              end if
          end if
      end do
      dispcheck = (sqrt(dispsqmax1) + sqrt(dispsqmax2) > neighbors%skin)
      
      end function checkDisp
************************************************************************
      subroutine genAtomBins()
      
C     Subroutine: genAtomBins

C     Inputs: None

C     Outputs: None

C     Purpose: Assign atoms to bins (see binlist, binarray, bincount)

C     Notes: Binning is efficient only when atoms are roughly
C     contiguous (otherwise, bins are too large).

      implicit none
      
C     local variables
      integer :: nbinsx, nbinsy
      integer :: i
      integer :: atom
      real(dp) :: x, y
      integer :: binx, biny, bcount
      
C     get x, y bounds, store in atombox structure
      call getXYAtomBounds(atombox%xmin,atombox%xmax,
     &                     atombox%ymin,atombox%ymax)
      
C     figure out # of bins, fudge end points slightly b/c of finite prec.
      atombox%xmin = atombox%xmin - tolconst
      atombox%ymin = atombox%ymin - tolconst
      nbinsx = ceiling((atombox%xmax - atombox%xmin)/neighbors%rneigh)
      nbinsy = ceiling((atombox%ymax - atombox%ymin)/neighbors%rneigh)
      
C     (re)allocate/initialize arrays
      if (allocated(neighbors%binarray)) then
          deallocate(neighbors%binarray)
      end if 
      allocate(neighbors%binarray(neighbors%nmaxbin,nbinsx,nbinsy))
      neighbors%binarray = 0
      if (allocated(neighbors%bincount)) then
          deallocate(neighbors%bincount)
      end if    
      allocate(neighbors%bincount(nbinsx,nbinsy))
      neighbors%bincount = 0
      
C     assign atoms to bins
      neighbors%binlist = 0
      do i = 1, nodes%natoms
          atom = nodes%atomlist(i)
          x = nodes%posn(1,atom)
          y = nodes%posn(2,atom)
          binx = ceiling((x - atombox%xmin)/neighbors%rneigh)
          biny = ceiling((y - atombox%ymin)/neighbors%rneigh)
          neighbors%binlist(1,i) = binx
          neighbors%binlist(2,i) = biny
          bcount = neighbors%bincount(binx,biny) + 1
          if (bcount > neighbors%nmaxbin) then
              write(*,*) 'nmaxbin is too small'
              write(*,*) 'increase nmaxbinfac'
              stop
          end if    
          neighbors%bincount(binx,biny) = bcount
          neighbors%binarray(bcount,binx,biny) = i
      end do
      
      end subroutine genAtomBins
************************************************************************
      subroutine genNeighList()
      
C     Subroutine: genNeighList

C     Inputs: None

C     Outputs: None

C     Purpose: Generate neighbor list, given atom positions
      
C     Algorithm: Loop through binlist to get neighbors.
C     Exploit symmetry to avoid duplicate work (i.e. use fact that if i
C     is a neighbor of j, j is a neighbor of i)
      
C     Notes: The result (neighbors%neighlist) contains entries
C     that are indices of nodes%atomlist. To get the actual node number,
C     use nodes%atomlist(i).
      
      implicit none
      
C     local variables
      integer :: i, j, kx, ky
      integer :: atom, binx, biny
      integer :: nbinsx, nbinsy
      integer :: neigh, neighatom, binxneigh, binyneigh
      integer :: neighcount1, neighcount2
      real(dp) :: atompos(2), neighpos(2)
      real(dp) :: distsq     
      
C     use bins to generate neighborlist      
      nbinsx = size(neighbors%bincount,1)
      nbinsy = size(neighbors%bincount,2)
      neighbors%neighlist = 0
      neighbors%neighcount = 0
      do i = 1,nodes%natoms
          atom = nodes%atomlist(i)
          binx = neighbors%binlist(1,i)
          biny = neighbors%binlist(2,i)   
          atompos = nodes%posn(1:2,atom) 
C         loop over neighboring/current bin
          do ky = -1, 1
          binyneigh = biny + ky
          do kx = -1, 1
          binxneigh = binx + kx
C         not outside box
          if ((binxneigh>=1).and.(binxneigh<=nbinsx).and.
     &        (binyneigh>=1).and.(binyneigh<=nbinsy)) then
C             loop through possible neighbors
              do j = 1, neighbors%bincount(binxneigh,binyneigh)
                  neigh = neighbors%binarray(j,binxneigh,binyneigh)
C                 use symmetry (if i is neigh. of j, j is neigh. of i)
                  if (i < neigh) then
                      neighatom = nodes%atomlist(neigh)
                      neighpos = nodes%posn(1:2,neighatom)
                      distsq = sum((atompos-neighpos)**2)
                      if (distsq < neighbors%rneighsq) then
                          neighcount1 = neighbors%neighcount(i) + 1
                          neighcount2 = neighbors%neighcount(neigh) + 1
                          if ((neighcount1 > neighbors%nmaxneigh).or.
     &                        (neighcount2 > neighbors%nmaxneigh)) then
                              write(*,*) 'nmaxneigh is too small'
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
      subroutine getAtomsInBoxGroupTemp(xmin,xmax,ymin,ymax)
      
C     Subroutine: getAtomsInBoxGroup

C     Inputs: xmin, xmax --- x-bounds of box
C             ymin, ymax --- y-bounds of box

C     Outputs: None

C     Purpose: Get atoms within box, assign them to temp group (used by damped minimization
C     routine in dislocation passing)
      
      implicit none
      
C     input variables
      real(dp) :: xmin, xmax
      real(dp) :: ymin, ymax
      
C     local variables
      integer :: binxmin, binxmax
      integer :: binymin, binymax
      integer :: nbinsx, nbinsy
      integer :: j, kx, ky
      integer :: listtemp(nodes%natoms)
      integer :: counter
      integer :: atom
      integer :: node
      real(dp) :: posn(2)
      integer :: gnum
      
      ! regenerate atom bins...
      call genAtomBins()
      
      binxmin = ceiling((xmin - atombox%xmin)/neighbors%rneigh)
      binymin = ceiling((ymin - atombox%ymin)/neighbors%rneigh)
      binxmax = ceiling((xmax - atombox%xmin)/neighbors%rneigh)
      binymax = ceiling((ymax - atombox%ymin)/neighbors%rneigh)
      
      ! make sure box doesn't extend past atom box
      nbinsx = size(neighbors%bincount,1)
      nbinsy = size(neighbors%bincount,2)
      binxmin = max(binxmin,1)
      binymin = max(binymin,1)
      binxmax = min(binxmax,nbinsx)
      binymax = min(binymax,nbinsy)
      
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
      gnum = getGroupNum(tempgroupname)
      if (allocated(groups(gnum)%members)) then
          deallocate(groups(gnum)%members)
      end if
      
C     assign atoms to group (automatic allocation), sort
      groups(gnum)%members = listtemp(1:counter)
      call mergeSort(groups(gnum)%members,counter)
      
C     get/regenerate masks
      call genGroupMaskAtoms(gnum)
      
      end subroutine getAtomsInBoxGroupTemp
************************************************************************
      end module mod_neighbors