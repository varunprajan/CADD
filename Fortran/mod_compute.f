      module mod_compute
      
C     Purpose: Contains routines for computing quantities; for instance,
C     centrosymmetry of every atom. Each quantity can be computed for a particular
C     group and using user-defined parameters (real numbers)
      
C     NOTES/TODO: Implement additional computes
C     (only centrosymmetry has been implemented so far)
      
C     NOTES/TODO: Somewhat type-limited: params/result are real arrays,
C     so integer parameters must be converted to reals first. Also,
C     result that is a number (not an array) can be stored as first element of array...

      use mod_types, only: dp
      use mod_neighbors, only: neighbors, updateNeighborsNoCheck
      use mod_nodes, only: nodes
      use mod_sort, only: mergeSortReal, mergeSortColsReal
      use mod_utils, only: readVecSize, writeVecSize
      use mod_math, only: intToLogical, logicalToInt
      use mod_groups, only: groups, getGroupNum
      
      implicit none
      
      private
      public :: initComputeData, getCentroAtoms, getCentroAtom, compute,
     &          readComputeData, writeComputeData, processComputeData
      
      type cdata
C     read-in
      logical :: active
      character(len=60) :: gname
      real(dp), allocatable :: params(:) ! parameters are real, so integer conversion must be done explicitly
C     processed
      real(dp), allocatable :: res(:)
      end type
      
      type computedata
      type(cdata) :: centro
C     can add more computes here as needed; need to add corresponding line to read, process, write
      end type
      
      type(computedata) :: compute
      
      contains
************************************************************************
      subroutine initComputeData(computefile)
      
C     Subroutine: initComputeData

C     Inputs: computefile --- filename where compute data is stored
C     (should be something like '[filepref]_compute')

C     Outputs: None

C     Purpose: Read, initialize data in "compute" structure, which holds
C     information related to computed quantities (such as centrosymmetry)

      implicit none
      
C     input variables
      character(len=*) :: computefile
      
      call readComputeData(computefile)
      call processComputeData()
      
      end subroutine initComputeData
************************************************************************
      subroutine readComputeData(computefile)
      
C     Subroutine: readComputeData

C     Inputs: dislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read all compute data (group name, parameters, etc.) from file

      implicit none
      
C     input variables
      character(len=*) :: computefile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=computefile)      
      
      call readComputeDataSub(iunit,compute%centro)
C     add additional line for each additional compute
      
      close(iunit)      
      
      end subroutine readComputeData
************************************************************************
      subroutine readComputeDataSub(iunit,comp)
      
C     Subroutine: readComputeDataSub

C     Inputs: iunit --- integer file specifier
C             comp --- structure (type cdata) containing data for individual compute

C     Outputs: None

C     Purpose: Read data for a specific compute from file
      
C     input variables
      integer :: iunit
      type(cdata) :: comp
      
C     local variables
      integer :: temp
      
      read(iunit,*) temp
      comp%active = intToLogical(temp)
      read(iunit,*) comp%gname
      call readVecSize(iunit,comp%params)
      
      end subroutine readComputeDataSub
************************************************************************
      subroutine processComputeData()

C     Subroutine: processComputeData

C     Inputs: None

C     Outputs: None

C     Purpose: Allocate results vector for computes
      
      allocate(compute%centro%res(nodes%natoms)) ! centro is only applicable to atoms
C     add additional line for each additional compute
      
      end subroutine processComputeData
************************************************************************
      subroutine writeComputeData(computefile)
      
C     Subroutine: writeComputeData

C     Inputs: computefile --- filename where compute data is stored
C     (should be something like '[filepref]_compute')

C     Outputs: None

C     Purpose: Write compute data to file (essentially
C     inverse of readComputeData). Useful in creating "restart" file

      implicit none
      
C     input variables
      character(len=*) :: computefile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=computefile)
      
      call writeComputeDataSub(iunit,compute%centro)
C     add additional line for each additional compute
      
      close(iunit)
      
      end subroutine writeComputeData
************************************************************************
      subroutine writeComputeDataSub(iunit,comp)
      
C     Subroutine: writeComputeDataSub

C     Inputs: iunit --- integer file specifier
C             comp --- structure (type cdata) containing data for individual compute

C     Outputs: None

C     Purpose: Write data for a specific compute to file
      
C     input variables
      integer :: iunit
      type(cdata) :: comp
      
C     local variables
      integer :: temp
      
      temp = logicalToInt(comp%active)
      write(iunit,*) temp
      write(iunit,*) comp%gname
      call writeVecSize(iunit,comp%params)      
      
      end subroutine writeComputeDataSub
************************************************************************
      subroutine getCentroAtoms()
      
C     Subroutine: getCentroAtoms

C     Inputs: None

C     Outputs: None

C     Purpose: Compute centrosymmetry for all atoms, after updating
C     neighbor list first.
      
      implicit none
      
C     local variables
      integer :: i
      integer :: gnum
      integer :: n
      
C     update neighbor list
      call updateNeighborsNoCheck()
      
      compute%centro%res = 0.0_dp
      gnum = getGroupNum(compute%centro%gname)
      n = nint(compute%centro%params(1)) ! number of nearest neighbors
C     loop over neighborlist
      do i = 1, nodes%natoms
          if (groups(gnum)%maskatoms(i)) then
              compute%centro%res(i) = getCentroAtom(i,n)
          end if    
      end do    
      
      end subroutine getCentroAtoms
************************************************************************
      function getCentroAtom(i,n) result(centro)
      
C     Function: getCentroAtom

C     Inputs: i --- index of atom in nodes%atomlist
C             n --- number of nearest neighbors to use in calculation

C     Outputs: centro --- centrosymmetry of atom

C     Purpose: Compute centrosymmetry for an individual atom, using
C     n nearest neighbors. Uses algorithm described in LAMMPS manual.
C     (http://lammps.sandia.gov/doc/compute_centro_atom.html)

C     input variables
      integer :: i
      integer :: n
      
C     output variables
      real(dp) :: centro
      
C     local variables
      integer :: j, k
      real(dp), allocatable :: temp(:)
      integer :: neighcount
      integer :: atom, atomj
      integer :: neighj
      real(dp) :: atompos(2), atomposj(2)
      real(dp), allocatable :: rmat(:,:)
      real(dp) :: rj(2), rk(2)
      real(dp) :: rjnormsq
      integer :: counter
      integer :: ntemp
      
C     check to see if there are enough neighbors
      neighcount = neighbors%neighcount(i)
      if (neighcount < n) then ! if not enough neighbors, = 0
          centro = 0.0_dp
          return
      end if
      
C     create rvecs
      allocate(rmat(3,neighcount))
      atom = nodes%atomlist(i)
      atompos = nodes%posn(1:2,atom)
      do j = 1, neighcount
          neighj = neighbors%neighlist(j,i)
          atomj = nodes%atomlist(neighj)
          atomposj = nodes%posn(1:2,atomj)
          rj = atomposj - atompos
          rjnormsq = sum(rj**2)
          rmat(1,j) = rjnormsq
          rmat(2:3,j) = rj
      end do
      
C     get n shortest (n closest neighbors) by sorting
      call mergeSortColsReal(rmat,neighcount,1,3)

C     allocate, initialize temp (stores intermediate results)
      ntemp = n*(n-1)/2
      allocate(temp(ntemp))
      temp = huge(0.0_dp) ! see sort, below
      
C     loop over pairs of rvecs
      counter = 0
      do j = 1, n
          rj = rmat(2:3,j)
          do k = 1, j - 1 ! symmetry
              rk = rmat(2:3,k)
              counter = counter + 1
              temp(counter) = sum((rj + rk)**2)
          end do
      end do
      
C     sort and add n/2 smallest
      call mergeSortReal(temp,ntemp)
      centro = sum(temp(1:n/2))
      
      end function getCentroAtom
************************************************************************    
      end module mod_compute