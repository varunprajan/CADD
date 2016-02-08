      module mod_groups
      
C     Purpose: Reads/writes/stores information about groups, including
C     their members and "masks", which are logical arrays that tell subroutines
C     which atoms/nodes/fenodes should be included. For example, typical usage
C     of a mask for group gnum for an atom loop is:
C     do i, nodes%natoms
C         if (groups(gnum)%maskatoms(i)) then
C             [do stuff]
C         end if
C     end do
      
      use mod_types, only: dp
      use mod_nodes, only: nodes
      implicit none
      
      private
      public :: groups, ngroups, initGroupData, writeGroupData,
     &          readGroupData, processGroupData, getMaskSortedIntersect,
     &          genGroupMaskAtoms, genGroupMaskFENodes, genGroupMaskAll,
     &          tempgroupnum, nextragroups
      
      type groupdata
C     read-in
      character(len=20) :: gname
      integer, allocatable :: members(:)
      integer :: nmembers
C     processed
      logical, allocatable :: maskatoms(:)
      logical, allocatable :: maskfenodes(:)
      logical, allocatable :: maskall(:)
      end type
      
C     module variables (global)
      type(groupdata), allocatable :: groups(:)
      integer :: ngroups
      integer, parameter :: nextragroups = 2
      integer, parameter :: allgroupnum = 1
      integer, parameter :: tempgroupnum = 2
      
      contains
************************************************************************          
      subroutine initGroupData(groupfile)

C     Subroutine: initGroupData

C     Inputs: groupfile --- filename where group data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read, initialize data in "groups" structure, which holds
C     information about members in groups, and group "masks."
      
      implicit none
      
C     input variables
      character(len=*) :: groupfile
      
      call readGroupData(groupfile)
      call processGroupData()
      
      end subroutine initGroupData
************************************************************************          
      subroutine readGroupData(groupfile)
      
C     Subroutine: readGroupData

C     Inputs: groupfile --- filename where group data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read group data (lists of members of groups) from file.
      
      implicit none
      
C     input variables
      character(len=*) :: groupfile
      
C     local variables
      integer :: iunit, i, j
      
      open(newunit=iunit,file=groupfile)
      
      read(iunit,*) ngroups
C     extra groups for "all" and "temp" groups
      ngroups = ngroups + 2
      allocate(groups(ngroups))
C     read non-"all" groups
      do i = 3, ngroups
          read(iunit,*) groups(i)%gname
          read(iunit,*) groups(i)%nmembers
          allocate(groups(i)%members(groups(i)%nmembers))
          read(iunit,*) (groups(i)%members(j),j=1,groups(i)%nmembers)
      end do
      
      close(iunit)
      
      end subroutine readGroupData
************************************************************************      
      subroutine processGroupData()

C     Subroutine: processGroupData

C     Inputs: None

C     Outputs: None

C     Purpose: Create "all" group (group containing all nodes), and
C     create "masks" for each group.
      
C     local variables
      integer :: i
      logical :: createoption(3)
      
C     first, create "all" group
      groups(allgroupnum)%gname = 'all'
      allocate(groups(allgroupnum)%members(nodes%nnodes))
      do i = 1, nodes%nnodes
          groups(allgroupnum)%members(i) = i
      end do
      
C     second, initialize "temp" group
      groups(tempgroupnum)%gname = 'temp'
      allocate(groups(tempgroupnum)%members(1))
      
C     then, create masks
      createoption = [.true.,.true.,.true.]
      do i = 1, ngroups
          allocate(groups(i)%maskatoms(nodes%natoms))
          call genGroupMaskAtoms(i)
          allocate(groups(i)%maskfenodes(nodes%nfenodes))
          call genGroupMaskFENodes(i)
          allocate(groups(i)%maskall(nodes%nnodes))
          call genGroupMaskAll(i)
      end do
      
      end subroutine processGroupData
************************************************************************
      subroutine genGroupMaskAtoms(gnum)
      
C     Subroutine: genGroupMasksAtoms

C     Inputs: gnum --- number of group

C     Outputs: None

C     Purpose: Create "mask" for atom group
      
C     input variables
      integer :: gnum
      
      groups(gnum)%maskatoms =
     &     getMaskSortedIntersect(groups(gnum)%members,nodes%atomlist)
     
      end subroutine genGroupMaskAtoms
************************************************************************      
      subroutine genGroupMaskFENodes(gnum)

C     Subroutine: genGroupMaskFENodes

C     Inputs: gnum --- number of group

C     Outputs: None

C     Purpose: Create "mask" for fe node group
      
C     input variables
      integer :: gnum
      
      groups(gnum)%maskfenodes =
     &     getMaskSortedIntersect(groups(gnum)%members,nodes%fenodelist)
     
      end subroutine genGroupMaskFENodes
************************************************************************
      subroutine genGroupMaskAll(gnum)

C     Subroutine: genGroupMaskAll

C     Inputs: gnum --- number of group

C     Outputs: None

C     Purpose: Create "mask" for all atoms/nodes
      
C     input variables
      integer :: gnum
      
      groups(gnum)%maskall =
     &    getMaskSortedIntersect(groups(gnum)%members,groups(1)%members)
     
      end subroutine genGroupMaskAll
************************************************************************
      subroutine writeGroupData(groupfile)
      
C     Subroutine: writeGroupData

C     Inputs: groupfile --- filename where group data is stored
C     (should be something like '[filepref]_groups')

C     Outputs: None

C     Purpose: Write group data to file (essentially
C     inverse of readGroupData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: groupfile      
      
C     local variables
      integer :: iunit, i, j
      
      open(newunit=iunit,file=groupfile)
      
C     write non-"all"/"temp" groups
      write(iunit,*) ngroups-2
      do i = 3, ngroups
          write(iunit,*) groups(i)%gname
          write(iunit,*) groups(i)%nmembers
          do j = 1, groups(i)%nmembers
              write(iunit,*) groups(i)%members(j)
          end do
      end do
      
      close(iunit)
      
      end subroutine writeGroupData
************************************************************************
      function getMaskSortedIntersect(list1,list2) result(mask2)
      
C     Function: getMaskSortedIntersect

C     Inputs: list1 - sorted m by 1 integer array
C             list2 - sorted n by 1 integer array

C     Outputs: mask2 - n by 1 logical array, where list2(mask2) =
C     intersect(list1,list2) in MATLAB-style "index" notation.

C     Purpose: Generate "masks" for groups. Typical usage is:
C     maskatoms = getMaskSortedIntersect(group(gnum)%members,nodes%atomlist)

C     Notes: Can be easily extended to get the intersection list, and/or to
C     get mask1
      
      implicit none
      
C     input variables
      integer :: list1(:)
      integer :: list2(:)
      
C     output variables
      logical, allocatable :: mask2(:)
      
C     local variables
      integer :: i, j
      integer :: m, n
      integer :: list1el, list2el
      
      i = 1
      j = 1
      m = size(list1)
      n = size(list2)
      allocate(mask2(n))
      mask2 = .false.
      do while ((i <= m).and.(j <= n))
          list1el = list1(i)
          list2el = list2(j)
          if (list1el == list2el) then
              mask2(j) = .true.
              i = i + 1
              j = j + 1
          else if (list1el > list2el) then
              j = j + 1
          else
              i = i + 1
          end if
      end do
      
      end function getMaskSortedIntersect
************************************************************************      
      end module