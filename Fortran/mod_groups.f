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

C     Note: Members of group are assumed to be sorted!!
      
      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_utils, only: readVecSize, writeVecSize
      implicit none
      
      private
      public :: groups, ngroups, initGroupData, writeGroupData,
     &          readGroupData, processGroupData, getMaskSortedIntersect,
     &          genGroupMaskAtoms, genGroupMaskFENodes, genGroupMaskAll,
     &          tempgroupname, getGroupNum, allgroupname
      
      type groupdata
C     read-in
      character(len=20) :: gname
      integer, allocatable :: members(:)
C     processed
      integer :: nmembers
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
      character(len=*), parameter :: tempgroupname = 'temp'
      character(len=*), parameter :: allgroupname = 'all'
      
      contains
************************************************************************          
      subroutine initGroupData(groupfile)

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

C     Inputs: groupfile --- filename where group data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read group data (lists of members of groups) from file.
      
      implicit none
      
C     input variables
      character(len=*) :: groupfile
      
C     local variables
      integer :: iunit, i
      
      open(newunit=iunit,file=groupfile)
      
      read(iunit,*) ngroups
C     extra groups for "all" and "temp" groups
      ngroups = ngroups + nextragroups
      allocate(groups(ngroups))
C     read non-"all" groups
      do i = nextragroups+1, ngroups
          read(iunit,*) groups(i)%gname
          call readVecSize(iunit,groups(i)%members)
          groups(i)%nmembers = size(groups(i)%members)
      end do
      
      close(iunit)
      
      end subroutine readGroupData
************************************************************************      
      subroutine processGroupData()

C     Inputs: None

C     Outputs: None

C     Purpose: Create "all" group (group containing all nodes), and
C     create "masks" for each group.
      
C     local variables
      integer :: i
      
C     first, create "all" group
      groups(allgroupnum)%gname = allgroupname
      allocate(groups(allgroupnum)%members(nodes%nnodes))
      do i = 1, nodes%nnodes
          groups(allgroupnum)%members(i) = i
      end do
      
C     second, initialize "temp" group
      groups(tempgroupnum)%gname = tempgroupname
      allocate(groups(tempgroupnum)%members(1))
      
C     then, create masks
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
      function getGroupNum(gname) result(gnum)

C     Inputs: gname --- name of group

C     Outputs: gnum --- index/number of group (1 <= gnum <= ngroups)

C     Purpose: From group name, get group index/number

      implicit none      
      
C     input variables
      character(len=*) :: gname
      
C     output variables
      integer :: gnum
      
C     local variables
      integer :: i
      
      do i = 1, ngroups
          if (trim(gname) == trim(groups(i)%gname)) then
              gnum = i
              return
          end if
      end do
      
      write(*,*) 'No group with name ', gname, ' was found.'
      stop
      
      end function getGroupNum
************************************************************************
      subroutine genGroupMaskAtoms(gnum)

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

C     Inputs: groupfile --- filename where group data is stored
C     (should be something like '[filepref]_groups')

C     Outputs: None

C     Purpose: Write group data to file (essentially
C     inverse of readGroupData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: groupfile      
      
C     local variables
      integer :: iunit, i
      
      open(newunit=iunit,file=groupfile)
      
C     write non-"all"/"temp" groups
      write(iunit,*) ngroups-2
      do i = 3, ngroups
          write(iunit,*) groups(i)%gname
          call writeVecSize(iunit,groups(i)%members)
      end do
      
      close(iunit)
      
      end subroutine writeGroupData
************************************************************************
      function getMaskSortedIntersect(list1,list2) result(mask2)

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