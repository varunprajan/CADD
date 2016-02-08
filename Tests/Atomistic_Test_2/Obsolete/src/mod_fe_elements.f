      module mod_fe_elements
      
C     Purpose: Reads/writes/stores information about FE elements.
      
      use mod_utils, only: writeIntMat, prettyPrintIntMat
      use mod_nodes, only: nodes
      use mod_materials, only: nmaterials
      use mod_sort, only: mergeSortCols
      use mod_types, only: dp
      use mod_math, only: getDuplicates, getUniqueInts
      use mod_fe_el_2d, only: getGauss_2d, getGaussEdge_2d
      implicit none
      
      private
      public :: initFEElementData, writeFEElementData, processGauss,
     & readFEElementData, processFEElementData, processEdges,
     & processNodeLists, processMaterialList, feelements, fematerials,
     & nfematerials, processInterfaceEdges, interfaceedges
      
      type feelementdata
C     read-in
      integer, allocatable :: connect(:,:)
      character(len=20) :: elname
      integer :: mnum
C     processed
C     processFEElementData
      integer :: nfeelements
      integer :: nedges
      integer :: nedgenodes
      integer :: nelnodes
      integer :: neldof
      integer :: nelip
      integer :: nedgeip
C     processGaussData
      real(dp), allocatable :: rgauss(:)
      real(dp), allocatable :: sgauss(:)
      real(dp), allocatable :: wgauss(:)
      real(dp), allocatable :: sedgegauss(:)
      real(dp), allocatable :: wedgegauss(:,:)
C     processEdges
      integer, allocatable :: bdedges(:,:)
      integer, allocatable :: neighbors(:,:)
C     processNodeLists
      integer, allocatable :: nodelist(:)
      integer, allocatable :: nodeinvlist(:)
      integer, allocatable :: bdnodelist(:)
      end type
      
      type fematerialdata
C     processed
      integer, allocatable :: list(:)
      integer, allocatable :: invlist(:)
      end type
      
      type interfaceedgesdata
C     processed
      integer, allocatable :: array(:,:)
      end type
      
C     module variables (public)
      type(feelementdata), allocatable :: feelements(:)
      type(fematerialdata) :: fematerials
      integer :: nfematerials
      type(interfaceedgesdata) :: interfaceedges
      
C     module variables (private)
      
      contains
************************************************************************ 
      subroutine initFEElementData(feelementfile)

C     Subroutine: initFEElementData

C     Inputs: feelementfile --- filename where finite element data is stored
C     (should be something like '[filepref]_feelements')

C     Outputs: None

C     Purpose: Read, initialize data in "feelements" structure, which holds
C     information about finite element connectivity, boundary nodes/edges,
C     etc., all of which is used by fe_main_2d to assemble/solve the FE equations

      implicit none
      
C     input variables
      character(len=*) :: feelementfile

      call readFEElementData(feelementfile)
      call processFEElementData()
      call processGauss()
      call processEdges()
      call processInterfaceEdges()
      call processNodeLists()
      call processMaterialList()

      end subroutine initFEElementData
************************************************************************
      subroutine readFEElementData(feelementfile)

C     Subroutine: readFEElementData

C     Inputs: feelementfile --- filename where feelement data is stored
C     (should be something like '[filepref]_feelements')

C     Outputs: None

C     Purpose: Read feelement data (connectivity, etc.) from file,
C     for *each* continuum material, initialize/allocate structures/arrays

      implicit none
      
C     input variables
      character(len=*) :: feelementfile

C     local variables
      integer :: iunit, i, j, k
      integer :: m, nrow
      
      open(newunit=iunit,file=feelementfile)

      read(iunit,*) nfematerials
      allocate(feelements(nfematerials))
      
      do i = 1, nfematerials          
          read(iunit,*) m, nrow
          allocate(feelements(i)%connect(nrow,m))
          read(iunit,*) ((feelements(i)%connect(k,j),k=1,nrow),j=1,m)
          
          read(iunit,*) feelements(i)%elname
          read(iunit,*) feelements(i)%mnum
      end do    
      
      close(iunit)

      end subroutine readFEElementData
************************************************************************
      subroutine processFEElementData()

C     Subroutine: processFEElementData

C     Inputs: None

C     Outputs: None

C     Purpose: Assign properties (# nodes per element, # integration points
C     per element, etc.) to feelements, based on element name

C     Notes: Uses ABAQUS naming scheme for elname

      implicit none

C     local variables
      integer :: i
      
      do i = 1, nfematerials
          feelements(i)%nfeelements = size(feelements(i)%connect,2)
          feelements(i)%nelnodes = size(feelements(i)%connect,1)
          if (feelements(i)%elname == 'CPE3') then
              feelements(i)%nedges = 3
              feelements(i)%nelip = 1
              feelements(i)%neldof = 6
              feelements(i)%nedgenodes = 2
              feelements(i)%nedgeip = 2
          elseif (feelements(i)%elname == 'CPE4') then
              feelements(i)%nedges = 4
              feelements(i)%nelip = 4
              feelements(i)%neldof = 8
              feelements(i)%nedgenodes = 2
              feelements(i)%nedgeip = 2
          else
              write(*,*) 'Unknown element type'
              write(*,*) 'Currently defined elements:'
              write(*,*) '1) CPE3'
              write(*,*) '2) CPE4'
              stop
          end if
      end do

      end subroutine processFEElementData
************************************************************************
      subroutine processGauss()

C     Subroutine: processGauss

C     Inputs: None

C     Outputs: None

C     Purpose: Assign weights, integration points for Gaussian quadrature
C              of line and area integrals
      
C     local variables
      integer :: i
      integer :: nelip, nedgenodes, nedgeip
      
      do i = 1, nfematerials
C         Gauss points for area integral
          nelip = feelements(i)%nelip
          allocate(feelements(i)%rgauss(nelip))
          allocate(feelements(i)%sgauss(nelip))
          allocate(feelements(i)%wgauss(nelip))
          call getGauss_2d(feelements(i)%elname,nelip,
     &  feelements(i)%rgauss,feelements(i)%sgauss,feelements(i)%wgauss)
     
C         Gauss points for line integral
          nedgenodes = feelements(i)%nedgenodes
          nedgeip = feelements(i)%nedgeip
          allocate(feelements(i)%sedgegauss(nedgeip))
          allocate(feelements(i)%wedgegauss(nedgenodes,nedgeip))
          call getGaussEdge_2d(feelements(i)%elname,nedgenodes,nedgeip,
     &      feelements(i)%sedgegauss,feelements(i)%wedgegauss)
      end do
      
      end subroutine processGauss
************************************************************************
      subroutine processEdges()

C     Subroutine: processEdges

C     Inputs: None

C     Outputs: None

C     Purpose: Populate neighbors and boundary edges (bdedges) arrays,
C     using element connectivity

C     Algorithm: Loop over all elements, generate list of nedges*nelements
C     edges. Then, determine unique and duplicate entries of this list by
C     sorting/looping. The unique entries are boundary edges (they are
C     not shared by two elements); the duplicate entries are interior edges,
C     from which we can get the neighbor connectivity
      
      implicit none

C     local variables
      integer :: i, j, k
      integer, allocatable :: edgestemp(:,:)
      integer, allocatable :: intedgestemp(:,:), bdedgestemp(:,:)
      integer :: nfeelements, nedgenodes, nedges, nedgestot, nelnodes
      integer :: counter
      integer :: index2
      integer :: node1, node2
      integer :: reverse
      integer :: rowcheck, rowtot, rowtotnew
      integer :: Nintedges, Nbdedges
      integer :: element1, edge1, element2, edge2
      
      do i = 1, nfematerials
      
C         loop over all edges, populate edgestemp
          nfeelements = feelements(i)%nfeelements
          nedgenodes = feelements(i)%nedgenodes
          nedges = feelements(i)%nedges
          nedgestot = nfeelements*nedges
          nelnodes = feelements(i)%nelnodes
          rowtot = 5
          allocate(edgestemp(rowtot,nedgestot))
          counter = 0
          do j = 1, nfeelements
          do k = 1, nedges
              counter = counter + 1
              if ((feelements(i)%elname == 'CPE3').or.
     &            (feelements(i)%elname == 'CPE4')) then
     
                  node1 = feelements(i)%connect(k,j)                  
                  index2 = mod(k,nelnodes) + 1 ! use trick to loop nelnodes+1 back into 1
                  node2 = feelements(i)%connect(index2,j)
                  
                  if (node1 > node2) then ! reverse nodes
                      edgestemp(1,counter) = node2
                      edgestemp(2,counter) = node1
                      edgestemp(5,counter) = 1 ! flag for reversal
                  else
                      edgestemp(1,counter) = node1
                      edgestemp(2,counter) = node2
                      edgestemp(5,counter) = 0
                  end if
                  edgestemp(3,counter) = j
                  edgestemp(4,counter) = k
              end if
          end do
          end do
          
C         sort edges
          rowcheck = 2
          call mergeSortCols(edgestemp,nedgestot,rowcheck,rowtot)
          
C         get duplicates
          rowtotnew = 2*rowtot - rowcheck
          call getDuplicates(edgestemp,nedgestot,rowcheck,rowtot,
     &            intedgestemp,bdedgestemp,Nintedges,Nbdedges,rowtotnew)
          
C         loop through internal edges, populate neighbor array
          allocate(feelements(i)%neighbors(nedges,nfeelements))
          feelements(i)%neighbors = 0
          do j = 1, Nintedges
              element1 = intedgestemp(3,j)
              edge1 = intedgestemp(4,j)
              element2 = intedgestemp(6,j)
              edge2 = intedgestemp(7,j)
              feelements(i)%neighbors(edge1,element1) = element2
              feelements(i)%neighbors(edge2,element2) = element1
          end do
          
C         loop through bdedges, populate boundary array
          allocate(feelements(i)%bdedges(3,Nbdedges))
          do j = 1, Nbdedges
              node1 = bdedgestemp(1,j)
              node2 = bdedgestemp(2,j)
              reverse = bdedgestemp(5,j)
              if (reverse == 0) then
                  feelements(i)%bdedges(1,j) = node1
                  feelements(i)%bdedges(2,j) = node2
              else
                  feelements(i)%bdedges(1,j) = node2
                  feelements(i)%bdedges(2,j) = node1
              end if    
              feelements(i)%bdedges(3,j) = bdedgestemp(3,j)
          end do
          
          deallocate(edgestemp)
          deallocate(intedgestemp)
          deallocate(bdedgestemp)
      end do
      
      end subroutine processEdges
************************************************************************
      subroutine processInterfaceEdges()

C     Subroutine: processInterfaceEdges

C     Inputs: None

C     Outputs: None

C     Purpose: Generate array of boundary facets for interface between
C     atomistic and continuum region: each column is a facet;
C     rows are [node i, node j, mnumfe, element]

      implicit none
      
C     local variables
      integer :: i, j
      integer :: nmaxbdedges
      integer :: counter
      integer, allocatable :: interfaceedgestemp(:,:)
      integer :: node1, node2
      integer :: element
      integer :: nodetype1, nodetype2
      
C     figure out storage requirements
      nmaxbdedges = 0
      do i = 1, nfematerials
          nmaxbdedges = nmaxbdedges + size(feelements(i)%bdedges,2)
      end do
      
C     loop over bdedges
      allocate(interfaceedgestemp(4,nmaxbdedges))
      counter = 0
      do i = 1, nfematerials
          do j = 1, size(feelements(i)%bdedges,2)
              node1 = feelements(i)%bdedges(1,j)
              node2 = feelements(i)%bdedges(2,j)
              element = feelements(i)%bdedges(3,j)
              nodetype1 = nodes%types(2,node1)
              nodetype2 = nodes%types(2,node2)
              if ((nodetype1==2).and.(nodetype2==2)) then ! both are interface nodes
                  counter = counter + 1
                  interfaceedgestemp(1,counter) = node1
                  interfaceedgestemp(2,counter) = node2
                  interfaceedgestemp(3,counter) = i
                  interfaceedgestemp(4,counter) = element
              end if
          end do    
      end do
      
C     assign to array
      allocate(interfaceedges%array(4,counter))
      interfaceedges%array = interfaceedgestemp(:,1:counter)
      
      end subroutine processInterfaceEdges
************************************************************************
      subroutine processNodeLists()

C     Subroutine: processNodeLists

C     Inputs: None

C     Outputs: None

C     Purpose: Generate lists of boundary nodes, all nodes, as well as
C     "inverse" lists (which basically allow us to go "backwards")

      implicit none

C     local variables
      integer :: i
      integer :: counter
      integer :: invlisttemp(nodes%nnodes)
      
C     This is not efficient from the viewpoint of storage, but makes
C     assembly of FE matrices easier. Would be better to use
C     numbering with nodes%fenodes
      
      do i = 1, nfematerials
          allocate(feelements(i)%nodeinvlist(nodes%nnodes))
          call getUniqueInts(feelements(i)%connect,nodes%nnodes,
     &        counter,feelements(i)%nodelist,feelements(i)%nodeinvlist)
          
          call getUniqueInts(feelements(i)%bdedges,nodes%nnodes,
     &                     counter,feelements(i)%bdnodelist,invlisttemp)
          
      end do

      end subroutine processNodeLists
************************************************************************
      subroutine processMaterialList()

C     Subroutine: processMaterialList

C     Inputs: None

C     Outputs: None

C     Purpose: Generate forwards and backwards lists, relating the
C     numbering of the fematerials to that of the overall materials
C     (in the materials struct)
      
C     local variables
      integer :: i, mnum
      
      allocate(fematerials%list(nfematerials))
      allocate(fematerials%invlist(nmaterials))
      fematerials%invlist = 0
      do i = 1, nfematerials
          mnum = feelements(i)%mnum
          fematerials%list(i) = mnum
          fematerials%invlist(mnum) = i
      end do
      
      end subroutine processMaterialList
************************************************************************
      subroutine writeFEElementData(feelementfile)

C     Subroutine: writeFEElementData

C     Inputs: feelementfile --- filename where finite element data is stored
C     (should be something like '[filepref]_feelements')

C     Outputs: None

C     Purpose: Write finite element data to file (essentially
C     inverse of readFEElementData). Useful in creating "restart" file

      implicit none
      
C     input variables
      character(len=*) :: feelementfile

C     local variables
      integer :: iunit, i
      
      open(newunit=iunit,file=feelementfile)
      
      write(iunit,*) nfematerials
      
      do i = 1, nfematerials
          call writeIntMat(feelements(i)%connect,iunit,.true.)
          write(iunit,*) feelements(i)%elname
          write(iunit,*) feelements(i)%mnum
      end do
      
      close(iunit)

      end subroutine writeFEElementData
************************************************************************
      end module