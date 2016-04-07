      module mod_find_crack_atomistic
      
C     Purpose: Finds crack position in atomistic region in a CADD simulation.
C     (Can be used in an atomistic only simulation, but the outer nodes need to be
C     labelled as interface nodes.)
C     Strategy: We detect the crack using geometry:
C     we triangulate the atoms, and determine the locations of
C     "large" triangles that are not at the edge of the region. The edge of the region
C     is associated with triangles that comprise one or more interface nodes. The "rightmost"
C     large triangle not on the edge corresponds to the crack tip

C     TODO: Modify for 3D (?)

      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_misc, only: misc
      use mod_delaunay, only: identifyLargeTri, delaunaydata,
     &  CIRCUMSQFACHEX, getTriCenter, regenDelaunay
      use mod_materials, only: materials
      use mod_neighbors, only: neighbors
      implicit none
      
      type crackfindingdata
C     (read-in)
      integer :: mnum
C     (processed)
      real(dp) :: burgers
      type(delaunaydata) :: delaunay
      character(len=20) :: lattice
      end type
      
C     module variables
      type(crackfindingdata) :: crackfinding
      
      private
      public :: findCrack, initAtomFindCrackData, readAtomFindCrackData,
     & processAtomFindCrackData, writeAtomFindCrackData, crackfinding,
     & isTriangleOnEdge

C     HARD-CODED CONSTANTS     
      real(dp), parameter :: FUDGECIRCUM = 0.75_dp ! fudge for circumradius**2 cutoff
                                                   ! should be at least 0.5, and probably less than 1
                                                   ! if too small, defects  (such as dislocations) may be found
                                                   ! if too large, the position of the crack tip will be inaccurate (slightly behind the actual crack)

      contains
************************************************************************
      subroutine initAtomFindCrackData(findcrackfile)
 
C     Inputs: findcrackfile --- filename where crack finding information is stored
C     (should be something like '[filepref]_atomfindcrack')
 
C     Outputs: None
 
C     Purpose: Read, initialize data in "findcrack" structure, which holds
C     information about finding the crack in the atomistic region
           
      implicit none
      
C     input variables
      character(len=*) :: findcrackfile
      
      call readAtomFindCrackData(findcrackfile)
      call processAtomFindCrackData()
      
      end subroutine initAtomFindCrackData
************************************************************************
      subroutine readAtomFindCrackData(findcrackfile)
 
C     Inputs: findcrackfile --- filename where crack finding information is stored
C     (should be something like '[filepref]_atomfindcrack')
 
C     Outputs: None
 
C     Purpose: Read crack finding data from file
           
      implicit none
      
C     input variables
      character(len=*) :: findcrackfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=findcrackfile)
      
      read(iunit,*) crackfinding%mnum
      
      close(iunit)      
      
      end subroutine readAtomFindCrackData
************************************************************************
      subroutine processAtomFindCrackData()
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Set up crack finding burgers vector, lattice, delaunay triangulation, etc.
      
C     Notes/TODO: Assumes detection band and adjacent FE continuum material all has same lattice, burgers vectors

C     Notes/TODO: Needs to be reworked for 3D...

      implicit none
      
C     local variables
      integer :: mnum
      
      mnum = crackfinding%mnum
      crackfinding%burgers = materials(mnum)%burgers
      crackfinding%lattice = materials(mnum)%lattice
      
      if (crackfinding%lattice == 'hex') then
          crackfinding%delaunay%circumradiussqcutoff =
     &      FUDGECIRCUM*CIRCUMSQFACHEX*(crackfinding%burgers)**2/3.0_dp ! 1/sqrt(3) is factor for circumradius for equilateral triangle
          crackfinding%delaunay%nodenums = nodes%realatomlist
          allocate(crackfinding%delaunay%xy(2,nodes%nrealatoms))
          neighbors%delaunayregencrack = .true. ! need to generate triangulation
      end if
      
      end subroutine processAtomFindCrackData
************************************************************************  
      subroutine writeAtomFindCrackData(findcrackfile)
 
C     Inputs: findcrackfile --- filename where crack finding information is stored
C     (should be something like '[filepref]_atomfindcrack')
 
C     Outputs: None
 
C     Purpose: Write crack finding data to file (essentially
C     inverse of readAtomFindCrackData). Useful in creating "restart" file
           
      implicit none
      
C     input variables
      character(len=*) :: findcrackfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=findcrackfile)
      
      write(iunit,*) crackfinding%mnum
      
      close(iunit)      
      
      end subroutine writeAtomFindCrackData
************************************************************************      
      function findCrack() result(crackpos)
 
C     Inputs: None
 
C     Outputs: crackpos --- position of crack
 
C     Purpose: Find the position of the crack within the atomistic region.
C     Strategy is to triangulate the atoms, and determine the locations of
C     "large" triangles that are not at the edge of the region. The "rightmost"
C     large triangle corresponds to the crack tip.
      
      implicit none
      
C     output variables
      real(dp) :: crackpos(2)
      
C     local variables
      integer :: i
      real(dp) :: tricenter(2)
      
C     construct or update triangulation
      call regenDelaunay(crackfinding%delaunay,
     &                                     neighbors%delaunayregencrack)
     
      crackpos = [-huge(0.0_dp),0.0_dp]
C     loop over large triangles, excluding ones near edge
      do i = 1, size(crackfinding%delaunay%trigood)
      if (.not.(crackfinding%delaunay%trigood(i))) then
          if (.not.isTriangleOnEdge(crackfinding%delaunay,i)) then
              tricenter = getTriCenter(crackfinding%delaunay,i)
              if (tricenter(1) > crackpos(1)) then
                  crackpos = tricenter
              end if
          end if
      end if        
      end do
      
      end function findCrack
************************************************************************
      function isTriangleOnEdge(delaunay,trinum) result(isonedge)
 
C     Inputs: delaunay --- structure containing information about delaunay triangulation
C             trinum --- index of triangle in triangulation
 
C     Outputs: isonedge --- logical indicating whether triangle is on the edge of the triangulation
 
C     Purpose: Defected (large) triangles are either associated with the crack,
C     or with the outside edge of the triangulated region. This function
C     determines whether the triangle is on the edge by evaluating whether
C     any of its nodes are interface nodes
      
      implicit none
      
C     input variables
      type(delaunaydata) :: delaunay
      integer :: trinum
      
C     output variables
      logical :: isonedge
      
C     local variables
      integer :: i
      integer :: v, node, nodetype
      
      do i = 1, 3
          v = delaunay%trivert(i,trinum)
          node = delaunay%nodenums(v)
          nodetype = nodes%types(2,node)
          isonedge = (nodetype == 2)
          if (isonedge) then
              return
          end if
      end do
      
      end function isTriangleOnEdge
************************************************************************ 
      end module mod_find_crack_atomistic