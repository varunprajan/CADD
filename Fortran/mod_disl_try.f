      module mod_disl_try
      
C     Purpose: Reads/writes/stores information about dislocations
C     obstacles, and sources for each fe (i.e. continuum) material.
C     Also contains methods for adding, deleting, updating dislocations,
C     etc.
      
      use mod_types, only: dp
      use mod_math, only: searchSortedBinary
      use mod_sort, only: insertionSortPlane
      use mod_mesh_find, only: findInAllInitiallyDef, 
     &        findInOneMatInitiallyDef, findInAllWithGuessDef
      use mod_slip_sys, only: getSlipPlane, slipsys
      use mod_materials, only: getNucleationLength
      use mod_disl_misc, only: dislmisc
      implicit none
      
      private
      public :: initDislData, writeDislData, disl, readDislData,
     &          addDislocation, sortPlaneCheck,
     &          deleteDislocation, checkTooManyDisl, dislt,
     &          initDislObsData, processDislObsData,
     &          writeDislObsData, initDislSourceData,
     &          processDislSourceData, writeDislSourceData,
     &          sources, obstacles, readDislObsData, readDislSourceData,
     &          initDislSortedPlanes, assignSourcesSortedPlanes,
     &          initObsSortedPlanes, assignDislSortedPlanes,
     &          assignObsSortedPlanes, sortPlanes, sortPlane,
     &          sortObsPlanes, sortDislPlanes,
     &          sortedplanedata, addDislocationSub, addObjSub,
     &          checkTooManyObj, processDislData, findEmptyDislSlot,
     &          initSourcesSortedPlanes, zeroDislDisp,
     &          setupSources, assignDislLocalPos, activateDislocations,
     &          assignSourcesLocalPos, assignNucleationLength,
     &          deleteDislocationSub, zeroObstacles, countActiveDisl,
     &          assignObsLocalPos, deleteDislocationSub2, sourcet,
     &          obstaclet

      type sortedplanedata
      integer :: nmax
      integer :: ncount
      integer, allocatable :: objnum(:)
      real(dp), allocatable :: relpos(:)
      logical :: resort
      end type

      type sortedplanesdata
      type(sortedplanedata), allocatable :: splane(:)
      end type

      type dislt
C     read-in
      integer :: cut
      real(dp) :: posn(2)
      integer :: sgn      
      integer :: slipsys
C     processed
      logical :: active
      real(dp) :: disp
      integer :: element
      real(dp) :: localpos(2)
      end type

      type sourcet
C     read-in
      real(dp) :: posn(2)
      integer :: slipsys
      real(dp) :: taucr
      real(dp) :: tnuc
C     processed
      integer :: element
      real(dp) :: localpos(2)
      real(dp) :: lnuc
      real(dp) :: tauprev
      real(dp) :: time
      end type
      
      type obstaclet
C     read-in
      real(dp) :: posn(2)
      integer :: slipsys
      real(dp) :: taucr
C     processed
      logical :: active
      logical :: computed
      integer :: element
      real(dp) :: localpos(2)
      end type
     
      type disldata
C     read-in (mostly)
      type(dislt), allocatable :: list(:)
      integer :: ndisl
C     processed
      type(sortedplanesdata), allocatable :: splanes(:)
      end type

      type obstacledata
C     read-in (mostly)
      type(obstaclet), allocatable :: list(:)
C     processed
      type(sortedplanesdata), allocatable :: splanes(:)
      end type
      
      type sourcedata
C     read-in (mostly)
      type(sourcet), allocatable :: list(:)
C     processed
      type(sortedplanesdata), allocatable :: splanes(:)
      end type
      
C     module variables
      type(disldata), allocatable :: disl(:)
      type(obstacledata), allocatable :: obstacles(:)
      type(sourcedata), allocatable :: sources(:)
      
      contains
************************************************************************
      subroutine initDislData(dislfile)

C     Inputs: dislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read, initialize data in "disl" structure, which holds
C     information about dislocation positions, orientations, etc. for
C     *each* continuum material

      implicit none
      
C     input variables
      character(len=*) :: dislfile
      
      call readDislData(dislfile)
      call processDislData()
      
      end subroutine initDislData
************************************************************************
      subroutine initDislSourceData(dislsourcefile)

C     Inputs: dislsourcefile --- filename where dislocation source data is stored
C     (should be something like '[filepref]_dislsource')

C     Outputs: None

C     Purpose: Read, initialize data in "sources" structure, which holds
C     information about dislocation sources

      implicit none
      
C     input variables
      character(len=*) :: dislsourcefile
      
      call readDislSourceData(dislsourcefile)
      call processDislSourceData()
      
      end subroutine initDislSourceData
************************************************************************
      subroutine initDislObsData(dislobsfile)

C     Inputs: dislobsfile --- filename where dislocation obstacle data is stored
C     (should be something like '[filepref]_dislobs')

C     Outputs: None

C     Purpose: Read, initialize data in "obstacles" structure, which holds
C     information about dislocation obstacles

      implicit none
      
C     input variables
      character(len=*) :: dislobsfile
      
      call readDislObsData(dislobsfile)
      call processDislObsData()
      
      end subroutine initDislObsData
************************************************************************
      subroutine readDislData(dislfile)

C     Inputs: dislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Read dislocation data (positions, etc.) from file,
C     for *each* continuum material, initialize/allocate structures/arrays

C     Notes: Dislocation arrays are allocated with # of rows equal to
C     nmaxdisl, which is a hard-coded constant. If ndisl > nmaxdisl,
C     an error is thrown. An alternative is to reallocate the array
C     every time ndisl > nmaxdisl...

      implicit none
      
C     input variables
      character(len=*) :: dislfile
      
C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: m, n
      integer :: nfematerials
      
      open(newunit=iunit,file=dislfile)
      
      read(iunit,*) nfematerials
      allocate(disl(nfematerials))
C     this should be identical to nfematerials from mod_fe_elements!
      do i = 1, nfematerials
          allocate(disl(i)%list(dislmisc%nmaxdisl(i)))

          read(iunit,*) m          
          do j = 1, m
              read(iunit,*) disl(i)%list(j)%cut
          end do
          
          read(iunit,*) m, n
          do j = 1, m
              read(iunit,*) (disl(i)%list(j)%posn(k), k = 1,n)
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) disl(i)%list(j)%sgn
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) disl(i)%list(j)%slipsys
          end do
          
          disl(i)%ndisl = m
      end do
      
      close(iunit)
      
      end subroutine readDislData
************************************************************************
      subroutine readDislSourceData(dislsourcefile)

C     Inputs: dislsourcefile --- filename where dislocation source data is stored
C     (should be something like '[filepref]_sources')

C     Outputs: None

C     Purpose: Read dislocation source data (positions, etc.) from file,
C     initialize/allocate structures/arrays
      
      implicit none
      
C     input variables
      character(len=*) :: dislsourcefile

C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: m, n
      integer :: nfematerials
      
      open(newunit=iunit,file=dislsourcefile)
      
      read(iunit,*) nfematerials
      allocate(sources(nfematerials))
C     this should be identical to nfematerials from mod_fe_elements!
      do i = 1, nfematerials
          read(iunit,*) m, n
          allocate(sources(i)%list(m))
          
          do j = 1, m
              read(iunit,*) (sources(i)%list(j)%posn(k), k = 1,n)
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) sources(i)%list(j)%slipsys
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) sources(i)%list(j)%taucr
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) sources(i)%list(j)%tnuc
          end do
      end do 
      
      close(iunit)
      
      end subroutine readDislSourceData
************************************************************************
      subroutine readDislObsData(dislobsfile)

C     Inputs: dislobsfile --- filename where dislocation obstacle data is stored
C     (should be something like '[filepref]_obstacles')

C     Outputs: None

C     Purpose: Read dislocation obstacle data (positions, etc.) from file,
C     initialize/allocate structures/arrays

      implicit none

C     input variables
      character(len=*) :: dislobsfile
      
C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: m, n
      integer :: nfematerials
      
      open(newunit=iunit,file=dislobsfile)
      
      read(iunit,*) nfematerials
      allocate(obstacles(nfematerials))
C     this should be identical to nfematerials from mod_fe_elements!
      do i = 1, nfematerials
          read(iunit,*) m, n
          allocate(obstacles(i)%list(m))
          
          do j = 1, m
              read(iunit,*) (obstacles(i)%list(j)%posn(k), k = 1,n)
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) obstacles(i)%list(j)%slipsys
          end do
          
          read(iunit,*) m
          do j = 1, m
              read(iunit,*) obstacles(i)%list(j)%taucr
          end do
      end do
      
      close(iunit)      
      
      end subroutine readDislObsData
************************************************************************
      subroutine writeDislData(dislfile)

C     Inputs: dislfile --- filename where dislocation data is stored
C     (should be something like '[filepref]_disl')

C     Outputs: None

C     Purpose: Write dislocation data to file (essentially
C     inverse of readDislData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: dislfile
      
C     local variables
      integer :: iunit
      integer :: i, j, k
      integer :: m, n, nfematerials
      
      open(newunit=iunit,file=dislfile)
      
      nfematerials = size(disl)
      write(iunit,*) nfematerials
      do i = 1, nfematerials
          m = countActiveDisl(i)
          write(iunit,*) m
          do j = 1, disl(i)%ndisl
              if (disl(i)%list(j)%active) then
                  write(iunit,*) disl(i)%list(j)%cut
              end if    
          end do
          write(iunit,*) ''
          
          n = size(disl(i)%list(1)%posn)
          write(iunit,*) m, n
          do j = 1, disl(i)%ndisl
              if (disl(i)%list(j)%active) then
                  write(iunit,*) (disl(i)%list(j)%posn(k), k = 1,n)
              end if    
          end do
          write(iunit,*) ''

          write(iunit,*) m          
          do j = 1, disl(i)%ndisl
              if (disl(i)%list(j)%active) then
                  write(iunit,*) disl(i)%list(j)%sgn
              end if    
          end do
          write(iunit,*) ''

          write(iunit,*) m
          do j = 1, disl(i)%ndisl
              if (disl(i)%list(j)%active) then
                  write(iunit,*) disl(i)%list(j)%slipsys
              end if    
          end do
          write(iunit,*) ''
          
      end do
      
      close(iunit)
      
      end subroutine writeDislData
************************************************************************
      function countActiveDisl(mnumfe) result(counter)

C     Inputs: mnumfe --- fe material number of dislocation(s)

C     Outputs: counter --- number of "active" dislocations within fe material
C     Some dislocations are "inactive", having escaped the mesh

C     Purpose: Count number of active dislocations within fe material   
      
      implicit none
      
C     input variables
      integer :: mnumfe
      
C     output variables
      integer :: counter
      
C     local variables
      integer :: i
      
      counter = 0
      do i = 1, disl(mnumfe)%ndisl
          if (disl(mnumfe)%list(i)%active) then
              counter = counter + 1
          end if    
      end do
      
      end function countActiveDisl
************************************************************************
      subroutine writeDislObsData(dislobsfile)

C     Inputs: dislobsfile --- filename where dislocation obstacle data is stored
C     (should be something like '[filepref]_dislobs')

C     Outputs: None

C     Purpose: Write dislocation obstacle data to file (essentially
C     inverse of readDislObsData). Useful in creating "restart" file
      
      implicit none

C     input variables
      character(len=*) :: dislobsfile
      
C     local variables
      integer :: iunit
      integer :: nfematerials
      integer :: m, n
      integer :: i, j, k
      
      open(newunit=iunit,file=dislobsfile)
      
      nfematerials = size(obstacles)
      write(iunit,*) nfematerials      
      do i = 1, nfematerials
          m = size(obstacles(i)%list)
          
          n = size(obstacles(i)%list(1)%posn)
          write(iunit,*) m, n
          do j = 1, m
              write(iunit,*) (obstacles(i)%list(j)%posn(k), k = 1,n)
          end do
          write(iunit,*) ''

          write(iunit,*) m          
          do j = 1, m
              write(iunit,*) obstacles(i)%list(j)%slipsys
          end do
          write(iunit,*) ''

          write(iunit,*) m
          do j = 1, m
              write(iunit,*) obstacles(i)%list(j)%taucr
          end do
          write(iunit,*) ''
          
      end do
      
      close(iunit)   
      
      end subroutine writeDislObsData
************************************************************************
      subroutine writeDislSourceData(dislsourcefile)

C     Inputs: dislsourcefile --- filename where dislocation source data is stored
C     (should be something like '[filepref]_dislsource')

C     Outputs: None

C     Purpose: Write dislocation source data to file (essentially
C     inverse of readDislSourceData). Useful in creating "restart" file
      
      implicit none

C     input variables
      character(len=*) :: dislsourcefile
      
C     local variables
      integer :: iunit
      integer :: nfematerials
      integer :: m, n
      integer :: i, j, k
      
      open(newunit=iunit,file=dislsourcefile)
      
      nfematerials = size(sources)
      write(iunit,*) nfematerials      
      do i = 1, nfematerials          
          m = size(sources(i)%list)
          
          n = size(sources(i)%list(1)%posn)
          write(iunit,*) m, n
          do j = 1, m
              write(iunit,*) (sources(i)%list(j)%posn(k), k = 1,n)
          end do
          write(iunit,*) ''

          write(iunit,*) m          
          do j = 1, m
              write(iunit,*) sources(i)%list(j)%slipsys
          end do
          write(iunit,*) ''

          write(iunit,*) m
          do j = 1, m
              write(iunit,*) sources(i)%list(j)%taucr
          end do
          write(iunit,*) ''
          
          write(iunit,*) m
          do j = 1, m
              write(iunit,*) sources(i)%list(j)%tnuc
          end do
          write(iunit,*) ''
          
      end do
      
      close(iunit)   
      
      end subroutine writeDislSourceData
************************************************************************ 
      subroutine processDislData()

C     Inputs: None

C     Outputs: None

C     Purpose: Assign dislocations to sorted planes structure, sort dislocations
C     according to relative position
      
      implicit none
      
      call activateDislocations()
      call zeroDislDisp()
      call assignDislLocalPos()
      call initDislSortedPlanes()
      call assignDislSortedPlanes()
      call sortDislPlanes()
      
      end subroutine processDislData
************************************************************************
      subroutine processDislSourceData()

C     Inputs: None

C     Outputs: None

C     Purpose: Assign sources to sorted planes structure, sort sources
C     according to relative position
      
      implicit none
      
      call setupSources()
      call assignNucleationLength()
      call assignSourcesLocalPos()
      call initSourcesSortedPlanes()
      call assignSourcesSortedPlanes()
      
      end subroutine processDislSourceData
************************************************************************
      subroutine processDislObsData()

C     Inputs: None

C     Outputs: None

C     Purpose: Assign obstacles to sorted planes structure, sort obstacles
C     according to relative position
      
      implicit none
      
      call zeroObstacles()
      call assignObsLocalPos()
      call initObsSortedPlanes()
      call assignObsSortedPlanes()
      call sortObsPlanes()
      
      end subroutine processDislObsData
************************************************************************
      subroutine activateDislocations()

C     Inputs: None

C     Outputs: None

C     Purpose: Activate all dislocations that have been read-in 
      
      implicit none
      
C     local variables
      integer :: i, j
      integer :: ndisl
      
      do i = 1, size(disl)
          ndisl = disl(i)%ndisl
          do j = 1, size(disl(i)%list)
              disl(i)%list(j)%active = (j <= ndisl)  
          end do
      end do    
      
      end subroutine activateDislocations
************************************************************************
      subroutine zeroDislDisp()

C     Inputs: None

C     Outputs: None

C     Purpose: Zero displacement vector for dislocations (i.e.
C     displacement along slip plane)
      
      implicit none
      
C     local variables
      integer :: i, j
      
      do i = 1, size(disl)
          do j = 1, disl(i)%ndisl
              disl(i)%list(j)%disp = 0.0_dp
          end do    
      end do
      
      end subroutine zeroDislDisp
************************************************************************
      subroutine setupSources()

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize timer and previous tau attributes for sources
      
      implicit none
      
C     local variables
      integer :: i, j
      
      do i = 1, size(sources)
          do j = 1, size(sources(i)%list)
              sources(i)%list(j)%time = 0.0_dp
              sources(i)%list(j)%tauprev = 0.0_dp
          end do    
      end do
      
      end subroutine setupSources
************************************************************************
      subroutine assignNucleationLength()

C     Inputs: None

C     Outputs: None

C     Purpose: Set nucleation length for each source

      implicit none

C     local variables
      integer:: i, j
      real(dp) :: taucr

      do i = 1, size(sources)
          do j = 1, size(sources(i)%list)
              taucr = sources(i)%list(j)%taucr
              sources(i)%list(j)%lnuc = getNucleationLength(i,taucr)
          end do
      end do  
      
      end subroutine assignNucleationLength
************************************************************************
      subroutine zeroObstacles()

C     Inputs: None

C     Outputs: None

C     Purpose: Zero "computed" attribute for obstacles (i.e. flag
C     indicating whether shear stress on obstacle has already been computed)
C     and "active" attribute (i.e. flag indicating whether obstacle is active --- i.e. if tau < taucr)   
      
      implicit none
      
C     local variables
      integer :: i, j
      
      do i = 1, size(obstacles)
          do j = 1, size(obstacles(i)%list)
              obstacles(i)%list(j)%computed = .false.
              obstacles(i)%list(j)%active = .false.
          end do
      end do    
      
      end subroutine zeroObstacles
************************************************************************
      subroutine assignDislLocalPos()

C     Inputs: None

C     Outputs: None

C     Purpose: Get FE element number, local position for each dislocation
C     using mesh_find routine
      
      implicit none

C     local variables
      integer :: i, j
      real(dp) :: x, y
      real(dp) :: r, s
      integer :: element
      logical :: badflip
      
      do i = 1, size(disl)
          do j = 1, disl(i)%ndisl
              x = disl(i)%list(j)%posn(1)
              y = disl(i)%list(j)%posn(2)
              call findInOneMatInitiallyDef(i,x,y,element,r,s,badflip)
              call checkImproperAssignment(badflip,x,y)
              disl(i)%list(j)%element = element
              disl(i)%list(j)%localpos(1) = r
              disl(i)%list(j)%localpos(2) = s
          end do    
      end do    
          
      end subroutine assignDislLocalPos   
************************************************************************
      subroutine assignSourcesLocalPos()

C     Inputs: None

C     Outputs: None

C     Purpose: Get FE element number, local position within element for each source
      
      implicit none

C     local variables
      integer :: i, j
      real(dp) :: x, y
      real(dp) :: r, s
      integer :: element
      logical :: badflip
      
      do i = 1, size(sources)
          do j = 1, size(sources(i)%list)
              x = sources(i)%list(j)%posn(1)
              y = sources(i)%list(j)%posn(2)
              call findInOneMatInitiallyDef(i,x,y,element,r,s,badflip)
              call checkImproperAssignment(badflip,x,y)
              sources(i)%list(j)%element = element
              sources(i)%list(j)%localpos = [r,s]
          end do    
      end do    
          
      end subroutine assignSourcesLocalPos
************************************************************************
      subroutine assignObsLocalPos()

C     Inputs: None

C     Outputs: None

C     Purpose: Get FE element number, local position within element for each obstacle
      
      implicit none

C     local variables
      integer :: i, j
      real(dp) :: x, y
      real(dp) :: r, s
      integer :: element
      logical :: badflip
      
      do i = 1, size(obstacles)
          do j = 1, size(obstacles(i)%list)
              x = obstacles(i)%list(j)%posn(1)
              y = obstacles(i)%list(j)%posn(2)
              call findInOneMatInitiallyDef(i,x,y,element,r,s,badflip)
              call checkImproperAssignment(badflip,x,y)
              obstacles(i)%list(j)%element = element
              obstacles(i)%list(j)%localpos = [r,s]
          end do    
      end do    
          
      end subroutine assignObsLocalPos
************************************************************************
      subroutine checkImproperAssignment(badflip,x,y)

C     Inputs: badflip --- logical indicating that object position
C                         within mesh could not be found
C             x, y --- position of object

C     Outputs: None

C     Purpose: If badflip, print helpful error message indicating where
C     object was when it could not be found within mesh
      
      implicit none
      
C     input variables
      logical :: badflip
      real(dp) :: x, y
      
      if (badflip) then
          write(*,*) 'Could not properly assign object'
          write(*,*) 'x', x, 'y', y
          stop
      end if
      
      end subroutine checkImproperAssignment
************************************************************************ 
      subroutine initDislSortedPlanes()

C     Inputs: None

C     Outputs: None

C     Purpose: Allocates sorted planes structures for disl, which contains
C     slip planes that hold dislocations in order of their relative position
C     along the slip plane
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, size(disl)
          call initSortedPlanes(i,disl(i)%splanes,
     &                          dislmisc%nmaxdislslip(i))
      end do
      
      end subroutine initDislSortedPlanes
************************************************************************
      subroutine initObsSortedPlanes()

C     Inputs: None

C     Outputs: None

C     Purpose: Allocates sorted planes structures for obstacles, which contains
C     slip planes that hold obstacles in order of their relative position
C     along the slip plane
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, size(obstacles)
          call initSortedPlanes(i,obstacles(i)%splanes,
     &                          dislmisc%nmaxobsslip(i))
      end do
      
      end subroutine initObsSortedPlanes
************************************************************************
      subroutine initSourcesSortedPlanes()

C     Inputs: None

C     Outputs: None

C     Purpose: Allocates sorted planes structures for sources, which contains
C     slip planes that hold sources (sources are unsorted...)
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, size(sources)
          call initSortedPlanes(i,sources(i)%splanes,
     &                          dislmisc%nmaxsrcslip(i))
      end do
      
      end subroutine initSourcesSortedPlanes
************************************************************************
      subroutine initSortedPlanes(mnumfe,splanes,nmaxobj)

C     Inputs: mnumfe --- fe material number of object
C             splanes --- sorted plane structure
C             nmaxobj --- number of objects allowed on a single slip plane

C     Outputs: None

C     Purpose: Allocates sorted planes structure (helper for initObjSortedPlanes,
C     initDislSortedPlanes)

      implicit none
      
C     input variables
      integer :: mnumfe
      type(sortedplanesdata), allocatable :: splanes(:)
      integer :: nmaxobj
      
C     local variables
      integer :: j, k
      integer :: nslipplanes, nslipsys
      
      nslipsys = size(slipsys(mnumfe)%theta)
      allocate(splanes(nslipsys))
      do j = 1, nslipsys
          nslipplanes = slipsys(mnumfe)%nslipplanes(j)
          allocate(splanes(j)%splane(nslipplanes))
          do k = 1, nslipplanes
              allocate(splanes(j)%splane(k)%objnum(nmaxobj))
              allocate(splanes(j)%splane(k)%relpos(nmaxobj))
              splanes(j)%splane(k)%ncount = 0
              splanes(j)%splane(k)%nmax = 0
              splanes(j)%splane(k)%objnum = 0
              splanes(j)%splane(k)%resort = .true.
          end do
      end do
      
      end subroutine initSortedPlanes
************************************************************************
      subroutine assignDislSortedPlanes()

C     Inputs: None

C     Outputs: None

C     Purpose: Assigns all dislocations to their respective slip planes, without sorting
      
      implicit none
      
C     local variables
      integer :: i, j
      integer :: isys, iplane
      real(dp) :: relpos
      
      do i = 1, size(disl)
          do j = 1, disl(i)%ndisl
              isys = disl(i)%list(j)%slipsys
              call getSlipPlane(disl(i)%list(j)%posn,
     &                          i,isys,iplane,relpos) ! adjusts pt to lie on slip plane
              call addObjSub(disl(i)%splanes(isys)%splane(iplane),
     &                        relpos,j)
          end do
      end do
      
      end subroutine assignDislSortedPlanes
************************************************************************
      subroutine assignObsSortedPlanes()
      
C     Inputs: None

C     Outputs: None

C     Purpose: Assigns all obstacles to their respective slip planes, without sorting
      
      implicit none
      
C     local variables
      integer :: i, j
      integer :: isys, iplane
      real(dp) :: relpos
      
      do i = 1, size(obstacles)
          do j = 1, size(obstacles(i)%list)
              isys = obstacles(i)%list(j)%slipsys
              call getSlipPlane(obstacles(i)%list(j)%posn,
     &                          i,isys,iplane,relpos) ! adjusts pt to lie on slip plane
              call addObjSub(obstacles(i)%splanes(isys)%splane(iplane),
     &                        relpos,j)
          end do
      end do
      
      end subroutine assignObsSortedPlanes
************************************************************************
      subroutine assignSourcesSortedPlanes()

C     Inputs: None

C     Outputs: None

C     Purpose: Assigns all sources to their respective slip planes, without sorting
            
      implicit none
      
C     local variables
      integer :: i, j
      integer :: isys, iplane
      real(dp) :: relpos
      
      do i = 1, size(sources)
          do j = 1, size(sources(i)%list)
              isys = sources(i)%list(j)%slipsys
              call getSlipPlane(sources(i)%list(j)%posn,
     &                          i,isys,iplane,relpos) ! adjusts pt to lie on slip plane
              call addObjSub(sources(i)%splanes(isys)%splane(iplane),
     &                        relpos,j)
          end do
      end do
      
      end subroutine assignSourcesSortedPlanes
************************************************************************
      subroutine sortDislPlanes()

C     Inputs: None

C     Outputs: None

C     Purpose: For each slip plane in disl, sorts dislocations according to
C     relative position along slip plane, keeping track of dislocation number (objnum)
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, size(disl)
          call sortPlanes(disl(i)%splanes)
      end do
      
      end subroutine sortDislPlanes
************************************************************************
      subroutine sortObsPlanes()

C     Purpose: For each slip plane in obstacles, sorts obstacles according to
C     relative position along slip plane, keeping track of obstacle number (objnum)
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, size(obstacles)
          call sortPlanes(obstacles(i)%splanes)
      end do
      
      end subroutine sortObsPlanes
************************************************************************
      subroutine sortPlanes(splanes)

C     In/Out: splanes --- sorted planes structure

C     Purpose: For each slip plane, sorts objects according to
C     relative position along slip plane, keeping track of object number (objnum)
      
      implicit none
      
C     in/out variables
      type(sortedplanesdata) :: splanes(:)
      
C     local variables
      integer :: j, k
      
      do j = 1, size(splanes)
          do k = 1, size(splanes(j)%splane)
              call sortPlaneCheck(splanes(j)%splane(k))
          end do
      end do
      
      end subroutine sortPlanes
************************************************************************
      subroutine sortPlaneCheck(splane)

C     In/Out: splane --- sorted plane structure

C     Purpose: If plane needs to be resorted, resort it using sortPlane
           
      implicit none
      
C     in/out variables
      type(sortedplanedata) :: splane
      
      if (splane%resort) then
          call sortPlane(splane)
          splane%resort = .false.
      end if
      
      end subroutine sortPlaneCheck
************************************************************************
      subroutine sortPlane(splane)

C     In/Out: splane --- sorted plane structure

C     Purpose: Sorts objects within slip plane according to
C     relative position, keeping track of object number (objnum)
           
      implicit none
      
C     in/out variables
      type(sortedplanedata) :: splane
      
C     local variables
      integer :: nmax
      
      nmax = splane%nmax
      call insertionSortPlane(splane%relpos(1:nmax),
     &                        splane%objnum(1:nmax),nmax)
      splane%nmax = splane%ncount
      
      end subroutine sortPlane
************************************************************************
      subroutine addDislocation(mnumfe,element,x,y,isys,bsgn,bcut)

C     Inputs: x, y --- global coordinates of dislocation
C             isys --- slip system for dislocation (in mnumfe)
C             bsgn --- sign of dislocation (+1 or -1)
C             bcut --- branch cut of dislocation (0 if to the left, 1 if to the right)

C     In/Out: mnumfe --- fe material to search in/material point was found in
C             element --- guess for starting element/element dislocation was found in
C                         (if zero, then no guess)

C     Outputs: dislnum --- number of dislocation in disl structure
      
C     Purpose: Adds dislocation to disl structure. First determines
C     empty slot for dislocation in disl structure. Then, assigns
C     dislocation its attributes. Finally, updates ndisl if necessary.
      
      implicit none
      
C     input variables
      real(dp) :: x, y
      integer :: isys
      integer :: bsgn, bcut
      
C     in/out variables
      integer :: mnumfe
      integer :: element      
      
C     local variables
      integer :: iplane
      real(dp) :: relpos
      integer :: dislnum
      logical :: badflip
      real(dp) :: r, s
      real(dp) :: posn(2)
      
      if (element == 0) then ! no/bad guess
          call findInAllInitiallyDef(x,y,mnumfe,element,r,s,badflip)
      else
          call findInAllWithGuessDef(x,y,mnumfe,element,r,s,badflip)
      end if
      
      call checkImproperAssignment(badflip,x,y)
      posn = [x,y]
      call getSlipPlane(posn,mnumfe,isys,iplane,relpos) ! adjusts posn slightly so that disl lies on slip plane
      call addDislocationSub(mnumfe,element,posn(1),posn(2),
     &                       isys,bsgn,bcut,r,s,dislnum)
      call addObjSub(disl(mnumfe)%splanes(isys)%splane(iplane),
     &               relpos,dislnum)
      
      end subroutine addDislocation
************************************************************************
      subroutine addDislocationSub(mnumfe,element,x,y,
     &                             isys,bsgn,bcut,r,s,dislnum)

C     Inputs: mnumfe --- fe material number of dislocation
C             element --- element number of dislocation
C             x, y --- global coordinates of dislocation
C             isys --- slip system for dislocation (in mnumfe)
C             bsgn --- sign of dislocation (+1 or -1)
C             bcut --- branch cut of dislocation (0 if to the left, 1 if to the right)
C             r, s --- local coordinates of dislocation in element

C     Outputs: dislnum --- number of dislocation in disl structure
      
C     Purpose: Adds dislocation to disl structure. First determines
C     empty slot for dislocation in disl structure. Then, assigns
C     dislocation its attributes. Finally, activates dislocation and 
C     updates ndisl if necessary.
      
      implicit none
      
C     input variables
      integer :: mnumfe, element
      real(dp) :: x, y
      integer :: isys
      integer :: bsgn, bcut
      real(dp) :: r, s
      
C     local variables
      integer :: dislnum
      
      call findEmptyDislSlot(mnumfe,dislnum)
      
      disl(mnumfe)%list(dislnum)%active = .true. ! activate disl
      disl(mnumfe)%list(dislnum)%cut = bcut
      disl(mnumfe)%list(dislnum)%posn(1) = x
      disl(mnumfe)%list(dislnum)%posn(2) = y
      disl(mnumfe)%list(dislnum)%slipsys = isys
      disl(mnumfe)%list(dislnum)%element = element
      disl(mnumfe)%list(dislnum)%localpos(1) = r
      disl(mnumfe)%list(dislnum)%localpos(2) = s
      disl(mnumfe)%list(dislnum)%sgn = bsgn
      disl(mnumfe)%list(dislnum)%disp = 0.0_dp
      
      if (dislnum > disl(mnumfe)%ndisl) then ! update dislnum
          disl(mnumfe)%ndisl = dislnum
      end if
      
      end subroutine addDislocationSub
************************************************************************
      subroutine findEmptyDislSlot(mnumfe,dislnum)

C     Inputs: mnumfe --- fe material number of dislocation 

C     Outputs: dislnum --- number of new dislocation in disl structure
      
C     Purpose: Finds number of empty slot for new dislocation in disl structure.
C     If no slot is available, return error
      
      implicit none
      
C     input variables
      integer:: mnumfe
      
C     output variables
      integer :: dislnum
      
C     local variables
      integer :: i
      integer :: nmax

      nmax = size(disl(mnumfe)%list)
      dislnum = nmax + 1
      do i = 1, nmax
          if (.not.disl(mnumfe)%list(i)%active) then
              dislnum = i
              return
          end if
      end do
      
      call checkTooManyDisl(dislnum,nmax,'disl')
      
      end subroutine findEmptyDislSlot
************************************************************************
      subroutine addObjSub(splane,relpos,objnum)

C     In/out: splane --- sorted plane structure

C     Inputs: relpos --- relative position of object along splane
C             objnum --- number of object in posn array
      
C     Purpose: Adds object to sorted plane structure
      
      implicit none
      
C     in/out variables
      type(sortedplanedata) :: splane
      
C     input variables
      real(dp) :: relpos
      integer :: objnum
      
C     local variables
      integer :: iobj
      
      iobj = splane%nmax + 1
      call checkTooManyObj(iobj,size(splane%relpos))
      splane%relpos(iobj) = relpos
      splane%objnum(iobj) = objnum
      splane%nmax = iobj
      splane%ncount = splane%ncount + 1
      splane%resort = .true.
      
      end subroutine addObjSub
************************************************************************
      subroutine deleteDislocation(mnumfe,isys,iplane,iobj)

C     Inputs: mnumfe --- material number that dislocation belongs to
C             isys --- number of slip system
C             iplane --- number of slip plane
C             iobj --- index of dislocation within slip plane

C     Outputs: None
      
C     Purpose: Remove dislocation from all disl structures

      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: isys, iplane, iobj
      
C     local variables
      integer :: dislnum
      
      dislnum = disl(mnumfe)%splanes(isys)%splane(iplane)%objnum(iobj)
      call deleteDislocationSub(mnumfe,dislnum)
      call deleteDislocationSub2(mnumfe,
     &                   disl(mnumfe)%splanes(isys)%splane(iplane),iobj)
      
      end subroutine deleteDislocation
************************************************************************
      subroutine deleteDislocationSub(mnumfe,dislnum)

C     Inputs: mnumfe --- material number that dislocation belongs to
C             dislnum --- number of dislocation

C     Outputs: None
      
C     Purpose: Remove dislocation from all disl structures

      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: dislnum
      
      disl(mnumfe)%list(dislnum)%active = .false. ! deactivate
      if (dislnum == disl(mnumfe)%ndisl) then ! decrement ndisl, if necessary
          disl(mnumfe)%ndisl = countActiveDisl(mnumfe)
      end if
      
      end subroutine deleteDislocationSub
************************************************************************
      subroutine deleteDislocationSub2(mnumfe,splane,iobj)

C     Inputs: mnumfe --- material number that dislocation belongs to
C             splane --- plane on which dislocation lies
C             iobj --- index of dislocation within slip plane

C     Outputs: None
      
C     Purpose: Remove object from sorted planes structure, decrement ncount, nmax
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: iobj
      
C     in/out variables
      type(sortedplanedata) :: splane
      
C     local variables
      integer :: i
      integer :: dislnum
      
      splane%relpos(iobj) = huge(0.0_dp) ! this ensures it will appear at end of list
C     don't delete objnum, since we might loop over object later...
      splane%ncount = splane%ncount - 1
      splane%resort = .true.
      
C     update nmax
      if (iobj == splane%nmax) then
          do i = splane%nmax, 1, -1
              dislnum = splane%objnum(i)
              if (disl(mnumfe)%list(dislnum)%active) then
                  splane%nmax = i
                  exit
              end if    
          end do
      end if
      
      end subroutine deleteDislocationSub2
************************************************************************
      subroutine checkTooManyDisl(dislnum,dislmax,disltype)

C     Inputs: dislnum --- number of dislocation
C             dislmax --- max. allowable number of total dislocations (i.e. columns of posn)

C     Outputs: None
      
C     Purpose: Checks to see if dislocation can be added to posn array. If not,
C     throws error and prints (hopefully) helpful message.

C     Notes: May want to replace this with re-allocation

      implicit none
      
C     input variables
      integer :: dislnum, dislmax
      character(len=*) :: disltype
      
      if (dislnum > dislmax) then
          write(*,*) 'Number of dislocations is too large'
          if (disltype == 'disl') then
              write(*,*) 'Increase nmaxdisl'
          else if (disltype == 'ghostdisl') then
              write(*,*) 'Increase nmaxghostdisl'
          else if (disltype == 'escapeddisl') then
              write(*,*) 'Increase nmaxescapeddisl'
          end if    
          stop          
      end if
      
      end subroutine checkTooManyDisl
************************************************************************
      subroutine checkTooManyObj(iobj,objmax)

C     Inputs: iobj --- index of object on slip plane
C             objmax --- max allowable index of object on slip plane (i.e. length of splane%relpos)

C     Outputs: None
      
C     Purpose: Checks to see if object can be added to a slip plane
C     If not, throws error and prints (hopefully) helpful message

C     Notes: May want to replace this with re-allocation

      implicit none
      
C     input variables
      integer :: iobj, objmax
      
      if (iobj > objmax) then
          write(*,*) 'Number of objects on a single plane'
          write(*,*) 'is too large.'
          stop          
      end if
      
      end subroutine checkTooManyObj
************************************************************************      
      end module mod_disl_try