      module mod_slip_sys

      use mod_types, only: dp
      use mod_math, only: projectVec, TOLCONST
      use mod_utils, only: writeMatTransposeSize, readMatTransposeSize,
     &                     writeVecSize, readVecSize
      implicit none
      
      private
      public :: initSlipSysData, readSlipSysData, processSlipSysData,
     &          getSlipPlane, writeSlipSysData, slipsys, invResolveDisp,
     &          resolveStress

      type slipsystemdata
C     read-in
      integer, allocatable :: nslipplanes(:)
      real(dp), allocatable :: origin(:,:)
      real(dp), allocatable :: space(:)
      real(dp), allocatable :: theta(:)
C     processed
      real(dp), allocatable :: trig(:,:)
      end type

      type(slipsystemdata), allocatable :: slipsys(:)

      contains
************************************************************************
      subroutine initSlipSysData(slipsysfile)

C     Inputs: slipsysfile --- filename where slip system data is stored
C     (should be something like '[filepref]_slipsys')

C     Outputs: None

C     Purpose: Read, initialize data in "slipsys" structure, which holds
C     information about slip systems and their orientations

      implicit none
      
C     input variables
      character(len=*) :: slipsysfile
      
      call readSlipSysData(slipsysfile)
      call processSlipSysData()
      
      end subroutine initSlipSysData
************************************************************************
      subroutine readSlipSysData(slipsysfile)

C     Inputs: slipsysfile --- filename where slip system data is stored
C     (should be something like '[filepref]_slipsys')

C     Outputs: None

C     Purpose: Read slip system data (angles of slip system, etc.) from file,
C     for *each* continuum material, initialize/allocate structures/arrays

      implicit none
      
C     input variables
      character(len=*) :: slipsysfile
      
C     local variables
      integer :: iunit
      integer :: i
      integer :: nfematerials

      open(newunit=iunit,file=slipsysfile)      

      read(iunit,*) nfematerials
      allocate(slipsys(nfematerials))

      do i = 1, nfematerials
          call readVecSize(iunit,slipsys(i)%nslipplanes)
          call readMatTransposeSize(iunit,slipsys(i)%origin)
          call readVecSize(iunit,slipsys(i)%space)
          call readVecSize(iunit,slipsys(i)%theta)
      end do
      
      close(iunit)
      
      end subroutine readSlipSysData
************************************************************************
      subroutine processSlipSysData()

C     Inputs: None

C     Outputs: None

C     Purpose: Computes and stores sin(theta), cos(theta) for each slip system,
C     to avoid recomputing each time tilde field is needed
      
      implicit none
      
C     local variables
      integer :: i, j
      integer :: nslipsys
      real(dp) :: theta
      
      do i = 1, size(slipsys)
          nslipsys = size(slipsys(i)%theta)
          allocate(slipsys(i)%trig(2,nslipsys)) ! cos and sin
          do j = 1, nslipsys
              theta = slipsys(i)%theta(j)
              slipsys(i)%trig(1,j) = cos(theta)
              slipsys(i)%trig(2,j) = sin(theta)
          end do
      end do
      
      end subroutine processSlipSysData
************************************************************************
      subroutine getSlipPlane(pt,mnumfe,isys,iplane,relpos)

C     Inputs: mnumfe --- material number that object belongs to
C             isys --- index of slip system within mnumfe that object belongs to

C     In/out: pt --- vector, length 2, object position. Position is adjusted
C             by subroutine to lie exactly on slip plane

C     Outputs: iplane --- index of slip plane that object belongs to (i.e. is closest to)
C              relpos --- relative position of object along that slip plane (relative to origin)

C     Purpose: Finds closest slip plane for an object at position pt, and its relative
C              position along that slip plane (with respect to origin)
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: isys
      
C     in/out variables
      real(dp) :: pt(2)
      
C     output variables
      integer :: iplane
      real(dp) :: relpos
      
C     local variables
      real(dp) :: origin(2), linevec(2), ptvec(2)
      real(dp), allocatable :: ptvecproj(:)
      real(dp) :: ptvecrej(2)
      real(dp) :: space
      real(dp) :: pdist
      real(dp) :: posperp(2), posparallel(2)
      
      origin = slipsys(mnumfe)%origin(:,isys)
      space = slipsys(mnumfe)%space(isys)
      linevec = slipsys(mnumfe)%trig(:,isys)
      ptvec = pt - origin
      ptvecproj = projectVec(ptvec,linevec)
      ptvecrej = ptvec - ptvecproj
      pdist = sqrt(sum(ptvecrej**2))
      iplane = nint(pdist/space) + 1 ! first slip plane is 1, not zero
      relpos = dot_product(ptvecproj,linevec)
      if (pdist < TOLCONST) then
          posperp = [0.0_dp,0.0_dp]
      else    
          posperp = ((iplane-1)*space)*(ptvecrej/pdist)
      end if    
      posparallel = relpos*linevec
      pt = origin + posperp + posparallel
      
      end subroutine getSlipPlane
************************************************************************
      function invResolveDisp(mnumfe,isys,relposdisp) result(disp)

C     Inputs: mnumfe --- material number that object belongs to
C             isys --- index of slip system within mnumfe that object belongs to
C             relposdisp --- displacement of object along slip plane

C     Outputs: disp --- displacement of object in global coordinate system (vector, 2 by 1)

C     Purpose: Goes from local displacement along slip plane to global displacement
C              ("inverse" of resolving global displacement into local coords)    
      
C     input variables
      integer :: mnumfe
      integer :: isys
      real(dp) :: relposdisp
      
C     output variables
      real(dp) :: disp(2)
      
      disp = relposdisp*slipsys(mnumfe)%trig(:,isys)
      
      end function invResolveDisp
************************************************************************
      function resolveStress(mnumfe,isys,stress) result(tau)

C     Inputs: mnumfe --- material number that object belongs to
C             isys --- index of slip system within mnumfe that object belongs to
C             stress --- stress in global coordinate system (vector, 3 by 1)

C     Outputs: tau --- shear stress along slip plane

C     Purpose: Resolves global stress along slip plane to determine shear stress 
C     (see Segurado and Llorca, 2007, Eqn. 3)
      
C     input variables
      integer :: mnumfe
      integer :: isys
      real(dp) :: stress(3)
      
C     output variables
      real(dp) :: tau
      
C     local variables
      real(dp) :: cost, sint
      real(dp) :: tx, ty
      
      cost = slipsys(mnumfe)%trig(1,isys)
      sint = slipsys(mnumfe)%trig(2,isys)
      tx = stress(1)*cost + stress(3)*sint
      ty = stress(3)*cost + stress(2)*sint
      tau = ty*cost - tx*sint
      
      end function resolveStress
************************************************************************
      subroutine writeSlipSysData(slipsysfile)

C     Inputs: slipsysfile --- filename where slip system data is stored
C     (should be something like '[filepref]_slipsys')

C     Outputs: None

C     Purpose: Write slip system data to file (essentially
C     inverse of readSlipSysData). Useful in creating "restart" file

      implicit none
      
C     input variables
      character(len=*) :: slipsysfile
      
C     local variables
      integer :: iunit
      integer :: i
      integer :: nfematerials
      
      open(newunit=iunit,file=slipsysfile)
      
      nfematerials = size(slipsys)
      write(iunit,*) nfematerials
      do i = 1, nfematerials
          call writeVecSize(iunit,slipsys(i)%nslipplanes)
          call writeMatTransposeSize(iunit,slipsys(i)%origin)
          call writeVecSize(iunit,slipsys(i)%space)
          call writeVecSize(iunit,slipsys(i)%theta)
      end do
      
      close(iunit)
      
      end subroutine writeSlipSysData
************************************************************************
      end module