      module mod_materials
      
C     Purpose: Stores information about materials, in an array of
C     structures (each structure contains information about a different
C     material).
      
C     Possible extensions: Add more fields to structure, to include
C     other material properties that are relevant (e.g. anisotropic
C     constants for dislocation fields --- p_i, etc.).
      
      use mod_types, only: dp
      use mod_utils, only: writeMatSize, readMatSize
      use mod_math, only: piconst
      implicit none
      
      private
      public :: materials, nmaterials,
     &          initMaterialData, writeMaterialData,
     &          readMaterialData, processMaterialData,
     &          getNucleationLength, getMuNuApprox
      
      type materialdata
C     read-in
      real(dp) :: burgers
      real(dp) :: disldrag
      real(dp) :: dislvmax
      real(dp), allocatable :: elconst(:,:)
      real(dp) :: lannih
      character(len=20) :: lattice
      real(dp) :: mass
      character(len=20) :: mname
      real(dp) :: rho
C     processed
      real(dp) :: rcoresq
      real(dp) :: mu
      real(dp) :: nu
      real(dp) :: dispprefac
      real(dp) :: stressprefac
      end type
      
C     module variables (global)
      type(materialdata), allocatable :: materials(:)
      integer :: nmaterials
      
C     module variables (private)
      real(dp), parameter :: corefac = 2.0_dp
      
      contains
************************************************************************      
      subroutine initMaterialData(materialfile)
      
C     Subroutine: initMaterialData

C     Inputs: materialfile --- filename where material data is stored
C     (should be something like '[filepref]_materials')

C     Outputs: None

C     Purpose: Read, initialize data in "materials" structure
      
      implicit none

C     input variables
      character(len=*) :: materialfile
      
      call readMaterialData(materialfile)
      call processMaterialData()
      
      end subroutine initMaterialData
************************************************************************      
      subroutine readMaterialData(materialfile)
      
C     Subroutine: readMaterialData

C     Inputs: materialfile --- filename where material data is stored
C     (should be something like '[filepref]_materials')

C     Outputs: None

C     Purpose: Read material data from file into "materials" structure
      
      implicit none

C     input variables
      character(len=*) :: materialfile
      
C     local variables
      integer :: iunit, i
      
      open(newunit=iunit,file=materialfile)
      
      read(iunit,*) nmaterials
      allocate(materials(nmaterials))
      do i = 1, nmaterials
C         read
          read(iunit,*) materials(i)%burgers
          read(iunit,*) materials(i)%disldrag
          read(iunit,*) materials(i)%dislvmax
          call readMatSize(iunit,materials(i)%elconst)
          read(iunit,*) materials(i)%lannih
          read(iunit,*) materials(i)%lattice
          read(iunit,*) materials(i)%mass
          read(iunit,*) materials(i)%mname
          read(iunit,*) materials(i)%rho
      end do
      
      close(iunit)
      
      end subroutine readMaterialData
************************************************************************
      subroutine processMaterialData()
      
C     Subroutine: processMaterialData

C     Inputs: None

C     Outputs: None

C     Purpose: Computes auxiliary quantities, e.g. mu, nu, using material
C     data

C     Notes: Requires material data to already be read, using readMaterialData
      
      implicit none
      
C     local variables
      real(dp) :: mu, nu
      integer :: i
      
      do i = 1, nmaterials
          call getMuNuApprox(materials(i)%elconst,mu,nu)
          materials(i)%mu = mu
          materials(i)%nu = nu
          materials(i)%rcoresq = (corefac*materials(i)%burgers)**2
          materials(i)%dispprefac =
     &               0.25_dp*materials(i)%burgers/(piconst*(1.0_dp-nu))
          materials(i)%stressprefac = 
     &             0.5_dp*materials(i)%burgers*mu/(piconst*(1.0_dp-nu))
      end do
      
      end subroutine processMaterialData
************************************************************************
      subroutine writeMaterialData(materialfile)

C     Subroutine: writeMaterialData

C     Inputs: materialfile --- filename where material data is stored
C     (should be something like '[filepref]_materials')

C     Outputs: None

C     Purpose: Write material data to file (essentially
C     inverse of readMaterialData). Useful in creating "restart" file
C     (i.e. resulting file should be able to be read in by readMaterialData)
      
      implicit none

C     input variables
      character(len=*) :: materialfile
      
C     local variables
      integer :: iunit, i
      
      open(newunit=iunit,file=materialfile)
      
      write(iunit,*) nmaterials
      do i = 1, nmaterials
          write(iunit,*) materials(i)%burgers
          write(iunit,*) materials(i)%disldrag
          write(iunit,*) materials(i)%dislvmax
          call writeMatSize(iunit,materials(i)%elconst)
          write(iunit,*) materials(i)%lannih
          write(iunit,*) materials(i)%lattice
          write(iunit,*) materials(i)%mass
          write(iunit,*) materials(i)%mname
          write(iunit,*) materials(i)%rho
          write(iunit,*) ''
      end do
      
      close(iunit)
      
      end subroutine writeMaterialData
************************************************************************
      subroutine getMuNuApprox(elconst,mu,nu)
      
C     Subroutine: getMuNuApprox

C     Inputs: elconst - 3 by 3 plane-strain stiffness matrix (Voigt notation)

C     Outputs: mu - shear modulus, assuming material is isotropic
C              nu - Poisson's ratio, assuming material is isotropic

C     Purpose: Computes mu, nu using stiffness matrix
      
C     Notes: Since elconst matrix is in Voigt notation,
C     elconst(3,3) = mu, not 2*mu. This subroutine does not work well
C     if the material is not isotropic, but in that case the anisotropic
C     expressions for dislocations fields, etc. should be used instead
C     of the isotropic ones that rely on mu, nu
      
      implicit none
      
C     input variables
      real(dp) :: elconst(3,3)
      
C     output variables
      real(dp) :: mu, nu
      
C     local variables
      
      mu = elconst(3,3)
      nu = elconst(1,2)/(elconst(1,1)+elconst(1,2))
      
      end subroutine getMuNuApprox
************************************************************************
      function getNucleationLength(mnumfe,taucr) result(lnuc)
      
C     Function: getNucleationLength

C     Inputs: mnumfe --- fe material number of source
C             taucr --- critical shear stress for source

C     Outputs: lnuc --- separation length between the dislocations
      
C     Notes: Needs to be changed for anisotropic materials! For isotropic
C     formula, see Segurado et al., 2007, MSMSE, Eqn. 1

      implicit none

C     input variables
      integer :: mnumfe
      real(dp) :: taucr
      
C     output variables
      real(dp) :: lnuc

C     local variables
      real(dp) :: mu, nu, burgers

      mu = materials(mnumfe)%mu
      nu = materials(mnumfe)%nu
      burgers = materials(mnumfe)%burgers
      lnuc = mu*burgers/(2.0_dp*piconst*(1.0_dp-nu)*taucr)
      
      end function getNucleationLength
************************************************************************      
      
      end module