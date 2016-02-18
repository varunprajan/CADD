      module mod_disl_ident_simple

C     Purpose: Implements simple algorithm for dislocation identification
C     described in Stukowski, 2014, JMPS. Currently does not perform
C     extra correction step described in Section 4, although this is relatively easy
C     to implement and may be necessary for some applications (non-zero temperature, etc.).
      
C     Notes/TODO: Implement correction step!
      
      use mod_types, only: dp
      use mod_math, only: piconst
      use mod_delaunay, only: delaunay
      use mod_materials, only: materials
      
      implicit none
      
      type dislidentsimpledata
C     (processed)
      integer :: mnum
      character(len=20) :: lattice
      real(dp) :: burgers
      real(dp), allocatable :: ideal(:,:)
      end type
      
      type(dislidentsimpledata) :: identsimple
      
      private 
      public :: computeCircuits, initDislIdentData, identsimple,
     &          getIdealVector

      contains
************************************************************************
      subroutine initDislIdentData()
      
C     Subroutine: initDislIdentData

C     Inputs: None

C     Outputs: None
      
C     Purpose: Initialize data for "simple" identification of dislocations,
C     by generating ideal vectors for every possible combination
C     of slip system and dislocation sign (see Stukowski paper)

C     Notes/TODO: Assumes detection band and adjacent FE continuum material all has same lattice, burgers vectors

C     Notes/TODO: Maybe this should be merged with initDetectionData?  
      
      implicit none
      
C     local variables
      integer :: i, nvec
      integer :: mnum
      real(dp) :: theta
      real(dp) :: latticevec(2)

      mnum = identsimple%mnum
      identsimple%lattice = materials(mnum)%lattice
      identsimple%burgers = materials(mnum)%burgers
      if (identsimple%lattice == 'hex') then
          nvec = 6
          allocate(identsimple%ideal(2,nvec))
          do i = 1, nvec
              theta = (i - 1)*piconst/3.0_dp ! increments of 60 degrees
              latticevec = identsimple%burgers*[cos(theta),sin(theta)]
              identsimple%ideal(:,i) = latticevec
          end do
      end if    
      
      end subroutine initDislIdentData
************************************************************************
      function computeCircuits() result(circuits)
      
C     Function: computeCircuits

C     Inputs: None

C     Outputs: circuits --- array, 2 by numtri, containing burgers vector for each triangle
      
C     Purpose: Computes burgers vector using burgers circuit for each triangle in triangulation.
C     Uses symmetry to avoid duplicate computation (i.e. fact that Lba = -Lab).
      
C     Notes/TODO: This computes vectors for "bad" triangles, which is unnecessary,
C     except for those edges that neighbor good triangles.
C     
      implicit none
      
C     output variables
      real(dp), allocatable :: circuits(:,:)
      
C     local variables
      integer :: i, j, k
      integer :: nedges
      integer :: v1, v2, neighbor
      real(dp) :: p1(2), p2(2), L12(2)
      logical :: noneighbor
      
      allocate(circuits(2,delaunay%numtri))
      circuits = 0.0_dp      
      
      nedges = size(delaunay%trineigh,1)
      do i = 1, delaunay%numtri
          do j = 1, nedges
              k = mod(j,nedges) + 1 ! index of next node
              v1 = delaunay%trivert(j,i)
              v2 = delaunay%trivert(k,i)
              neighbor = delaunay%trineigh(j,i)
              noneighbor = (neighbor < 0)
              if ((v1 < v2).or.noneighbor) then ! use symmetry, except when there is no neighbor
                  p1 = delaunay%xy(:,v1)
                  p2 = delaunay%xy(:,v2)
                  L12 = getIdealVector(p1-p2) ! p1 - p2 because vertices are in ccw order,
                                              ! when they need to be in cw order (see Stukowski)
                  circuits(:,i) = circuits(:,i) + L12 
                  if (.not.noneighbor) then
                      circuits(:,neighbor) = circuits(:,neighbor) - L12
                  end if    
              end if
          end do    
      end do      
      
      end function computeCircuits
************************************************************************
      function getIdealVector(vec) result(idealvec)
      
C     Function: getIdealVector

C     Inputs: vec --- actual vector (length 2) connecting two atoms

C     Outputs: idealvec --- ideal vector (length 2)
      
C     Purpose: Computes ideal vector most closely corresponding to actual vector,
C     using 2-norm 

      implicit none
      
C     input variables
      real(dp) :: vec(2)
      
C     output variables
      real(dp) :: idealvec(2)
      
C     local variables
      real(dp) :: bestnorm, norm
      integer :: i
      real(dp) :: idealvectry(2)
      
      bestnorm = huge(0.0_dp)
      do i = 1, size(identsimple%ideal,2)
          idealvectry = identsimple%ideal(:,i)
          norm = sum((vec - idealvectry)**2)
          if (norm < bestnorm) then
              bestnorm = norm
              idealvec = idealvectry
          end if
      end do
      
      end function getIdealVector
************************************************************************
      end module mod_disl_ident_simple      
      