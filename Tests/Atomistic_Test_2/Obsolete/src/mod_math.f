      module mod_math
      
C     Purpose: Collection of utilities for "general purpose math"
C     (interpolation, distance finding, etc.) and for defining constants
C     (pi, tolerance, etc.).

C     Notes: Value of tol is used in various places,
C     typically to check before dividing by small numbers. 
C     Also, I have used the naming convention
C     of attaching const to the end of any parameters (e.g. pi -> piconst)
      
      use mod_types, only: dp
      implicit none
      
      private
      public :: tolconst, piconst, linearInterp,
     &          rotateVec2d, rotateStress2d, invertMat2,
     &          getUnitNormalRHR, getIntersectionTwoLines,
     &          getIdentityMatrix, getUniqueInts, getDuplicates,
     &          searchSortedBinary, searchSortedBrute,
     &          searchSortedSpecial, isMultiple
      
      real(dp), parameter :: tolconst = 1.0e-10_dp
      real(dp), parameter :: piconst = 3.14159265358979323846_dp
      
      contains

************************************************************************
      function linearInterp(x,xvec,yvec,n) result(y)
      
C     Inputs: x - value of independent variable
C             xvec - (float) vector of x-values of length n (assumed to be sorted)
C             yvec - (float) vector of y-values of length n
C             n - length of xvec, yvec

C     Outputs: value of y corresponding to x, using linear interpolation
C              over table (xvec, yvec). Returns error if x is outside
C              bounds of table
      
C     Purpose: Linear interpolation over table.

C     Notes: Does not check if x-values are sorted!
      
      implicit none
      
C     input variables
      real(dp) :: x
      real(dp) :: xvec(n), yvec(n)
      integer :: n
      
C     output variables
      real(dp) :: y
      
C     local variables
      integer :: idx
      real(dp) :: xlow, xhigh, ylow, yhigh
      
C     don't allow extrapolation
      if ((x < xvec(1)).or.(x > xvec(n))) then
          write(*,*) 'Bounds are bad: extrapolation'
          write(*,*) x
          stop
      end if
      idx = searchSortedBinary(x,xvec,n)
      
      xlow = xvec(idx)
      xhigh = xvec(idx+1)
      ylow = yvec(idx)
      yhigh = yvec(idx+1)
      y = (yhigh-ylow)/(xhigh-xlow)*(x-xlow) + ylow
      
      end function linearInterp
************************************************************************
      function searchSortedBrute(el,vec,n) result(idx)
      
C     Inputs: el - value for search
C             vec - vector of values, so that vec(1) <= el <= vec(n)
C             n - length of xvec, yvec

C     Outputs: index of el, where vec(index) <= el <= vec(index+1)

C     Purpose: Find index for table lookup (interpolation) using brute force search

C     Notes: Does not check if vec is sorted, or whether vec(1) <= el <= vec(n)!
      
      implicit none
      
C     input variables
      real(dp) :: el
      real(dp) :: vec(n)
      integer :: n
      
C     output variables
      integer :: idx
      
      idx = maxloc(vec,dim=1,mask=vec<=el)
      
      end function searchSortedBrute
************************************************************************
      function searchSortedBinary(el,vec,n) result(idxlow)

C     Inputs: el - value for search
C             vec - vector of values, so that vec(1) <= el <= vec(n)
C             n - length of xvec, yvec

C     Outputs: index of el, where vec(index) <= el <= vec(index+1)

C     Purpose: Find index for table lookup (interpolation) using binary search

C     Notes: Does not check if vec is sorted, or whether vec(1) <= el <= vec(n)!
      
      implicit none
      
C     input variables
      real(dp) :: el
      real(dp) :: vec(n)
      integer :: n
      
C     output variables
      integer :: idxlow
      
C     local variables
      integer :: idxhigh, idxmid
      
      idxlow = 1
      idxhigh = n
      
      do while (idxhigh - idxlow > 1)
          idxmid = (idxlow + idxhigh)/2
          if (el >= vec(idxmid)) then
              idxlow = idxmid
          else
              idxhigh = idxmid
          end if
      end do
      
      end function searchSortedBinary
************************************************************************
      function searchSortedSpecial(el,vec,n) result(idx)

C     Inputs: el - value for search
C             vec - vector of values, so that vec(1) <= el <= vec(n)
C             n - length of xvec, yvec

C     Outputs: index of el, where vec(index) <= el <= vec(index+1)

C     Purpose: Find index for table lookup (interpolation) assuming that
C     vec(2) ... vec(n-1) are evenly spaced (which is true for my tunable potentials)

C     Notes: Does not check if vec is sorted, or whether vec(1) <= el <= vec(n)!
      
      implicit none
      
C     input variables
      real(dp) :: el
      real(dp) :: vec(n)
      integer :: n
      
C     output variables
      integer :: idx
      
C     local variables
      integer :: idxlow, idxhigh
      real(dp) :: low, high, slope
      
      idxlow = 2
      low = vec(idxlow)
      idxhigh = n - 1
      high = vec(idxhigh)
      if (el <= low) then
          idx = 1
      else if (el >= high) then
          idx = n - 1
      else
          slope = (high - low)/(idxhigh - idxlow)
          idx = floor((el-low)/slope) + idxlow
      end if    
      
      end function searchSortedSpecial
************************************************************************
      subroutine rotateVec2d(cost,sint,xo,yo,xn,yn)
      
C     Inputs: cost, sint - cosine, sine of theta (rotation angle,
C                                                  counterclockwise)
C             xo, yo - coordinates of point in old coordinate system
      
C     Outputs: xn, yn - coordinates of point in new coordinate system

C     Purpose: Get new components of 2d vector, when coordinates are
C     rotated by an angle of theta, counterclockwise
      
      implicit none
      
C     input variables
      real(dp) :: cost, sint
      real(dp) :: xo, yo
      real(dp) :: xn, yn
      
      xn = xo*cost + yo*sint
      yn = -xo*sint + yo*cost
      
      end subroutine rotateVec2d
************************************************************************
      subroutine rotateStress2d(cost,sint,s11o,s22o,tauo,s11n,s22n,taun)

C     Inputs: cost, sint - cosine, sine of theta (rotation angle,
C                                                  counterclockwise)
C             s11o, s22o, tauo - old components of stress tensor
      
C     Outputs: s11n, s22n, taun - new components of stress tensor

C     Purpose: Get new components of 2d stress tensor, when coordinates are
C     rotated by an angle of theta, counterclockwise
      
      implicit none
      
C     input variables
      real(dp) :: cost, sint
      real(dp) :: cos2t, sin2t
      real(dp) :: s11o, s22o, tauo
      real(dp) :: s11n, s22n, taun
      real(dp) :: sigsum, sigdiffhalf
      
C     trig
      cos2t = cost**2 - sint**2
      sin2t = 2.0_dp*cost*sint
      
C     aux. constants      
      sigsum = s11o + s22o
      sigdiffhalf = 0.5_dp*(s11o - s22o)

C     rotation
      s11n = 0.5_dp*sigsum + sigdiffhalf*cos2t + tauo*sin2t
      taun =                -sigdiffhalf*sin2t + tauo*cos2t
      
C     use fact that trace of stress tensor is invariant
      s22n = sigsum - s11n
      
      end subroutine rotateStress2d
************************************************************************
      function invertMat2(mat) result(invmat)

C     Inputs: mat - real 2 by 2 array
      
C     Outputs: invmat - inverse of mat, real 2 by 2 array

C     Purpose: Get inverse of real 2 by 2 matrix, or return error if the
C     matrix is singular
      
      implicit none
      
C     input variables
      real(dp) :: mat(2,2)
      
C     output variables
      real(dp) :: invmat(2,2)
      
C     local variables
      real(dp) :: det, invdet
      
      det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
      
      if (abs(det) > tolconst) then
          invdet = 1.0_dp/det
          invmat(1,1) = invdet*mat(2,2)
          invmat(2,2) = invdet*mat(1,1)
          invmat(1,2) = -invdet*mat(1,2)
          invmat(2,1) = -invdet*mat(2,1)      
      else
          write(*,*) 'Matrix is singular'
          stop
      end if
      
      end function invertMat2
************************************************************************
      function getUnitNormalRHR(vec) result(normalvec)
      
C     Inputs: vec - real vector, length 2
      
C     Outputs: normalvec - unit normal to vec, real vector, length 2

C     Purpose: Get unit normal to a vector in x-y plane,
C     obeying the right hand rule (i.e. cross(vec/norm(vec),z) = normalvec)
      
      implicit none
      
C     input variables
      real(dp) :: vec(2)
      
C     output variables
      real(dp) :: normalvec(2)
      
C     local variables
      real(dp) :: fac

      normalvec(1) = vec(2)
      normalvec(2) = -vec(1)
      fac = 1.0_dp/sqrt(normalvec(1)**2 + normalvec(2)**2)
      normalvec = fac*normalvec
      
      end function getUnitNormalRHR
************************************************************************
      subroutine getIntersectionTwoLines(p0,p1,p2,p3,pint,isint)

C     Inputs: p0, p1, p2, p3 --- coordinates (vector, length 2) of points defining lines
C                                first line is from p0 to p1; second is from p2 to p3
      
C     Outputs: pint --- coordinates of intersection point (if it exists) (vector, length 2)
C              isint --- boolean indicating whether intersection exists

C     Purpose: Check whether two lines intersect, and, if so, where
      
C     Notes: Does not return correct results for all edge cases (e.g. coincident lines,
C     where, in fact, there are many intersections).
C     The only edge cases that should arise in the actual application
C     (detecting whether/where a dislocation has crossed the interface)
C     are: 1) a "T-Intersection" and 2) parallel, non-coincident lines.
C     I have checked these cases...
      
C     input variables
      real(dp) :: p0(2), p1(2)
      real(dp) :: p2(2), p3(2)
      
C     output variables
      real(dp) :: pint(2)
      logical :: isint
      
C     local variables
      real(dp) :: s10(2), s32(2), s02(2)
      real(dp) :: det, invdet
      real(dp) :: s, t
      
      s10 = p1 - p0
      s32 = p3 - p2
      s02 = p0 - p2

      det = s10(1)*s32(2) - s10(2)*s32(1)
      if (abs(det) < tolconst) then ! lines are essentially parallel
          isint = .false.
          return
      end if
      invdet = 1.0_dp/det
      s = (s10(1)*s02(2) - s10(2)*s02(1))*invdet
      t = (s32(1)*s02(2) - s32(2)*s02(1))*invdet

      isint = ((s >= 0).and.(s <= 1).and.(t >= 0).and.(t <= 1)) ! collision detected
      if (isint) then
          pint = p0 + t*s10
      end if
      
      end subroutine getIntersectionTwoLines
************************************************************************
      function getIdentityMatrix(N) result(idmat)

C     Function: getIdentityMatrix

C     Inputs: N --- integer, size of identity matrix

C     Outputs: idmat --- identity matrix of size N

C     Purpose: Generate identity matrix of size N
      
C     input variables
      integer :: N
      
C     output variables
      real(dp) :: idmat(N,N)
      
C     local variables
      integer :: i
      
      idmat = 0.0_dp
      do i = 1, N
          idmat(i,i) = 1.0_dp
      end do    
      
      end function getIdentityMatrix
************************************************************************
      subroutine getDuplicates(A,N,rowcheck,rowtot,Adup,Anondup,
     &                         Ndup,Nnondup,rowtotnew)

C     Subroutine: getDuplicates
      
C     Inputs: A --- integer array, rowtot by N, of sorted list of columns
C                   (e.g. generated by mergeSortCols)
C             N --- columns of A
C             rowcheck --- number of rows to check to determine duplicates
C                   (e.g. this equals rowsort in mergeSortcols)
C             rowtot --- rows of A

C     Outputs: Adup --- *allocatable* integer array of all the duplicates columns of A, rowtotnew by Ndup
C                       the array is larger because we may need information from the
C                       "extra rows" of A (i.e. rowcheck+1:rowtot for both entries)
C                       (see example below)
C              Anondup --- *allocatable* integer array of all the nonduplicate columns of A, rowtot by Nnondup
C              Ndup --- columns of Adup (number of duplicates)
C              Nnondup --- columns of Anondup (number of non-duplicates)
C              rowtotnew --- rows of Adup, equal to 2*rowtot - rowcheck
      
C     Purpose: Converts sorted list of columns into list with duplicates and
C     nonduplicates...useful for finding matching edges (interior)/non-matching edges (boundary)
C     for a mesh from the connectivity

C     Example: [[1,2,5],[1,3,5],[1,3,1],[4,1,1]] (A) -> [[1,2,5],[4,1,1]] (Anondup)
C              and [[1,3,5,1]] (Adup)

      implicit none
      
C     input variables
      integer :: rowtot
      integer :: rowcheck
      integer :: A(rowtot,N)
      integer :: N
      
C     output variables
      integer, allocatable :: Adup(:,:)
      integer, allocatable :: Anondup(:,:)
      integer :: Ndup
      integer :: Nnondup
      integer :: rowtotnew
      
C     local variables
      integer :: counter
      integer :: temp(rowtot), temp2(rowtot)
      logical :: match
      integer :: Aduptemp(2*rowtot - rowcheck,N)
      integer :: Anonduptemp(rowtot,N)
      
      Ndup = 0
      Nnondup = 0
      rowtotnew = 2*rowtot - rowcheck
      counter = 1
      do while (counter <= N)
          temp = A(:,counter)
          counter = counter + 1
          if (counter > N) then ! end of list
              match = .false.
          else
              temp2 = A(:,counter)
              match = all(temp(1:rowcheck)==temp2(1:rowcheck))
          end if
          if (match) then
C             we have a match; populate duplicate array
              Ndup = Ndup + 1
              Aduptemp(1:rowcheck,Ndup) = temp(1:rowcheck)
              Aduptemp(rowcheck+1:rowtot,Ndup) = temp(rowcheck+1:rowtot)
              Aduptemp(rowtot+1:rowtotnew,Ndup) =
     &                                          temp2(rowcheck+1:rowtot)
              counter = counter + 1
          else
C             we don't have a match; populate nonduplicate array
              Nnondup = Nnondup + 1
              Anonduptemp(:,Nnondup) = temp
          end if
      end do
      allocate(Adup(rowtotnew,Ndup))
      allocate(Anondup(rowtot,Nnondup))
      Adup = Aduptemp(:,1:Ndup)
      Anondup = Anonduptemp(:,1:Nnondup)
      
      end subroutine getDuplicates
************************************************************************
      subroutine getUniqueInts(mat,maxel,counter,intlist,invintlist)
 
C     Subroutine: getUniqueInts
 
C     Inputs: mat --- two-dimensional *integer* array
C             maxel --- largest element of mat

C     Outputs: counter --- number of unique elements in mat
C              intlist --- *allocatable*, sorted list, containing unique entries of 
C                          mat in locations 1:counter, and zeros thereafter
C                          dimensions maxel by 1
C              invintlist --- index list, where invintlist[intlist[i]] = i
C                             for i <= counter. Useful for getting node index
C                             from node number, for instance. Dimensions maxel by 1
C      
C     Purpose: Gets unique integers in an array, as well as index array

      implicit none
      
C     input variables
      integer :: mat(:,:)
      integer :: maxel
      
C     output variables
      integer :: counter
      integer, allocatable :: intlist(:)
      integer :: invintlist(maxel)
      
C     local variables
      integer :: matsize(2)
      integer :: i, j
      logical :: inlist(maxel)
      integer :: intlisttemp(maxel)
      
      matsize = shape(mat)
      inlist = .false.
      do i = 1, matsize(1)
          do j = 1, matsize(2)
              inlist(mat(i,j)) = .true.
          end do    
      end do
      
      counter = 0
      invintlist = 0
      do i = 1, maxel
          if (inlist(i)) then
              counter = counter + 1
              intlisttemp(counter) = i
              invintlist(i) = counter
          end if
      end do
      allocate(intlist(counter))
      intlist = intlisttemp(1:counter)
      
      end subroutine getUniqueInts
************************************************************************
      function isMultiple(x,n,def) result(res)
 
C     Subroutine: isMultiple
 
C     Inputs: x - integer value of interest
C             n - (integer) multiple
C             def - logical

C     Outputs: res - logical
C      
C     Purpose: Checks whether x is a multiple of n; returns def if n = 0
C              Useful for deciding whether to writeDump (based on dumpincrement),
C               etc.

      
C     input variables
      integer :: x
      integer :: n
      logical :: def
      
C     output variables
      logical :: res
      
      if (n/=0) then
          res = (mod(x,n) == 0)
      else
          res = def
      end if
      
      end function isMultiple
************************************************************************
      end module