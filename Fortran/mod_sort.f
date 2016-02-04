      module mod_sort
      
C     Purpose: Contains algorithms to sort integer lists or matrices
C     (by columns). Logic can be extended straightforwardly to arrays of other
C     types.
      
C     Possible extensions: Merge sorting might not be a good idea for
C     very large arrays; the large number of recursive calls might cause
C     memory issues. So, implement quicksort/timsort?
      
C     Notes: shamelessly stolen/adapted from
C     http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran

      use mod_types, only: dp
      
      implicit none
      
      private
      public :: mergeSort, mergeSub, mergeSortCols, mergeSubCols,
     &          insertionSortPlane, mergeSortReal, mergeSubReal,
     &          mergeSortColsReal, mergeSubColsReal
      
      contains
************************************************************************
      subroutine insertionSortPlane(relpos,objnum,N)

C     Subroutine: insertionSortPlane

C     Inputs: relpos --- unsorted list of floats (relative positions)
C             objnum --- object number corresponding to each relpos
C             N --- relpos, objnum

C     Outputs: A --- sorted list of integers

C     Purpose: Sort objects according to relative position within slip plane
C     (MATLAB equivalent: [relpos,idx] = sort(relpos); objnum = objnum(idx))

C     Notes: If we are dealing with small lists (N < ~100), insertionSort is quite fast
C     Furthermore, the relpos list, in most cases, will be mostly sorted already
C     In this case, insertionSort is very fast (O(N) for the case of already sorted)
      
      implicit none
      
C     input variables
      integer :: N
      
C     in/out variables
      real(dp) :: relpos(N)
      integer :: objnum(N)
      
C     local variables
      integer :: i, j
      real(dp) :: relpostemp
      integer :: objnumtemp
      
      do i = 2, N
          j = i - 1
          relpostemp = relpos(i)
          objnumtemp = objnum(i)
          do while (j >= 1) 
              if (relpos(j) <= relpostemp) then ! could combine with do while statement, but Fortran does not support short-circuiting
                  exit
              end if
              relpos(j+1) = relpos(j)
              objnum(j+1) = objnum(j)
              j = j - 1
          end do
          relpos(j+1) = relpostemp
          objnum(j+1) = objnumtemp
      end do
      
      end subroutine insertionSortPlane
************************************************************************
      recursive subroutine mergeSort(A,N)

C     Subroutine: mergeSort

C     Inputs: A --- (unsorted) list of integers
C             N --- length of A

C     Outputs: A --- sorted list of integers

C     Purpose: Sort list of integers using merge sort algorithm 

      implicit none
 
C     input variables
      integer :: N
      
C     in/out variables
      integer :: A(N)
      
C     local variables
      integer :: T(N)
      integer :: NA, NB
      integer :: V
 
      if (N == 1) then
          return
      else if (N == 2) then ! this branch is theoretically unnecessary,
C         but prevents empty vectors in mergeSub and unnecessary recursion
          if (A(1) > A(2)) then
              V = A(1)
              A(1) = A(2)
              A(2) = V
          end if
          return
      end if
      
C     integer division (implicit rounding)
      NA = (N+1)/2
      NB = N-NA
 
      call mergeSort(A(1:NA),NA)
      call mergeSort(A(NA+1:N),NB)
      call mergeSub(A(1:NA),NA,A(NA+1:N),NB,T,N)
      A = T
 
      end subroutine mergeSort
************************************************************************
      subroutine mergeSub(A,NA,B,NB,C,NC)
      
C     Subroutine: mergeSub

C     Inputs: A --- sorted list of integers
C             NA --- length of A
C             B --- sorted list of integers
C             NB --- length of B
C             NC --- = NA + NB, length of C

C     Outputs: C --- sorted list of integers, with elements from A and B

C     Purpose: Merge two sorted lists of integers, A and B, into a
C     sorted list of integers C. Helper routine for mergeSort
      
      implicit none
 
C     input variables
      integer :: NA, NB, NC
      integer :: A(NA)
      integer :: B(NB)
      
C     output variables
      integer :: C(NC)
 
C     local variables
      integer :: I, J, K
 
      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
          if (A(I) <= B(J)) then
              C(K) = A(I)
              I = I + 1
          else
              C(K) = B(J)
              J = J + 1
          endif
          K = K + 1
      end do
      do while (I <= NA)
          C(K) = A(I)
          I = I + 1
          K = K + 1
      enddo
      do while (J <= NB)
          C(K) = B(J)
          J = J + 1
          K = K + 1
      enddo
 
      end subroutine mergeSub  
************************************************************************
      recursive subroutine mergeSortReal(A,N)

C     Subroutine: mergeSortReal

C     Inputs: A --- (unsorted) list of real numbers
C             N --- length of A

C     Outputs: A --- sorted list of real numbers

C     Purpose: Sort list of real numbers using merge sort algorithm 

      implicit none
 
C     input variables
      integer :: N
      
C     in/out variables
      real(dp) :: A(N)
      
C     local variables
      real(dp) :: T(N)
      integer :: NA, NB
      real(dp) :: V
 
      if (N == 1) then
          return
      else if (N == 2) then ! this branch is theoretically unnecessary,
C         but prevents empty vectors in mergeSub and unnecessary recursion
          if (A(1) > A(2)) then
              V = A(1)
              A(1) = A(2)
              A(2) = V
          end if
          return
      end if
      
C     integer division (implicit rounding)
      NA = (N+1)/2
      NB = N-NA
 
      call mergeSortReal(A(1:NA),NA)
      call mergeSortReal(A(NA+1:N),NB)
      call mergeSubReal(A(1:NA),NA,A(NA+1:N),NB,T,N)
      A = T
 
      end subroutine mergeSortReal
************************************************************************
      subroutine mergeSubReal(A,NA,B,NB,C,NC)
      
C     Subroutine: mergeSubReal

C     Inputs: A --- sorted list of real numbers
C             NA --- length of A
C             B --- sorted list of real numbers
C             NB --- length of B
C             NC --- = NA + NB, length of C

C     Outputs: C --- sorted list of real numbers, with elements from A and B

C     Purpose: Merge two sorted lists of real numbers, A and B, into a
C     sorted list of real numbers C. Helper routine for mergeSortReal
      
      implicit none
 
C     input variables
      integer :: NA, NB, NC
      real(dp) :: A(NA)
      real(dp) :: B(NB)
      
C     output variables
      real(dp) :: C(NC)
 
C     local variables
      integer :: I, J, K
 
      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
          if (A(I) <= B(J)) then
              C(K) = A(I)
              I = I + 1
          else
              C(K) = B(J)
              J = J + 1
          endif
          K = K + 1
      end do
      do while (I <= NA)
          C(K) = A(I)
          I = I + 1
          K = K + 1
      enddo
      do while (J <= NB)
          C(K) = B(J)
          J = J + 1
          K = K + 1
      enddo
 
      end subroutine mergeSubReal 
************************************************************************
      recursive subroutine mergeSortCols(A,N,ROWSORT,ROWTOT)
 
C     Subroutine: mergeSortCols

C     Inputs: A --- (unsorted) array of integers, ROWTOT by N
C             N --- number of columns of A
C             ROWSORT --- number of rows to sort over to determine ordering
C                         e.g. [1,4] > [1,3] if we sort by two rows
C                         but not if we sort by one row
C             ROWTOT --- rows of A

C     Outputs: A --- sorted array of integers, according to first ROWSORT rows

C     Purpose: Sort array of integers according to first ROWSORT rows.
C     This is essentially equivalent to sortrows(A',1:ROWSORT)
C     in MATLAB (note the transpose)
 
      integer :: N
      integer :: ROWSORT, ROWTOT
      integer :: A(ROWTOT,N)
      
C     local variables
      integer :: T(ROWTOT,N)
      integer :: NA, NB
      logical :: test, flip
      integer :: L
      integer :: V(ROWTOT)
 
      if (N == 1) then
          return
      else if (N == 2) then  ! this branch is theoretically unnecessary,
C         but prevents empty matrices in mergeSubCols and unnecessary recursion
          test = .true.
          flip = .false.
          L = 1
          do while (test)
              if (A(L,1) < A(L,2)) then
                  test = .false.
              else if (A(L,1) > A(L,2)) then
                  test = .false.
                  flip = .true.
              else
                  if (L == ROWSORT) then
                      test = .false.
                  else
                      L = L + 1
                  end if    
              end if
          end do
          if (flip) then
              V = A(:,1)
              A(:,1) = A(:,2)
              A(:,2) = V
          end if
          return
      end if
      
      NA = (N+1)/2
      NB = N-NA
 
      call mergeSortCols(A(:,1:NA),NA,ROWSORT,ROWTOT)
      call mergeSortCols(A(:,NA+1:N),NB,ROWSORT,ROWTOT)
      call mergeSubCols(A(:,1:NA),NA,A(:,NA+1:N),NB,
     &                  T,N,ROWSORT,ROWTOT)
      A = T
 
      end subroutine mergeSortCols
************************************************************************
      subroutine mergeSubCols(A,NA,B,NB,C,NC,ROWSORT,ROWTOT)
 
C     Subroutine: mergeSubCols

C     Inputs: A --- sorted array of integers, ROWTOT by NA
C             NA --- columns of A
C             B --- sorted array of integers, ROWTOT by NB
C             NB --- columns of B
C             NC --- = NA + NB, columns of C
C             ROWSORT --- number of rows to sort over to determine ordering
C                         e.g. [1,4] > [1,3] if we sort by two rows
C                         but not if we sort by one row
C             ROWTOT --- rows of A, B, and C

C     Outputs: C --- sorted array of integers, with elements from A and B,
C                    ROWTOT by NC

C     Purpose: Merge two sorted arrays of integers, A and B, into a
C     sorted array of integers C. Helper routine for mergeSortCols.
 
C     input variables
      integer :: NA, NB, NC
      integer :: ROWSORT, ROWTOT
      integer :: A(ROWTOT,NA)
      integer :: B(ROWTOT,NB)
      
C     output variables
      integer :: C(ROWTOT,NC)
 
C     local variables
      integer :: I, J, K, L
      logical :: test
 
      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
          test = .true.
          L = 1
          do while (test)
              if (A(L,I) < B(L,J)) then
                  C(:,K) = A(:,I)
                  I = I + 1
                  test = .false.
              else if (A(L,I) > B(L,J)) then
                  C(:,K) = B(:,J)
                  J = J + 1
                  test = .false.
              else
                  if (L == ROWSORT) then
                      C(:,K) = A(:,I)
                      I = I + 1
                      test = .false.
                  else
                      L = L + 1
                  end if    
              end if    
          end do
          K = K + 1
      end do
      do while (I <= NA)
          C(:,K) = A(:,I)
          I = I + 1
          K = K + 1
      enddo
      do while (J <= NB)
          C(:,K) = B(:,J)
          J = J + 1
          K = K + 1
      enddo
 
      end subroutine mergeSubCols
************************************************************************
      recursive subroutine mergeSortColsReal(A,N,ROWSORT,ROWTOT)
 
C     Subroutine: mergeSortColsReal

C     Inputs: A --- (unsorted) array of real numbers, ROWTOT by N
C             N --- number of columns of A
C             ROWSORT --- number of rows to sort over to determine ordering
C                         e.g. [1,4] > [1,3] if we sort by two rows
C                         but not if we sort by one row
C             ROWTOT --- rows of A

C     Outputs: A --- sorted array of real numbers, according to first ROWSORT rows

C     Purpose: Sort array of real numbers according to first ROWSORT rows.
C     This is essentially equivalent to sortrows(A',1:ROWSORT)
C     in MATLAB (note the transpose)
 
      integer :: N
      integer :: ROWSORT, ROWTOT
      real(dp) :: A(ROWTOT,N)
      
C     local variables
      real(dp) :: T(ROWTOT,N)
      integer :: NA, NB
      logical :: test, flip
      integer :: L
      real(dp) :: V(ROWTOT)
 
      if (N == 1) then
          return
      else if (N == 2) then  ! this branch is theoretically unnecessary,
C         but prevents empty matrices in mergeSubCols and unnecessary recursion
          test = .true.
          flip = .false.
          L = 1
          do while (test)
              if (A(L,1) < A(L,2)) then
                  test = .false.
              else if (A(L,1) > A(L,2)) then
                  test = .false.
                  flip = .true.
              else
                  if (L == ROWSORT) then
                      test = .false.
                  else
                      L = L + 1
                  end if    
              end if
          end do
          if (flip) then
              V = A(:,1)
              A(:,1) = A(:,2)
              A(:,2) = V
          end if
          return
      end if
      
      NA = (N+1)/2
      NB = N-NA
 
      call mergeSortColsReal(A(:,1:NA),NA,ROWSORT,ROWTOT)
      call mergeSortColsReal(A(:,NA+1:N),NB,ROWSORT,ROWTOT)
      call mergeSubColsReal(A(:,1:NA),NA,A(:,NA+1:N),NB,
     &                  T,N,ROWSORT,ROWTOT)
      A = T
 
      end subroutine mergeSortColsReal
************************************************************************
      subroutine mergeSubColsReal(A,NA,B,NB,C,NC,ROWSORT,ROWTOT)
 
C     Subroutine: mergeSubColsReal

C     Inputs: A --- sorted array of real numbers, ROWTOT by NA
C             NA --- columns of A
C             B --- sorted array of real numbers, ROWTOT by NB
C             NB --- columns of B
C             NC --- = NA + NB, columns of C
C             ROWSORT --- number of rows to sort over to determine ordering
C                         e.g. [1,4] > [1,3] if we sort by two rows
C                         but not if we sort by one row
C             ROWTOT --- rows of A, B, and C

C     Outputs: C --- sorted array of real numbers, with elements from A and B,
C                    ROWTOT by NC

C     Purpose: Merge two sorted arrays of real numbers, A and B, into a
C     sorted array of real numbers C. Helper routine for mergeSortCols.
 
C     input variables
      integer :: NA, NB, NC
      integer :: ROWSORT, ROWTOT
      real(dp) :: A(ROWTOT,NA)
      real(dp) :: B(ROWTOT,NB)
      
C     output variables
      real(dp) :: C(ROWTOT,NC)
 
C     local variables
      integer :: I, J, K, L
      logical :: test
 
      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
          test = .true.
          L = 1
          do while (test)
              if (A(L,I) < B(L,J)) then
                  C(:,K) = A(:,I)
                  I = I + 1
                  test = .false.
              else if (A(L,I) > B(L,J)) then
                  C(:,K) = B(:,J)
                  J = J + 1
                  test = .false.
              else
                  if (L == ROWSORT) then
                      C(:,K) = A(:,I)
                      I = I + 1
                      test = .false.
                  else
                      L = L + 1
                  end if    
              end if    
          end do
          K = K + 1
      end do
      do while (I <= NA)
          C(:,K) = A(:,I)
          I = I + 1
          K = K + 1
      enddo
      do while (J <= NB)
          C(:,K) = B(:,J)
          J = J + 1
          K = K + 1
      enddo
 
      end subroutine mergeSubColsReal
************************************************************************
      end module