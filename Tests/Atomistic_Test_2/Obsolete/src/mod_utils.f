      module mod_utils
      
C     Purpose: Contains auxiliary routines for pretty printing and writing
C     of two-dimensional arrays
      
      use mod_types, only: dp
      implicit none
      
      private
      public :: prettyPrintRealMat, writeRealMat,
     &          prettyPrintIntMat, writeIntMat
      
      contains
      
************************************************************************     
      subroutine prettyPrintRealMat(mat,matname)
      
C     Inputs: mat --- two-dimensional *real* array
C             matname --- name of array (character string)

C     Outputs: None
      
C     Purpose: Pretty-prints matrix (two-dimensional array) and its name
C     to the console.

      implicit none

C     input variables
      real(dp) :: mat(:,:)
      character(len=*) :: matname
      
C     local variables
      integer :: matsize(2)
      integer :: i,j
      
      matsize = shape(mat)
      write(*,*) matname
      do i = 1,matsize(1)
            write(*,*) (mat(i,j), j=1,matsize(2))
      end do
      end subroutine prettyPrintRealMat
************************************************************************     
      subroutine writeRealMat(mat,iunit,transposeoption)
      
C     Inputs: mat --- two-dimensional *real* array
C             iunit --- specifier for file (to be written to)

C     Outputs: None
      
C     Purpose: Writes matrix or its transpose (two-dimensional array) and its size
C     to a file.

      implicit none

C     input variables
      real(dp) :: mat(:,:)
      integer :: iunit
      logical :: transposeoption
      
C     local variables
      integer :: matsize(2)
      integer :: i,j
      
      matsize = shape(mat)
      
      if (transposeoption) then
          write(iunit,*) matsize(2), matsize(1)
          do i = 1,matsize(2)
                write(iunit,*) (mat(j,i), j=1,matsize(1))
          end do
      else
          write(iunit,*) matsize(1), matsize(2)
          do i = 1,matsize(1)
                write(iunit,*) (mat(i,j), j=1,matsize(2))
          end do
      end if    
      
      end subroutine writeRealMat
************************************************************************
      subroutine prettyPrintIntMat(mat,matname)
      
C     Inputs: mat --- two-dimensional *integer* array
C             matname --- name of array (character string)

C     Outputs: None
      
C     Purpose: Pretty-prints matrix (two-dimensional array) and its name
C     to the console

      implicit none

C     input variables
      integer :: mat(:,:)
      character(len=*) :: matname
      
C     local variables
      integer :: matsize(2)
      integer :: i,j
      
      matsize = shape(mat)
      write(*,*) matname
      do i = 1,matsize(1)
            write(*,*) (mat(i,j), j=1,matsize(2))
      end do
      end subroutine prettyPrintIntMat
************************************************************************
      subroutine writeIntMat(mat,iunit,transposeoption)

C     Inputs: mat --- two-dimensional *integer* array
C             iunit --- specifier for file (to be written to)

C     Outputs: None
      
C     Purpose: Writes matrix or its transpose (two-dimensional array) and its size
C     to a file.

      implicit none

C     input variables
      integer :: mat(:,:)
      integer :: iunit
      logical :: transposeoption
      
C     local variables
      integer :: matsize(2)
      integer :: i,j
      
      matsize = shape(mat)

      if (transposeoption) then
          write(iunit,*) matsize(2), matsize(1)
          do i = 1,matsize(2)
                write(iunit,*) (mat(j,i), j=1,matsize(1))
          end do
      else
          write(iunit,*) matsize(1), matsize(2)
          do i = 1,matsize(1)
                write(iunit,*) (mat(i,j), j=1,matsize(2))
          end do
      end if 

      end subroutine writeIntMat
************************************************************************ 
      
      end module