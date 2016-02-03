      module mod_utils
      
C     Purpose: Contains auxiliary routines for pretty printing, reading, and writing
C     of one- and two-dimensional arrays
      
      use mod_types, only: dp
      implicit none
      
      private :: prettyPrintRealMat, writeRealMat,
     &          prettyPrintIntMat, writeIntMat
      public ::  prettyPrintVec, prettyPrintMat, readMatSize,
     &          readMatTransposeSize, readVecSize, writeMatSize,
     &          writeMatTransposeSize, writeMat, writeMatTranspose,
     &          writeVecSize, writeVec
      
************************************************************************
      interface prettyPrintVec
          module procedure prettyPrintRealVec, prettyPrintIntVec   
      end interface
************************************************************************       
      interface prettyPrintMat
          module procedure prettyPrintRealMat, prettyPrintIntMat    
      end interface
************************************************************************
      interface readMatSize
          module procedure readRealMatSize, readIntMatSize
      end interface
************************************************************************
      interface readMatTransposeSize
          module procedure readRealMatTransposeSize,
     &                     readIntMatTransposeSize
      end interface
************************************************************************
      interface readVecSize
          module procedure readRealVecSize, readIntVecSize
      end interface    
************************************************************************      
      interface writeMatSize
          module procedure writeRealMatSize, writeIntMatSize
      end interface
************************************************************************
      interface writeMatTransposeSize
          module procedure writeRealMatTransposeSize,
     &                     writeIntMatTransposeSize
      end interface
************************************************************************
      interface writeMat
          module procedure writeRealMat, writeIntMat
      end interface
************************************************************************
      interface writeMatTranspose
          module procedure writeRealMatTranspose, writeIntMatTranspose
      end interface
************************************************************************
      interface writeVecSize
          module procedure writeRealVecSize, writeIntVecSize
      end interface    
************************************************************************
      interface writeVec
          module procedure writeRealVec, writeIntVec
      end interface    
************************************************************************
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
      integer :: i, j
      
      matsize = shape(mat)
      write(*,*) matname
      do i = 1, matsize(1)
            write(*,*) (mat(i,j), j=1,matsize(2))
      end do
      
      end subroutine prettyPrintRealMat
************************************************************************
      subroutine prettyPrintIntMat(mat,matname)
      
C     Inputs: mat --- two-dimensional *integer* array
C             matname --- name of array (character string)

C     Outputs: None
      
C     Purpose: Pretty-prints matrix (two-dimensional array) and its name

      implicit none

C     input variables
      integer :: mat(:,:)
      character(len=*) :: matname
      
C     local variables
      integer :: matsize(2)
      integer :: i, j
      
      matsize = shape(mat)
      write(*,*) matname
      do i = 1, matsize(1)
            write(*,*) (mat(i,j), j=1,matsize(2))
      end do
      
      end subroutine prettyPrintIntMat
************************************************************************
      subroutine prettyPrintRealVec(vec,vecname)
      
C     Inputs: vec --- one-dimensional *real* array
C             vecname --- name of array (character string)

C     Outputs: None
      
C     Purpose: Pretty-prints vector *as a column* and its name

      implicit none

C     input variables
      real(dp) :: vec(:)
      character(len=*) :: vecname
      
C     local variables
      integer :: i
      
      write(*,*) vecname
      do i = 1, size(vec)
            write(*,*) vec(i)
      end do
      
      end subroutine prettyPrintRealVec
************************************************************************
      subroutine prettyPrintIntVec(vec,vecname)
      
C     Inputs: vec --- one-dimensional *integer* array
C             vecname --- name of array (character string)

C     Outputs: None
      
C     Purpose: Pretty-prints vector *as a column* and its name

      implicit none

C     input variables
      integer :: vec(:)
      character(len=*) :: vecname
      
C     local variables
      integer :: i
      
      write(*,*) vecname
      do i = 1, size(vec)
            write(*,*) vec(i)
      end do
      
      end subroutine prettyPrintIntVec
************************************************************************
      subroutine readRealVecSize(iunit,vec)
      
C     Inputs: iunit --- specifier for file (to be written to)

C     Outputs: vec --- one-dimensional *real* array
      
C     Purpose: Reads vector from file, using size.

      implicit none

C     input variables
      integer :: iunit
      
C     output variables
      real(dp), allocatable :: vec(:)
      
C     local variables
      integer :: i
      integer :: n
      
      read(iunit,*) n
      allocate(vec(n))
      read(iunit,*) (vec(i),i=1,n)
      
      end subroutine readRealVecSize
************************************************************************
      subroutine readIntVecSize(iunit,vec)
      
C     Inputs: iunit --- specifier for file (to be written to)

C     Outputs: vec --- one-dimensional *integer* array
      
C     Purpose: Reads vector from file, using size.

      implicit none

C     input variables
      integer :: iunit
      
C     output variables
      integer, allocatable :: vec(:)
      
C     local variables
      integer :: i
      integer :: n
      
      read(iunit,*) n
      allocate(vec(n))
      read(iunit,*) (vec(i),i=1,n)
      
      end subroutine readIntVecSize
************************************************************************
      subroutine readRealMatSize(iunit,mat)
      
C     Inputs: iunit --- specifier for file (to be written to)

C     Outputs: mat --- two-dimensional *real* array
      
C     Purpose: Reads matrix from file, using size.

      implicit none

C     input variables
      integer :: iunit
      
C     output variables
      real(dp), allocatable :: mat(:,:)
      
C     local variables
      integer :: j, k
      integer :: m, n
      
      read(iunit,*) m, n
      allocate(mat(m,n))
      read(iunit,*) ((mat(j,k),k=1,n),j=1,m)
      
      end subroutine readRealMatSize
************************************************************************
      subroutine readIntMatSize(iunit,mat)
      
C     Inputs: iunit --- specifier for file (to be written to)

C     Outputs: mat --- two-dimensional *integer* array
      
C     Purpose: Reads matrix from file, using size.

      implicit none

C     input variables
      integer :: iunit
      
C     output variables
      integer, allocatable :: mat(:,:)
      
C     local variables
      integer :: j, k
      integer :: m, n
      
      read(iunit,*) m, n
      allocate(mat(m,n))
      read(iunit,*) ((mat(j,k),k=1,n),j=1,m)
      
      end subroutine readIntMatSize
************************************************************************
      subroutine readRealMatTransposeSize(iunit,mat)
      
C     Inputs: iunit --- specifier for file (to be written to)

C     Outputs: mat --- two-dimensional *real* array
      
C     Purpose: Reads matrix from file, using size, storing it in transposed form.

      implicit none

C     input variables
      integer :: iunit
      
C     output variables
      real(dp), allocatable :: mat(:,:)
      
C     local variables
      integer :: j, k
      integer :: m, n
      
      read(iunit,*) m, n
      allocate(mat(n,m))
      read(iunit,*) ((mat(k,j),k=1,n),j=1,m)
      
      end subroutine readRealMatTransposeSize
************************************************************************
      subroutine readIntMatTransposeSize(iunit,mat)
      
C     Inputs: iunit --- specifier for file (to be written to)

C     Outputs: mat --- two-dimensional *integer* array
      
C     Purpose: Reads matrix from file, using size, storing it in transposed form.

      implicit none

C     input variables
      integer :: iunit
      
C     output variables
      integer, allocatable :: mat(:,:)
      
C     local variables
      integer :: j, k
      integer :: m, n
      
      read(iunit,*) m, n
      allocate(mat(n,m))
      read(iunit,*) ((mat(k,j),k=1,n),j=1,m)
      
      end subroutine readIntMatTransposeSize
************************************************************************
      subroutine writeRealVecSize(iunit,vec)
      
C     Inputs: iunit --- specifier for file (to be written to)
C             vec --- one-dimensional *real* array    

C     Outputs: None
      
C     Purpose: Writes size of vector and vector to a file, followed by space.

      implicit none

C     input variables
      real(dp) :: vec(:)
      integer :: iunit
      
C     local variables
      integer :: n
      
      n = size(vec)
      write(iunit,*) n
      call writeRealVec(iunit,vec)
      write(iunit,*) '' ! space for readability
      
      end subroutine writeRealVecSize
************************************************************************
      subroutine writeIntVecSize(iunit,vec)
      
C     Inputs: iunit --- specifier for file (to be written to)
C             vec --- one-dimensional *integer* array    

C     Outputs: None
      
C     Purpose: Writes size of vector and vector to a file, followed by space.

      implicit none

C     input variables
      integer :: vec(:)
      integer :: iunit
      
C     local variables
      integer :: n
      
      n = size(vec)
      write(iunit,*) n
      call writeIntVec(iunit,vec)
      write(iunit,*) '' ! space for readability
      
      end subroutine writeIntVecSize
************************************************************************
      subroutine writeRealVec(iunit,vec)
      
C     Inputs: iunit --- specifier for file (to be written to)
C             vec --- one-dimensional *real* array
            
C     Outputs: None
      
C     Purpose: Writes vector to a file.

      implicit none

C     input variables
      real(dp) :: vec(:)
      integer :: iunit
      
C     local variables
      integer :: i
      
      do i = 1, size(vec)
          write(iunit,*) vec(i)
      end do
      
      end subroutine writeRealVec
************************************************************************ 
      subroutine writeIntVec(iunit,vec)
      
C     Inputs: iunit --- specifier for file (to be written to)
C             vec --- one-dimensional *integer* array
            
C     Outputs: None
      
C     Purpose: Writes vector to a file.

      implicit none

C     input variables
      integer :: vec(:)
      integer :: iunit
      
C     local variables
      integer :: i
      
      do i = 1, size(vec)
          write(iunit,*) vec(i)
      end do
      
      end subroutine writeIntVec
************************************************************************
      subroutine writeRealMatSize(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *real* array    

C     Outputs: None
      
C     Purpose: Writes matrix (two-dimensional array) and its size
C     to a file, followed by space.

      implicit none

C     input variables
      integer :: iunit
      real(dp) :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      
      matsize = shape(mat)

      write(iunit,*) matsize(1), matsize(2)
      call writeRealMat(iunit,mat)
      write(iunit,*) '' ! space for readability

      end subroutine writeRealMatSize
************************************************************************ 
      subroutine writeIntMatSize(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *integer* array    

C     Outputs: None
      
C     Purpose: Writes matrix (two-dimensional array) and its size
C     to a file, followed by space.

      implicit none

C     input variables
      integer :: iunit
      integer :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      
      matsize = shape(mat)

      write(iunit,*) matsize(1), matsize(2)
      call writeIntMat(iunit,mat)
      write(iunit,*) '' ! space for readability

      end subroutine writeIntMatSize
************************************************************************
      subroutine writeRealMatTransposeSize(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *real* array

C     Outputs: None
      
C     Purpose: Writes transpose of matrix (two-dimensional array) and its size
C     to a file, followed by space.

      implicit none

C     input variables
      integer :: iunit
      real(dp) :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      
      matsize = shape(mat)

      write(iunit,*) matsize(2), matsize(1)
      call writeRealMatTranspose(iunit,mat)
      write(iunit,*) '' ! space for readability

      end subroutine writeRealMatTransposeSize
************************************************************************
      subroutine writeIntMatTransposeSize(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *integer* array

C     Outputs: None
      
C     Purpose: Writes transpose of matrix (two-dimensional array) and its size
C     to a file, followed by space.

      implicit none

C     input variables
      integer :: iunit
      integer :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      
      matsize = shape(mat)

      write(iunit,*) matsize(2), matsize(1)
      call writeIntMatTranspose(iunit,mat)
      write(iunit,*) '' ! space for readability

      end subroutine writeIntMatTransposeSize
************************************************************************
      subroutine writeRealMat(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *real* array   

C     Outputs: None
      
C     Purpose: Writes matrix (two-dimensional array) to a file.

      implicit none

C     input variables
      integer :: iunit
      real(dp) :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      integer :: i, j
      
      matsize = shape(mat)

      do i = 1,matsize(1)
          write(iunit,*) (mat(i,j), j=1,matsize(2))
      end do

      end subroutine writeRealMat
************************************************************************
      subroutine writeIntMat(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *integer* array   

C     Outputs: None
      
C     Purpose: Writes matrix (two-dimensional array) to a file.

      implicit none

C     input variables
      integer :: iunit
      integer :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      integer :: i, j
      
      matsize = shape(mat)

      do i = 1,matsize(1)
          write(iunit,*) (mat(i,j), j=1,matsize(2))
      end do

      end subroutine writeIntMat
************************************************************************
      subroutine writeRealMatTranspose(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *real* array   

C     Outputs: None
      
C     Purpose: Writes transpose of matrix (two-dimensional array) to a file.

      implicit none

C     input variables
      integer :: iunit
      real(dp) :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      integer :: i, j
      
      matsize = shape(mat)

      do i = 1,matsize(2)
          write(iunit,*) (mat(j,i), j=1,matsize(1))
      end do

      end subroutine writeRealMatTranspose
************************************************************************
      subroutine writeIntMatTranspose(iunit,mat)

C     Inputs: iunit --- specifier for file (to be written to)
C             mat --- two-dimensional *integer* array   

C     Outputs: None
      
C     Purpose: Writes transpose matrix (two-dimensional array) to a file.

      implicit none

C     input variables
      integer :: iunit
      integer :: mat(:,:)
      
C     local variables
      integer :: matsize(2)
      integer :: i, j
      
      matsize = shape(mat)

      do i = 1,matsize(2)
          write(iunit,*) (mat(j,i), j=1,matsize(1))
      end do

      end subroutine writeIntMatTranspose
************************************************************************        
      end module mod_utils