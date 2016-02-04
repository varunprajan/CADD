      module mod_types
      
C     Purpose: Sets the precision for all "float" type variables/arrays
C     used in the program, through the variable "dp".
      
      implicit none
      integer, parameter :: dp = kind(0.0d0)
      
      private
      public :: dp
      
      end module