      module mod_poten_assign
      
C     Purpose: Assigns pointers to getPotForcesAll based on style of potential
C     (misc%potstyle). This module will support implementation
C     of additional potential styles in modules (like mod_poten_pair_table).
      
      use mod_misc, only: misc
      use mod_poten_pair_table, only: getPotForcesAllTable
      
      implicit none
      private
      public :: getPotForcesAll_ptr, initPotentialStyle
      
      procedure(Dummy), pointer :: getPotForcesAll_ptr => NULL()      
      
      contains
************************************************************************
      subroutine initPotentialStyle()
 
C     Subroutine: initPotentialStyle
 
C     Inputs: None
 
C     Outputs: None
 
C     Purpose: Assigns getPotForcesAll_ptr to function, based on style
C     of potential (misc%potstyle). Currently only implemented for tabular
C     pair potential.
      
      select case (trim(misc%potstyle))
          case ('table')
              getPotForcesAll_ptr => getPotForcesAllTable
          case default
              write(*,*) 'Potential style has not yet been defined'
              stop
      end select
      
      end subroutine initPotentialStyle
************************************************************************
      subroutine Dummy()
      
      end subroutine Dummy
************************************************************************
      end module mod_poten_assign