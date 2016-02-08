      module mod_poten_main
      
      use mod_misc, only: misc
      use mod_poten_pair_table, only: getPotForcesAllTable
      
      implicit none
      private
      public :: getPotForcesAll_ptr, assignPotStyle
      
      procedure(Dummy), pointer :: getPotForcesAll_ptr => NULL()      
      
      contains
************************************************************************
      subroutine assignPotStyle()
      
      select case (misc%potstyle)
          case ('table')
              getPotForcesAll_ptr => getPotForcesAllTable
          case default
              write(*,*) 'Potential style has not yet been defined'
              stop
      end select
      
      end subroutine assignPotStyle
************************************************************************
      subroutine Dummy()
      
      end subroutine Dummy
************************************************************************      
      
      end module mod_poten_main