      module mod_fe_main_2d_assign
      
C     Purpose: Assigns pointers to FE subroutines depending on whether
C     simulation has no dislocations (fe, cadd_nodisl) or dislocations
C     (dd, cadd). This will change, for instance, how the right hand side 
C     is assembled, how the nodal displacements are obtained, etc.
C     Compare mod_fe_main_2d (with dislocations) and mod_fe_main_2d_no_disl
C     (no dislocations)
      
      use mod_types, only: dp
      use mod_misc, only: misc
      use mod_fe_main_2d, only: solveAll, assembleAndFactorAll,
     & getTotalDispAtPoint, updateFENodalPosnAll, initAssembly
      use mod_fe_main_2d_no_disl, only: solveAllNoDisl,
     & assembleAndFactorAllNoDisl, getTotalDispAtPointNoDisl,
     & updateFENodalPosnAllNoDisl, initAssemblyNoDisl
     
      private
      public :: updateFENodalPosnAll_ptr, getTotalDispAtPoint_ptr,
     &          assignFE, solveAll_ptr, assembleAndFactorAll_ptr,
     &          initAssembly_ptr

      procedure(Dummy), pointer :: solveAll_ptr => NULL()
      procedure(Dummy2), pointer :: assembleAndFactorAll_ptr => NULL()
      procedure(Dummy3), pointer :: updateFENodalPosnAll_ptr => NULL()
      procedure(Dummy4), pointer :: getTotalDispAtPoint_ptr => NULL()
      procedure(Dummy5), pointer :: initAssembly_ptr => NULL()
      
      contains
************************************************************************
      subroutine assignFE(simtype)

C     Inputs: simtype --- simulation type: atomistic, fe, dd, cadd_nodisl, or cadd

C     Outputs: None

C     Purpose: Assigns pointers to FE subroutines depending on whether
C     simulation has no dislocations (fe, cadd_nodisl) or dislocations
C     (dd, cadd).
      
      implicit none
      
C     input variables
      character(len=*) :: simtype
      
      select case (trim(simtype))
          case ('cadd','dd')
              solveAll_ptr => solveAll
              assembleAndFactorAll_ptr => assembleAndFactorAll
              updateFENodalPosnAll_ptr => updateFENodalPosnAll
              getTotalDispAtPoint_ptr => getTotalDispAtPoint
              initAssembly_ptr => initAssembly
          case ('cadd_nodisl','fe')
              solveAll_ptr => solveAllNoDisl
              assembleAndFactorAll_ptr => assembleAndFactorAllNoDisl
              updateFENodalPosnAll_ptr => updateFENodalPosnAllNoDisl
              getTotalDispAtPoint_ptr => getTotalDispAtPointNoDisl
              initAssembly_ptr => initAssemblyNoDisl
          case ('atomistic')
              return    
          case default
              write(*,*) 'Simulation type has not yet been defined'
              stop
      end select
      
      end subroutine
************************************************************************
      subroutine Dummy()
      
      end subroutine Dummy
************************************************************************ 
      subroutine Dummy2()
      
      end subroutine Dummy2
************************************************************************
      subroutine Dummy3()
      
      end subroutine Dummy3
************************************************************************
      function Dummy4(posn,mnumfe,eltypenum,element,r,s) result(disp)
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      integer :: eltypenum
      integer :: element
      real(dp) :: r, s
      
C     output variables
      real(dp) :: disp(2)
            
      end function Dummy4                                                     
************************************************************************
      subroutine Dummy5()
      
      end subroutine Dummy5
************************************************************************   
      end module