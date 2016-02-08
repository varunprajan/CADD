      module mod_pad_update
      
C     Purpose: Initialize, update positions of pad atoms
C     using total displacement field (FE (hat) + DD (tilde))

C     Notes: The presentation in Shilkrot et al., JMPS, 2004 is a bit misleading...
C     there, they appear to claim that just the FE displacements are needed.

      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_fe_elements, only: nfematerials, feelements
      use mod_fe_main_2d, only: getFEDispAtPoint
      use mod_disl_fields2, only: getTildeDispAtPointAll
      use mod_mesh_find, only: findInAllInitially
      implicit none
      
      private
      public :: padatoms
      
      type paddata
C     processed
      integer, allocatable :: element(:,:)
      real(dp), allocatable :: localpos(:,:)
      end type      

C     module variables
      type(paddata) :: padatoms
      
      contains
************************************************************************
      subroutine initPad()
      
C     local variables
      integer :: i
      integer :: node
      real(dp) :: xp, yp
      integer :: mnumfe
      integer :: element
      real(dp) :: r, s
      logical :: badflip
      
      allocate(padatoms%element(2,nodes%npadatoms))
      allocate(padatoms%localpos(2,nodes%npadatoms))
      do i = 1, nodes%npadatoms
          node = nodes%padatomlist(i)
          xp = nodes%posn(1,node)
          yp = nodes%posn(2,node)
          call findInAllInitially(xp,yp,mnumfe,element,r,s,badflip)
          if (badflip) then
              write(*,*) 'Unable to locate pad atom at', xp, yp
              stop
          end if    
          call assignPad(i,mnumfe,element,r,s)  
      end do
      
      end subroutine initPad
************************************************************************
      subroutine updatePad()
      
C     local variables
      integer :: i
      integer :: node
      integer :: mnumfe, element
      integer :: nelnodes
      real(dp) :: r, s
      real(dp) :: disphat(2), disptilde(2), disp(2)
      real(dp) :: posn(2)
      
      do i = 1, nodes%npadatoms
          node = nodes%padatomlist(i)
          mnumfe = padatoms%element(1,i)
          element = padatoms%element(2,i)
          r = padatoms%localpos(1,i)
          s = padatoms%localpos(2,i)
          nelnodes = feelements(mnumfe)%nelnodes
          posn = nodes%posn(1:2,node)
          disphat = getFEDispAtPoint(mnumfe,element,nelnodes,r,s)
          disptilde = getTildeDispAtPointAll(posn,mnumfe)
          disp = disphat + disptilde ! use total displacement to update pad
          nodes%posn(1:2,node) = posn + (disp - nodes%posn(4:5,node))
          nodes%posn(4:5,node) = disp
      end do    
      
      end subroutine updatePad
************************************************************************
      subroutine assignPad(padnum,mnumfe,element,r,s)
      
C     input variables
      integer :: padnum
      integer :: mnumfe
      integer :: element
      real(dp) :: r, s
      
      padatoms%element(1,padnum) = mnumfe
      padatoms%element(2,padnum) = element
      padatoms%localpos(1,padnum) = r
      padatoms%localpos(2,padnum) = s
      
      end subroutine assignPad
************************************************************************      
      
      end module