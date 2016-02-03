      module mod_pad_atoms
      
C     Purpose: Initialize, update positions of pad atoms
C     using total displacement field (FE (hat) + DD (tilde))

C     Notes: The presentation in Shilkrot et al., JMPS, 2004 is a bit misleading...
C     there, they appear to claim that just the FE displacements are needed.

      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_fe_elements, only: nfematerials, feelements
      use mod_fe_main_2d_assign, only: getTotalDispAtPoint_ptr
      use mod_mesh_find, only: findInAllInitially
      implicit none
      
      private
      public :: padatoms, initPad, updatePad, assignPad
      
      type padt
      integer :: mnumfe
      integer :: element
      real(dp) :: localpos(2)
      end type

C     module variables
      type(padt), allocatable :: padatoms(:)
      
      contains
************************************************************************
      subroutine initPad()
      
C     Subroutine: initPad

C     Inputs: None

C     Outputs: None

C     Purpose: Read, initialize data in "padatoms" structure, which holds
C     information about pad atoms (atoms within FE region, whose positions
C     are adjusted during FE/DD step, and held fixed during atomistic step)
      
C     local variables
      integer :: i
      integer :: node
      real(dp) :: xp, yp
      integer :: mnumfe
      integer :: element
      real(dp) :: r, s
      logical :: badflip
      
      allocate(padatoms(nodes%npadatoms))
      do i = 1, nodes%npadatoms
          node = nodes%padatomlist(i)
          xp = nodes%posn(1,node)
          yp = nodes%posn(2,node)
          call findInAllInitially(xp,yp,mnumfe,element,r,s,badflip)
          if (badflip) then
              write(*,*) 'Unable to locate pad atom at', xp, yp
              stop
          end if    
          call assignPad(padatoms(i),mnumfe,element,r,s)  
      end do
      
      end subroutine initPad
************************************************************************
      subroutine assignPad(padatom,mnumfe,element,r,s)
      
C     Subroutine: assignPad

C     Inputs: padatom --- container for a particular pad atom
C             mnumfe --- fe material number for pad atom
C             element --- fe element number for pad atom
C             r, s --- local position of pad atom within fe element

C     Outputs: None

C     Purpose: Assign attributes to pad atom
      
C     input variables
      integer :: mnumfe
      integer :: element
      real(dp) :: r, s
      
C     in/out variables
      type(padt) :: padatom
      
      padatom%mnumfe = mnumfe
      padatom%element = element
      padatom%localpos(1) = r
      padatom%localpos(2) = s
      
      end subroutine assignPad
************************************************************************
      subroutine updatePad()
      
C     Subroutine: updatePad

C     Inputs: None

C     Outputs: None

C     Purpose: Update positions of all pad atoms according to combined FE/DD displacement
      
C     local variables
      integer :: i
      integer :: node
      integer :: mnumfe, eltypenum, element
      real(dp) :: r, s
      real(dp) :: disp(2)
      real(dp) :: posnundef(2)
      
      do i = 1, nodes%npadatoms
          node = nodes%padatomlist(i)
          mnumfe = padatoms(i)%mnumfe
          element = padatoms(i)%element
          r = padatoms(i)%localpos(1)
          s = padatoms(i)%localpos(2)
          posnundef = nodes%posn(1:2,node) - nodes%posn(4:5,node)
          eltypenum = feelements(mnumfe)%eltypenum
          disp = getTotalDispAtPoint_ptr(posnundef,
     &                                   mnumfe,eltypenum,element,r,s)
          nodes%posn(1:2,node) = posnundef + disp ! undeformed position + disp
          nodes%posn(4:5,node) = disp
      end do    
      
      end subroutine updatePad
************************************************************************     
      
      end module