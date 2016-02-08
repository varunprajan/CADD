      module mod_damping
      
C     Purpose: Reads/writes/stores information about
C     damping (i.e. whether it's on or off, damping coefficient)
      
C     Possible extensions: Add more fields for more complicated types of
C     damping (stadium damping, finite temperature)
      
      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_groups, only: groups, nextragroups
      implicit none
      
      private 
      public :: getDampingForcesAll, initDampingData, writeDampingData,
     &          getDampingForce, readDampingData, damping
      
      type dampingdata
C     read-in
      logical :: flag
      real(dp) :: gamma
      integer :: gnum
      end type

C     module variables      
      type(dampingdata) :: damping
      
      contains
************************************************************************
      subroutine initDampingData(dampingfile)
      
C     Subroutine: initDampingData

C     Inputs: dampingfile --- filename where damping data is stored
C     (should be something like '[filepref]_damping')

C     Outputs: None

C     Purpose: Read, initialize data in "damping" structure.
      
      implicit none
      
C     input variables
      character(len=*) :: dampingfile
      
      call readDampingData(dampingfile)
      
      end subroutine initDampingData
************************************************************************
      subroutine readDampingData(dampingfile)

C     Subroutine: readDampingData

C     Inputs: dampingfile --- filename where damping data is stored
C     (should be something like '[filepref]_damping')

C     Outputs: None

C     Purpose: Read data from file into "damping" structure
      
      implicit none
      
C     input variables
      character(len=*) :: dampingfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=dampingfile)
      
      read(iunit,*) temp
C     use explicit integer -> logical conversion
      damping%flag = (temp /= 0)
      read(iunit,*) damping%gamma
      read(iunit,*) temp
      if (temp == 0) then ! no group specified
          damping%gnum = 1 ! all group
      else    
          damping%gnum = temp + nextragroups ! account for two extra groups
      end if    
      
      close(iunit)
      
      end subroutine readDampingData
************************************************************************
      subroutine writeDampingData(dampingfile)

C     Subroutine: writeDampingData

C     Inputs: dampingfile --- filename where damping data is stored
C     (should be something like '[filepref]_damping')

C     Outputs: None

C     Purpose: Write damping data to file (essentially
C     inverse of readDampingData). Useful in creating "restart" file
C     (i.e. resulting file should be able to be read in by readDampingData)

      implicit none
      
C     input variables
      character(len=*) :: dampingfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=dampingfile)
      
C     use explicit logical -> integer conversion
      if (damping%flag) then
          temp = 1
      else
          temp = 0
      end if    
      write(iunit,*) temp
      write(iunit,*) damping%gamma
      if (damping%gnum == 1) then
          temp = 0
      else
          temp = damping%gnum - nextragroups ! account for extra groups
      end if
      write(iunit,*) temp
      
      close(iunit)
      
      end subroutine writeDampingData
************************************************************************
      subroutine getDampingForcesAll(gnum,gamma,dampflag)

C     Subroutine: getDampingForcesAll

C     Inputs: gnum --- group number for which damping is implemented
C             gamma --- damping coefficient
C             dampflag --- flag indicating whether damping is on

C     Outputs: None

C     Purpose: Compute damping forces for all atoms in group gnum
      
      implicit none
      
C     input variables
      integer :: gnum
      real(dp) :: gamma
      logical :: dampflag
      
C     local variables
      integer :: i, atom
      real(dp) :: posn(2), vel(2)
      
      nodes%dampforces = 0.0_dp
      if (dampflag) then
          do i = 1, nodes%natoms
C         only do operation for atoms in group
          if (groups(gnum)%maskatoms(i)) then
              atom = nodes%atomlist(i)
              posn = nodes%posn(1:2,atom)
              vel = nodes%posn(6:7,atom)
              nodes%dampforces(:,i) = getDampingForce(posn,vel,gamma)
          end if
          end do
      end if
      
      end subroutine getDampingForcesAll
************************************************************************     
      function getDampingForce(posn,vel,gamma) result(force)
      
C     Function: getDampingForce

C     Inputs: posn --- vector, 2 by 1 [x, y]
C             vel --- vector, 2 by 1 [vx,vy]
C             gamma --- damping coefficient

C     Outputs: force --- vector, 2 by 1 [fx,fy]

C     Purpose: Compute damping force from damping coefficient, atom
C     velocity. Currently just for athermal situation, constant damping
C     coefficient
      
C     Notes: Would need some modification for finite temperature
C     to include random Langevin-type kicks (see Qu et al. 2005 paper)

C     TODO: Does nonconst damping need to be implemented for the pad?

      implicit none

C     input variables
      real(dp) :: posn(2)
      real(dp) :: vel(2)
      real(dp) :: gamma
      
C     output variables
      real(dp) :: force(2)

      force = -gamma*vel
      
      end function getDampingForce
************************************************************************
      end module