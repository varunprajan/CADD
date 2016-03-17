      module mod_damping
      
C     Purpose: Reads/writes/stores information about
C     damping (i.e. whether it's on or off, damping coefficient)
      
C     Possible extensions: Add more fields for more complicated types of
C     damping (stadium damping, finite temperature)
      
      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_groups, only: groups, getGroupNum
      use mod_math, only: logicalToInt, intToLogical
      implicit none
      
      private 
      public :: getDampingForcesAll, initDampingData, writeDampingData,
     &          getDampingForce, readDampingData, normaldamping,
     &          dampingdata, readDampingDataSub, writeDampingDataSub
      
      type dampingdata
C     read-in
      logical :: flag
      real(dp) :: gamma
      character(len=60) :: gname
      end type

C     module variables      
      type(dampingdata) :: normaldamping
      
      contains
************************************************************************
      subroutine initDampingData(dampingfile)

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

C     Inputs: dampingfile --- filename where damping data is stored
C     (should be something like '[filepref]_damping')

C     Outputs: None

C     Purpose: Read data from file into "normaldamping" structure
      
      implicit none
      
C     input variables
      character(len=*) :: dampingfile
      
C     local variables
      integer :: iunit
      
      open(newunit=iunit,file=dampingfile)
      
      call readDampingDataSub(iunit,normaldamping)
      
      close(iunit)
      
      end subroutine readDampingData
************************************************************************
      subroutine readDampingDataSub(iunit,damp)

C     Inputs: iunit --- integer file specifier
             
C     In/out: damp --- structure containing damping data

C     Outputs: None

C     Purpose: Read data from file into damp structure (note: also used in mod_disl_detect_pass)
      
      implicit none
      
C     input variables
      integer :: iunit
      
C     in/out variables
      type(dampingdata) :: damp
      
C     local variables
      integer :: temp
      
      read(iunit,*) temp
C     use explicit integer -> logical conversion
      damp%flag = intToLogical(temp)
      read(iunit,*) damp%gamma
      read(iunit,*) damp%gname
      
      end subroutine readDampingDataSub
************************************************************************
      subroutine writeDampingData(dampingfile)

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
      
      open(newunit=iunit,file=dampingfile)
      
      call writeDampingDataSub(iunit,normaldamping)
      
      close(iunit)
      
      end subroutine writeDampingData
************************************************************************
      subroutine writeDampingDataSub(iunit,damp)

C     Inputs: iunit --- integer file specifier
             
C     In/out: damp --- structure containing damping data

C     Outputs: None

C     Purpose: Write data from file into damp structure (note: also used in mod_disl_detect_pass)
      
      implicit none
      
C     input variables
      integer :: iunit
      type(dampingdata) :: damp
      
C     local variables
      integer :: temp
      
C     use explicit logical -> integer conversion
      temp = logicalToInt(damp%flag)    
      write(iunit,*) temp
      write(iunit,*) damp%gamma
      write(iunit,*) damp%gname
      
      end subroutine
************************************************************************
      subroutine getDampingForcesAll(damp)

C     Inputs: damp --- structure containing damping information (group name, gamma, etc.)

C     Outputs: None

C     Purpose: Compute damping forces for all atoms in group damp%gname
      
      implicit none
      
C     input variables
      type(dampingdata) :: damp
      
C     local variables
      integer :: i, atom
      real(dp) :: posn(2), vel(2)
      integer :: gnum
      real(dp) :: gamma
      
      nodes%dampforces = 0.0_dp
      gnum = getGroupNum(damp%gname)
      gamma = damp%gamma
      if (damp%flag) then
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