      module mod_integrate
      
C     Purpose: Has routines for integrating Newton's second law (sum of 
C     forces = ma) for atomistics. Currently, just the Velocity-Verlet
C     algorithm has been implemented (with possibility of damping).
C     
C     Possible extensions: Finite temperature (requires thermostatting),
C     non-uniform spatial damping (for pad).
      
      use mod_types, only: dp
      use mod_nodes, only: nodes
      use mod_materials, only: materials
      use mod_poten_assign, only: getPotForcesAll_ptr
      use mod_groups, only: groups, getGroupNum
      use mod_damping, only: getDampingForcesAll, dampingdata
      use mod_neighbors, only: neighbors, updateNeighborsCheck,
     &                 updateNeighborsNoCheck, updateNeighIncrementCurr
      implicit none
      
      private 
      public :: loopVerlet, velocityVerletSubAll,
     &          velocityVerletSub1, velocityVerletSub2
      
      contains
************************************************************************
      subroutine loopVerlet(increments,dt,gname,damp)

C     Subroutine: loopVerlet

C     Inputs: increments --- number of increments to run velocityVerlet
C             dt --- timestep
C             gname --- name of group for verlet update
C             damp --- structure containing damping information (group name, gamma, etc.)

C     Outputs: None

C     Purpose: Wrapper for velocityVerletSubAll. Used for "minimization"
C     by dislocation detection/passing module. Very similar to runNVE,
C     except no output.
      
      implicit none
      
C     input variables
      integer :: increments
      real(dp) :: dt
      character(len=*) :: gname
      type(dampingdata) :: damp
      
C     local variables
      integer :: i
      
C     initialize neighbors, forces      
      call updateNeighborsNoCheck()
      call getPotForcesAll_ptr()
      call getDampingForcesAll(damp)
      
C     run verlet loop
      do i = 1, increments
          call updateNeighIncrementCurr(1)
          call velocityVerletSubAll(dt,gname,damp)
          call updateNeighborsCheck()
      end do
      
      end subroutine loopVerlet
************************************************************************    
      subroutine velocityVerletSubAll(dt,gname,damp)
      
C     Subroutine: velocityVerletSubAll

C     Inputs: dt --- timestep
C             gname --- name of group for verlet update
C             damp --- structure containing damping information (group name, gamma, etc.)

C     Outputs: None

C     Purpose: Update positions and velocities of real atoms (not pad
C     atoms) in one increment

C     Notes: Only one pot. force call is necessary, since positions are
C     only updated once, and pot. forces don't depend on velocities.

C     Other notes: Technically, I think, the damping forces should be calculated
C     twice, once during step 3 and again during step 4, since velocity is changing.
C     However, this is not how lammps does it...so, the step 4 calculation is commented out
      
      implicit none
      
C     input variables
      real(dp) :: dt
      character(len=*) :: gname
      type(dampingdata) :: damp

C     local variables
      integer :: gnum

      gnum = getGroupNum(gname)

C     update velocities and positions (steps 1 and 2)
      call velocityVerletSub1(0.5_dp*dt,gnum)
      call velocityVerletSub2(dt,gnum)
      
C     get potential and damping forces (positions and velocities have changed)
C     (step 3)
      call getPotForcesAll_ptr()
      call getDampingForcesAll(damp)   

C     update velocities (step 4)     
      call velocityVerletSub1(0.5_dp*dt,gnum)
      
C     get damping forces only (velocities have changed)
C     call getDampingForcesAll(dgnum,gamma,dampflag)
      
      end subroutine velocityVerletSubAll
************************************************************************
      subroutine velocityVerletSub1(dt,gnum)

C     Subroutine: velocityVerletSub1

C     Inputs: dt --- timestep
C             gnum --- number of group (velocityVerlet is applied only
C             to atoms in group)

C     Outputs: None

C     Purpose: Update velocities of real atoms (not pad atoms) for
C     increment dt
      
C     Notes: Assumes forces due to other atoms (pot. forces) have
C     already been calculated. Also note that, in the velocity Verlet
C     algorithm, this subroutine is called with argument dt/2

      implicit none

C     input variables
      real(dp) :: dt
      integer :: gnum
      
C     local variables
      integer :: i, atom, atommat, bcflag
      real(dp) :: mass
      real(dp) :: dvel(2), force(2)
      
      do i = 1, nodes%natoms
C     only do operation for atoms in group
      if (groups(gnum)%maskatoms(i)) then
          atom = nodes%atomlist(i)
          atommat = nodes%types(1,atom)
          bcflag = nodes%types(3,atom)
C         check if atom is free to move
C         note: pad atom is expected to be completely fixed (bcflag = 3)
          if (bcflag /= 3) then
              mass = materials(atommat)%mass
              force = nodes%potforces(:,i) + nodes%dampforces(:,i)
              dvel = force*dt/mass
              
C             if x is fixed
              if (bcflag == 1) then
                  dvel(1) = 0.0_dp
C             if y is fixed
              else if (bcflag == 2) then
                  dvel(2) = 0.0_dp
              end if
              
              nodes%posn(6:7,atom) = nodes%posn(6:7,atom) + dvel
          end if
      end if
      end do
      
      end subroutine velocityVerletSub1
************************************************************************
      subroutine velocityVerletSub2(dt,gnum)
      
C     Subroutine: velocityVerletSub2

C     Inputs: dt --- timestep
C             gnum --- number of group (velocityVerlet is applied only
C             to atoms in group)

C     Outputs: None

C     Purpose: Update positions of real atoms (not pad atoms) for
C     increment dt using the velocities
      
      implicit none
      
C     input variables
      real(dp) :: dt
      integer :: gnum
      
C     local variables
      integer :: i, atom, bcflag
      real(dp) :: dpos(2)
      
      do i = 1, nodes%natoms
C     only do operation for atoms in group
      if (groups(gnum)%maskatoms(i)) then          
          atom = nodes%atomlist(i)
          bcflag = nodes%types(3,atom)
C         check if atom is free to move
C         note: pad atom is expected to be completely fixed (bcflag = 3)
          if (bcflag /= 3) then
              dpos = nodes%posn(6:7,atom)*dt
              
C             if x is fixed
              if (bcflag == 1) then
                  dpos(1) = 0.0_dp
C             if y is fixed
              else if (bcflag == 2) then
                  dpos(2) = 0.0_dp
              end if
              
              nodes%posn(1:2,atom) = nodes%posn(1:2,atom) + dpos
              nodes%posn(4:5,atom) = nodes%posn(4:5,atom) + dpos
          end if
      end if    
      end do

      end subroutine velocityVerletSub2
************************************************************************
      
      end module