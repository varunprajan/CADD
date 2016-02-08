      module mod_poten_pair_table

C     Purpose: Has routines for computing forces for pair potential from
C     linear interpolation of tabular data (force vs. interatomic separation)
C     
C     Possible extensions: Compute energies? Develop other poten_pair styles? (e.g. Morse/LJ)
      
      use mod_types, only: dp
      use mod_potentials, only: potentials
      use mod_neighbors, only: neighbors
      use mod_nodes, only: nodes
      use mod_math, only: linearInterp
      use mod_interactions, only: interactions
      implicit none
      
      private
      public :: getPotForceTable, getPotForcesAllTable
      
      contains

************************************************************************
      subroutine getPotForcesAllTable()
      
C     Subroutine: getPotForcesAllTable

C     Inputs: None

C     Outputs: None

C     Purpose: Loop through all atoms in group, compute forces using helper
C     function getForce, assign results to nodes%atomforces

C     Notes: Implementation of periodic images may need to be improved...
C     are both Lz and images necessary? (Of course, none of this matters
C     for 2D hex lattice, where Lz is arbitrary and images = 0)
      
      implicit none
      
C     local variables
      integer :: i, j, k, neigh
      integer :: atom, neighatom
      integer :: atommat, neighmat
      integer :: atomtype, neightype
      integer :: pnum
      real(dp) :: atompos(3), neighpos(3), neighposnew(3)
      real(dp) :: forcevec(3)
      logical :: getneighborforce, getatomforce
      
      nodes%potforces = 0.0_dp
      do i = 1, nodes%natoms
          atom = nodes%atomlist(i)
          atommat = nodes%types(1,atom)
          atomtype = nodes%types(2,atom)
          atompos = nodes%posn(1:3,atom)
          do j = 1, neighbors%neighcount(i) ! loop over neighbors
              neigh = neighbors%neighlist(j,i)
              neighatom = nodes%atomlist(neigh)
              neighmat = nodes%types(1,neighatom)
              neightype = nodes%types(2,neighatom)
              neighpos = nodes%posn(1:3,neighatom)
              pnum = interactions%mat(atommat,neighmat)
              do k = -neighbors%images, neighbors%images ! periodic images
                  getneighborforce = ((k==0).and.(i<neigh)) ! if neighbor is in central layer, and use symmetry (and no self-interaction)
                  getatomforce = (getneighborforce.or.
     &                            ((k/=0).and.(atomtype/=-1))) ! don't use out-of-plane images for pad atoms
                  if (getatomforce) then
                      neighposnew = neighpos
                      neighposnew(3) = neighpos(3) + k*neighbors%Lz
                      forcevec = getPotForceTable(pnum,atompos,
     &                                                      neighposnew)
                      if (atomtype/=-1) then ! not a pad
                          nodes%potforces(:,i) =
     &                    nodes%potforces(:,i) + forcevec(1:2)
                      end if
                      if (getneighborforce) then
                          if (neightype/=-1) then ! not a pad
                              nodes%potforces(:,neigh) =
     &                        nodes%potforces(:,neigh) - forcevec(1:2)
                          end if
                      end if
                  end if
              end do 
          end do 
      end do
          
      end subroutine getPotForcesAllTable
************************************************************************  
      function getPotForceTable(pnum,atompos,neighpos) result(forcevec)
      
C     Function: getPotForceTable

C     Inputs: None

C     Outputs: None

C     Purpose: Computes force vector between atom and its neighbor, using
C     linear interpolation (table lookup)

C     Algorithm: Using # of potential and distance between atoms,
C     consults appropriate potential table to get force magnitude,
C     Uses the fact that the force is a central force to get its direction.
      
      implicit none
      
C     input variables
      integer :: pnum
      real(dp) :: atompos(3), neighpos(3)
      
C     output variables
      real(dp) :: forcevec(3)
      
C     local variables
      real(dp) :: r, forcemag
      real(dp) :: delpos(3)
      integer :: n
      
      delpos = atompos - neighpos
      r = sqrt(sum(delpos**2))
      if (r > potentials(pnum)%forcecutoff) then
          forcevec = 0.0_dp
      else    
          n = size(potentials(pnum)%pottable,1)
          forcemag = linearInterp(r,potentials(pnum)%pottable(:,1),
     &                              potentials(pnum)%pottable(:,3),n)
          forcevec = forcemag/r*delpos
      end if    
      
      end function getPotForceTable
************************************************************************
      end module mod_poten_pair_table