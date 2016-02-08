      module mod_mesh_find_brute
      
C     Purpose: Module for finding in which element in the mesh a
C     particular point lies. Constructs  auxiliary information needed for
C     this task: array with neighbors of a particular element.

      use mod_types, only: dp
      use mod_fe_elements, only: feelements, nfematerials
      use mod_math, only: getDist
      use mod_nodes, only: nodes
      implicit none
      
      private
      
      contains
************************************************************************

      end module