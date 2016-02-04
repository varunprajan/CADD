      module mod_fe_main_2d

C     Purpose: Initializes, assembles, factors, and solves finite element
C     equations for "hat" fields. Subtracts off dislocation fields ("tilde" fields)
C     to get correct boundary conditions.
C     Relies on MA57 from Harwell Subroutine Library to solve K*u = f for
C     sparse, symmetric, not necessarily positive definite K.
C     Enforces displacement boundary conditions using Lagrange multipliers.
      
C     Possible extensions/TODO: Nonzero traction boundary conditions,
C     multimaterial

      use mod_types, only: dp
      use mod_math, only: getUnitNormalRHR
      use mod_fe_el_2d, only: getK_2d, getB_2d, felib
      use mod_nodes, only: nodes
      use mod_materials, only: materials
      use mod_fe_elements, only: feelements, nfematerials, fematerials
      use mod_disl_fields2, only: getTildeDispAtPointAll,
     &                            getTildeStressAtPointAll
      implicit none
      
      private
      public :: solveAll, assembleAndFactorAll, solveOneMat,
     &  factorOneMat, initAssembly, countEqns, assembleKNormal,
     &  assembleLagrange, getDislForceRHS, getDislForceSub, getDispRHS,
     &  getFEStrainAtPoint, getFEDispAtPoint, updateFENodalPosn,
     &  assembly, getFEStressAtPoint, getTotalDispAtPoint,
     &  updateFENodalPosnAll
      
      type assemblydata
C     processed
      integer, allocatable :: rowindex(:)
      integer, allocatable :: colindex(:)
      real(dp), allocatable :: Ksparse(:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: fact(:)
      integer, allocatable :: ifact(:)
      integer :: neqnsnormal
      integer :: nentriesnormal
      integer :: neqnstot
      integer :: nentriestot
      end type assemblydata
      
C     module variables (local)
      type(assemblydata), allocatable :: assembly(:)
      real(dp), parameter :: lfactfac = 3.0_dp
      real(dp), parameter :: lifactfac = 3.0_dp
      real(dp), parameter :: lkeepfac = 1.1_dp
      real(dp), parameter :: lworkfac = 1.1_dp
      
      contains
************************************************************************ 
      subroutine solveAll()

C     Subroutine: solveAll

C     Inputs: None

C     Outputs: None

C     Purpose: Loop over fe materials, assembling right-hand-side (f)
C     and solving K*u = f for each. (Displacements are not updated in nodes%posn,
C     since this is, in general, costly, because dislocation fields are needed at every node.)
      
C     local variables
      integer :: i
      integer :: eltypenum
      
      do i = 1, nfematerials          
          assembly(i)%rhs = 0.0_dp ! initialize/wipe previous solution
          eltypenum = feelements(i)%eltypenum
          call getDislForceRHS(i,eltypenum)
          call getDispRHS(i)
          call solveOneMat(i)
      end do
      
      end subroutine solveAll
************************************************************************
      subroutine initAssembly()

C     Subroutine: initAssembly

C     Inputs: None

C     Outputs: None

C     Purpose: Allocate memory for K, rhs, etc. for each continuum material,
C     using the results of countEqns. Need only be done once.
      
C     local variables
      integer :: i
      integer :: neqnsnormal, nentriesnormal
      integer :: neqnstot, nentriestot
      
      allocate(assembly(nfematerials))
      do i = 1, nfematerials
          call countEqns(i,neqnsnormal,nentriesnormal,
     &                     neqnstot,nentriestot)
          assembly(i)%neqnsnormal = neqnsnormal
          assembly(i)%nentriesnormal = nentriesnormal
          assembly(i)%neqnstot = neqnstot
          assembly(i)%nentriestot = nentriestot
          allocate(assembly(i)%rowindex(nentriestot))
          allocate(assembly(i)%colindex(nentriestot))
          allocate(assembly(i)%Ksparse(nentriestot))
          allocate(assembly(i)%rhs(neqnstot))
      end do    
      
      end subroutine
************************************************************************
      subroutine assembleAndFactorAll()

C     Subroutine: assembleAndFactorAll

C     Inputs: None

C     Outputs: None

C     Purpose: Loop over fe materials, assembling K ("normal" K and Lagrange)
C     and factoring K for each.  Requires prior assembly of K, using
C     assembleAndFactorAll
      
C     local variables
      integer :: i
      integer :: mnum, eltypenum
      real(dp) :: C(3,3)
      
      do i = 1, nfematerials       
          mnum = fematerials%list(i)
          eltypenum = feelements(i)%eltypenum
          C = materials(mnum)%elconst          
          call assembleKNormal(i,eltypenum,C)
          call assembleLagrange(i)
          call factorOneMat(i)
      end do
      
      end subroutine assembleAndFactorAll
************************************************************************
      subroutine updateFENodalPosnAll()

C     Subroutine: updateFENodalPosnAll

C     Inputs: None

C     Outputs: None

C     Purpose: Loop over fe materials, updating positions and displacements
C     of fe nodes for each
      
C     Notes/TODO: May not work for multi-material (double update?)
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, nfematerials       
          call updateFENodalPosn(i)
      end do
      
      end subroutine updateFENodalPosnAll
************************************************************************  
      subroutine solveOneMat(mnumfe)

C     Subroutine: solveOneMat

C     Inputs: mnumfe --- number of continuum (FE) material

C     Outputs: None

C     Purpose: Using factorized K for a single continuum material, and
C     right-hand-side vector (f), solve K*u = f by calling MA57C (see
C     Harwell Subroutine Library MA57).

C     Notes: This relies on certain working arrays (e.g. fact), whose memory
C     must be allocated explicitly (see documentation to MA57). I use
C     certain constants to do this allocation, but they might be too small
C     for some problems...

      implicit none
      
C     input variables
      integer :: mnumfe
      
C     local variables
      real(dp) :: cntl(5)
      integer :: icntl(20)
      integer :: job
      integer :: neqns
      integer :: nrhs
      integer :: lfact, lifact
      integer :: lwork, lrhs
      real(dp), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      integer :: info(40)
      
      call MA57ID(cntl,icntl)
      job = 1
      nrhs = 1
      neqns = assembly(mnumfe)%neqnstot
      lrhs = neqns
      lfact = size(assembly(mnumfe)%fact)
      lifact = size(assembly(mnumfe)%ifact)
      lwork = ceiling(lworkfac*neqns*nrhs)
      allocate(work(lwork))
      allocate(iwork(neqns))
      call MA57CD(job,neqns,assembly(mnumfe)%fact,lfact,
     &                      assembly(mnumfe)%ifact,lifact,
     &                 nrhs,assembly(mnumfe)%rhs,lrhs,work,lwork,
     &                 iwork,icntl,info)
      
      end subroutine solveOneMat
************************************************************************
      subroutine factorOneMat(mnumfe)

C     Subroutine: factorOneMat

C     Inputs: mnumfe --- number of continuum (FE) material

C     Outputs: None

C     Purpose: Factor K for a single continuum material, using call to
C     MA57B (see Harwell Subroutine Library MA57). Need only be done
C     once, at beginning of run.

C     Notes: This relies on certain working arrays (e.g. keep), whose memory
C     must be allocated explicitly (see documentation to MA57). I use
C     certain constants to do this allocation, but they might be too small
C     for some problems...

      implicit none
      
C     input variables
      integer :: mnumfe
      
C     local variables
      real(dp) :: cntl(5)
      integer :: icntl(20)
      integer :: lkeep
      integer :: neqns, nentries
      integer, allocatable :: keep(:)
      integer, allocatable :: iwork(:)
      integer :: info(40)
      real(dp) :: rinfo(20)
      integer :: lfact
      integer :: lifact
      
      call MA57ID(cntl,icntl)
      neqns = assembly(mnumfe)%neqnstot
      nentries = assembly(mnumfe)%nentriestot
      lkeep = ceiling(lkeepfac*(5*neqns + nentries +
     &                          max(neqns,nentries) + 42))
      allocate(keep(lkeep))
      allocate(iwork(5*neqns))
      call MA57AD(neqns,nentries,assembly(mnumfe)%rowindex,
     &                           assembly(mnumfe)%colindex,
     &           lkeep,keep,iwork,icntl,info,rinfo)
      lfact = ceiling(lfactfac*info(9))
      lifact = ceiling(lifactfac*info(10))
      allocate(assembly(mnumfe)%fact(lfact))
      allocate(assembly(mnumfe)%ifact(lifact))
      call MA57BD(neqns,nentries,assembly(mnumfe)%Ksparse,
     &                           assembly(mnumfe)%fact,lfact,
     &                           assembly(mnumfe)%ifact,lifact,
     &            lkeep,keep,iwork,icntl,cntl,info,rinfo)
      
      end subroutine factorOneMat
************************************************************************
      subroutine countEqns(mnumfe,neqnsnormal,nentriesnormal,
     &                            neqnstot,nentriestot)

C     Subroutine: countEqns

C     Inputs: mnumfe --- number of continuum (FE) material

C     Outputs: neqnsnormal --- number of normal (i.e. non-Lagrange) equations
C                              i.e. rows in stiffness matrix
C              nentriesnormal --- number of normal (i.e. non-Lagrange) entries
C                              i.e. non-zero entries in stiffness matrix
C              neqnstot --- total number of equations
C              nentriestot --- total number of non-zero stiffness matrix entries

C     Purpose: Count up the number of nonzero entries, equations in stiffness
C     matrix, including Lagrange BCs, for the purposes of memory allocation
C     (see initAssembly).

      implicit none
      
C     input variables
      integer :: mnumfe
      
C     output variables
      integer :: neqnsnormal, nentriesnormal
      integer :: neqnstot, nentriestot
      
C     local variables
      integer :: i, node
      integer :: eltypenum
      integer :: neldof, nelements
      integer :: bcflag, nodetype
      logical :: xfixed, yfixed
      integer :: neqns, nentries
      
C     normal contributions (i.e. excluding lagrange)
      neqns = 2*size(feelements(mnumfe)%nodelist)
      neqnsnormal = neqns
      
C     each element contributes neldof diagonal entries
C     and neldof*(neldof-1)/2 off-diagonal entries
C     (just upper triangle)
      nelements = size(feelements(mnumfe)%connect,2)
      eltypenum = feelements(mnumfe)%eltypenum
      neldof = felib(eltypenum)%neldof
      nentries = nelements*(neldof + neldof*(neldof-1)/2)
      nentriesnormal = nentries
      
C     lagrange multipliers (from fixed displacements)    
      do i = 1, size(feelements(mnumfe)%bdnodelist)
          node = feelements(mnumfe)%bdnodelist(i)
          nodetype = nodes%types(2,node)
          bcflag = nodes%types(3,node)
C         interface atom or disp. fixed
          xfixed = ((nodetype==2).or.(bcflag==1).or.(bcflag==3))
          yfixed = ((nodetype==2).or.(bcflag==2).or.(bcflag==3))
          if (xfixed) then
              neqns = neqns + 1
              nentries = nentries + 1
          end if
          if (yfixed) then
              neqns = neqns + 1
              nentries = nentries + 1
          end if
      end do
      
      neqnstot = neqns
      nentriestot = nentries
      
      end subroutine countEqns
************************************************************************      
      subroutine assembleKNormal(mnumfe,eltypenum,C)

C     Subroutine: assembleKNormal

C     Inputs: mnumfe --- number of continuum (FE) material
C             eltypenum --- number of element type in fe element library
C             C --- stiffness matrix of material

C     Outputs: None

C     Purpose: Assemble "normal" entries in stiffness matrix K (i.e.
C     non-Lagrange), going element-by-element

      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: eltypenum
      real(dp) :: C(3,3)
      
C     local variables
      integer :: i, j, k1, k2
      real(dp) :: posn(2,felib(eltypenum)%nelnodes)
      integer :: dofnum(felib(eltypenum)%neldof)
      integer :: node, nodeidx
      real(dp) :: K(felib(eltypenum)%neldof,felib(eltypenum)%neldof)
      integer :: centries
      
      centries = 0
      do i = 1, feelements(mnumfe)%nfeelements
          do j = 1, felib(eltypenum)%nelnodes
              node = feelements(mnumfe)%connect(j,i)
C             use undeformed positions
              posn(:,j) = nodes%posn(1:2,node) - nodes%posn(4:5,node)
              nodeidx = feelements(mnumfe)%nodeinvlist(node)
              dofnum(2*j-1) = 2*nodeidx - 1
              dofnum(2*j) = 2*nodeidx
          end do
          K = getK_2d(posn,C,eltypenum)
          do k1 = 1, felib(eltypenum)%neldof
          do k2 = 1, felib(eltypenum)%neldof
C             just upper triangle
              if (k1 <= k2) then
                  centries = centries + 1
                  assembly(mnumfe)%rowindex(centries) = dofnum(k1)
                  assembly(mnumfe)%colindex(centries) = dofnum(k2)
                  assembly(mnumfe)%Ksparse(centries) = K(k1,k2)
              end if    
          end do
          end do
      end do
      
      end subroutine assembleKNormal
************************************************************************
      subroutine assembleLagrange(mnumfe)

C     Subroutine: assembleLagrange

C     Inputs: mnumfe --- number of continuum (FE) material

C     Outputs: None

C     Purpose: Assemble lagrange entries in stiffness matrix K (for
C     displacement boundary conditions)

      implicit none 

C     input variables
      integer :: mnumfe
      
C     local variables
      integer :: i
      integer :: node, nodeidx
      integer :: bcflag, nodetype
      logical :: xfixed, yfixed
      integer :: ceqns, centries
      
      ceqns = assembly(mnumfe)%neqnsnormal
      centries = assembly(mnumfe)%nentriesnormal
      do i = 1, size(feelements(mnumfe)%bdnodelist)
          node = feelements(mnumfe)%bdnodelist(i)
          nodeidx = feelements(mnumfe)%nodeinvlist(node)
          nodetype = nodes%types(2,node)
          bcflag = nodes%types(3,node)
C         interface atom or disp. fixed
          xfixed = ((nodetype==2).or.(bcflag==1).or.(bcflag==3))
          yfixed = ((nodetype==2).or.(bcflag==2).or.(bcflag==3))
          if (xfixed) then
              ceqns = ceqns + 1
              centries = centries + 1
              assembly(mnumfe)%rowindex(centries) = 2*nodeidx - 1
              assembly(mnumfe)%colindex(centries) = ceqns
              assembly(mnumfe)%Ksparse(centries) = 1
          end if
          if (yfixed) then
              ceqns = ceqns + 1
              centries = centries + 1
              assembly(mnumfe)%rowindex(centries) = 2*nodeidx
              assembly(mnumfe)%colindex(centries) = ceqns
              assembly(mnumfe)%Ksparse(centries) = 1
          end if
      end do
      
      end subroutine assembleLagrange
************************************************************************    
      subroutine getDislForceRHS(mnumfe,eltypenum)
      
C     Subroutine: getDislForceSub

C     Inputs: mnumfe --- number of continuum (FE) material
C             eltypenum --- number of element type in fe element library

C     Outputs: None

C     Purpose: Assemble nodal forces (f in K*u = f) due to dislocations
C     (recall that FE tractions must counteract DD (tilde) tractions, so
C     these DD forces are subtracted off.)

      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: eltypenum
      
C     local variables
      integer :: i, j
      integer :: node
      integer :: nodeidx(felib(eltypenum)%nedgenodes)
      integer :: nodetype, bcflag, idx
      logical :: allfixed
      logical :: xfixed(felib(eltypenum)%nedgenodes)
      logical :: yfixed(felib(eltypenum)%nedgenodes)
      real(dp) :: posn(2,felib(eltypenum)%nedgenodes)
      real(dp) :: fx(felib(eltypenum)%nedgenodes)
      real(dp) :: fy(felib(eltypenum)%nedgenodes)
      real(dp) :: sedgegauss(felib(eltypenum)%nedgeip)
      real(dp) :: wedgegauss(felib(eltypenum)%nedgenodes,
     &                       felib(eltypenum)%nedgeip)
      integer :: nedgenodes
      integer :: nedgeip
      
      sedgegauss = felib(eltypenum)%sedgegauss
      wedgegauss = felib(eltypenum)%wedgegauss
      nedgenodes = felib(eltypenum)%nedgenodes
      do i = 1, size(feelements(mnumfe)%bdedges,1)
C         get positions, normal
          allfixed = .true.
          do j = 1, nedgenodes
              node = feelements(mnumfe)%bdedges(i,j)
              nodetype = nodes%types(2,node)
              bcflag = nodes%types(3,node)
              xfixed(j) = ((nodetype==2).or.(bcflag==1).or.(bcflag==3))
              yfixed(j) = ((nodetype==2).or.(bcflag==1).or.(bcflag==3))
              nodeidx(j) = feelements(mnumfe)%nodeinvlist(node)
C             use undeformed positions
              posn(:,j) = nodes%posn(1:2,node) - nodes%posn(4:5,node)
              allfixed = (allfixed.and.xfixed(j).and.yfixed(j)) 
          end do
          
          if (.not.allfixed) then
              call getDislForceSub(mnumfe,nedgenodes,nedgeip,
     &                             sedgegauss,wedgegauss,posn,fx,fy)
C             add *opposite* to RHS (to cancel disl. field)
C             only if dof is unconstrained
              do j = 1, nedgenodes
                  idx = nodeidx(j)
C                 x-direction is unconstrained
                  if (.not.xfixed(j)) then
                      assembly(mnumfe)%rhs(2*idx - 1) =
     &                assembly(mnumfe)%rhs(2*idx - 1) - fx(j)
                  end if
C                 y-direction is unconstrained
                  if (.not.yfixed(j)) then
                      assembly(mnumfe)%rhs(2*idx) =
     &                assembly(mnumfe)%rhs(2*idx) - fy(j)
                  end if
              end do
          end if
      end do    

      end subroutine getDislForceRHS
************************************************************************       
      subroutine getDislForceSub(mnumfe,nedgenodes,nedgeip,
     &                           sedgegauss,wedgegauss,posn,fx,fy)

C     Subroutine: getDislForceSub

C     Inputs: mnumfe --- number of continuum (FE) material
C             nedgenodes --- number of FE nodes per edge
C             nedgeip --- number of integration poitns per edge
C             sedgegauss --- locations of integration points along edge (0 <= s <= 1)
C             wedgegauss --- weights relating nodal forces to tractions at integration points,
C                        nedgenodes by nedgeip (f_i = matmul(wgauss,traction_i))
C             posn --- array of nodal positions, 2 by nedgenodes,
C                      each column is a coordinate pair

C     Outputs: fx, fy --- forces in x/y-direction at edge nodes, length nedgenodes

C     Purpose: Get forces due to dislocations on edge nodes in two steps:
C              1) Evaluate traction due to dislocation fields at integration points
C              2) Relate these tractions to nodal forces using Gaussian quadrature
     
C     Notes: Gaussian quadrature makes this more accurate than evaluating
C            tractions at nodes and doing "element-by-element lumping") 

      implicit none

C     input variables
      integer :: mnumfe
      integer :: nedgenodes, nedgeip
      real(dp) :: sedgegauss(nedgeip)
      real(dp) :: wedgegauss(nedgenodes,nedgeip)
      real(dp) :: posn(2,nedgenodes)
      
C     output variables
      real(dp) :: fx(nedgenodes)
      real(dp) :: fy(nedgenodes)

C     local variables
      integer :: i
      real(dp) :: posnip(2)
      real(dp) :: posndiff(2), normal(2)
      real(dp) :: stress(3)
      real(dp) :: tracx(nedgeip)
      real(dp) :: tracy(nedgeip)
      real(dp) :: edgelen
 
      posndiff = posn(:,nedgenodes) - posn(:,1) ! vector connecting two extreme points
      normal = getUnitNormalRHR(posndiff)
      edgelen = posndiff(1)**2 + posndiff(2)**2
      
      do i = 1, nedgeip
          posnip = posn(:,1) + posndiff*sedgegauss(i)
          stress = getTildeStressAtPointAll(posnip,mnumfe)
          tracx(i) = stress(1)*normal(1) + stress(3)*normal(2)
          tracy(i) = stress(3)*normal(1) + stress(2)*normal(2)
      end do
      
      fx = edgelen*matmul(wedgegauss,tracx)
      fy = edgelen*matmul(wedgegauss,tracy)
      
      end subroutine getDislForceSub
************************************************************************      
      subroutine getDispRHS(mnumfe)

C     Subroutine: getDispRHS

C     Inputs: mnumfe --- number of continuum (FE) material

C     Outputs: None

C     Purpose: Get right-hand-sides of (Lagrange) displacement boundary conditions:
C     uhat = utot - utilde (tilde - infinite medium dislocation field,
C     hat - complementary FE field)

      implicit none
      
C     input variables
      integer :: mnumfe
      
C     local variables
      integer :: i
      integer :: node, nodeidx
      integer :: bcflag, nodetype
      logical :: xfixed, yfixed
      real(dp) :: posnnew(2), disptot(2)
      real(dp) :: dispdd(2), dispbc(2)
      integer :: ceqns
      
      ceqns = assembly(mnumfe)%neqnsnormal
      do i = 1, size(feelements(mnumfe)%bdnodelist)
          node = feelements(mnumfe)%bdnodelist(i)
          nodeidx = feelements(mnumfe)%nodeinvlist(node)
          nodetype = nodes%types(2,node)
          bcflag = nodes%types(3,node)
          xfixed = ((nodetype==2).or.(bcflag==1).or.(bcflag==3))
          yfixed = ((nodetype==2).or.(bcflag==2).or.(bcflag==3))
          if ((xfixed).or.(yfixed)) then
              posnnew = nodes%posn(1:2,node)
              disptot = nodes%posn(4:5,node)
C             tilde field
              dispdd = getTildeDispAtPointAll(posnnew-disptot,mnumfe) ! use orig. position
C             hat field
              dispbc = disptot - dispdd
              if (xfixed) then
                  ceqns = ceqns + 1
                  assembly(mnumfe)%rhs(ceqns) = dispbc(1)
              end if
              if (yfixed) then
                  ceqns = ceqns + 1
                  assembly(mnumfe)%rhs(ceqns) = dispbc(2)
              end if
          end if    
      end do 
      
      end subroutine getDispRHS
************************************************************************
      function getTotalDispAtPoint(posn,mnumfe,eltypenum,element,r,s)
     &                                                      result(disp)

C     Function: getTotalDispAtPoint

C     Inputs: posn --- (undeformed) coordinates of point (vector of length 2)
C             mnumfe --- number of continuum (FE) material
C             eltypenum --- number of element type in fe element library              
C             element --- number of FE element
C             r, s --- local coordinates of point of interest w.r.t. element

C     Outputs: displacement --- vector, length 2, of displacements at point of interest

C     Purpose: Get total displacement at point: sum of FE and DD displacements.
C     Used to update pad atom positions, and for creating dump file
      
      implicit none
      
C     input variables
      real(dp) :: posn(2)
      integer :: mnumfe
      integer :: eltypenum
      integer :: element
      real(dp) :: r, s
      
C     output variables
      real(dp) :: disp(2)
      
      disp = getFEDispAtPoint(mnumfe,eltypenum,element,r,s)
      disp = disp + getTildeDispAtPointAll(posn,mnumfe)
      
      end function getTotalDispAtPoint                                                      
************************************************************************
      function getFEDispAtPoint(mnumfe,eltypenum,element,r,s)
     &                                                      result(disp)

C     Function: getFEDispAtPoint

C     Inputs: mnumfe --- number of continuum (FE) material
C             eltypenum --- number of element type in fe element library              
C             element --- number of FE element
C             r, s --- local coordinates of point of interest w.r.t. element

C     Outputs: displacement --- vector, length 2, of displacements at point of interest

C     Purpose: Get displacement at point, using FE solution for displacements
C     (from assembly(mnumfe)%rhs)

      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: eltypenum
      integer :: element
      real(dp) :: r, s
      
C     output variables
      real(dp) :: disp(2)
      
C     local variables
      integer :: i
      real(dp), allocatable :: N(:)
      real(dp) :: udisp(felib(eltypenum)%nelnodes)
      real(dp) :: vdisp(felib(eltypenum)%nelnodes)
      integer :: node, nodeidx

      N = felib(eltypenum)%getN_2d_ptr(r,s)
      do i = 1, felib(eltypenum)%nelnodes
          node = feelements(mnumfe)%connect(i,element)
          nodeidx = feelements(mnumfe)%nodeinvlist(node)
          udisp(i) = assembly(mnumfe)%rhs(2*nodeidx - 1)
          vdisp(i) = assembly(mnumfe)%rhs(2*nodeidx)
      end do
      disp(1) = dot_product(N,udisp)
      disp(2) = dot_product(N,vdisp)
      
      end function getFEDispAtPoint
************************************************************************
      function getFEStressAtPoint(mnumfe,element,r,s) result(stress)

C     Function: getFEStressAtPoint

C     Inputs: mnumfe --- number of continuum (FE) material 
C             element --- number of FE element
C             r, s --- local coordinates of point of interest w.r.t. element

C     Outputs: stress --- vector, length 3, of stresses at point of interest.

C     Purpose: Get stress at point, using strain at point and elasticity matrix
      
      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: element
      real(dp) :: r, s
      
C     output variables
      real(dp) :: stress(3)
      
C     local variables
      integer:: mnum
      real(dp) :: C(3,3)
      integer :: eltypenum
      real(dp) :: strain(3)
     
      mnum = fematerials%list(mnumfe)
      C = materials(mnum)%elconst
      eltypenum = feelements(mnumfe)%eltypenum
      strain = getFEStrainAtPoint(mnumfe,eltypenum,element,r,s)
      stress = matmul(C,strain)
      
      end function getFEStressAtPoint
************************************************************************
      function getFEStrainAtPoint(mnumfe,eltypenum,element,r,s)
     &                                                   result(strain)

C     Function: getFEStrainAtPoint

C     Inputs: mnumfe --- number of continuum (FE) material
C             eltypenum --- number of element type in fe element library 
C             element --- number of FE element
C             r, s --- local coordinates of point of interest w.r.t. element

C     Outputs: strain --- vector, length 3, of strains at point of interest.
C                         Uses *engineering* notation.

C     Purpose: Get strain at point, using FE solution for displacements
C     (from assembly(mnumfe)%rhs), and strain-displacement matrix (B).

      implicit none
      
C     input variables
      integer :: mnumfe
      integer :: eltypenum
      integer :: element
      real(dp) :: r, s
      
C     output variables
      real(dp) :: strain(3)
      
C     local variables
      integer :: i
      real(dp) :: J(2,2)
      real(dp) :: posn(2,felib(eltypenum)%nelnodes)
      real(dp) :: B(3,felib(eltypenum)%neldof)
      integer :: node, nodeidx
      real(dp) :: disp(felib(eltypenum)%neldof)

      do i = 1, felib(eltypenum)%nelnodes
          node = feelements(mnumfe)%connect(i,element)
          posn(:,i) = nodes%posn(1:2,node) - nodes%posn(4:5,node) ! use undeformed positions
          nodeidx = feelements(mnumfe)%nodeinvlist(node)
          disp(2*i-1) = assembly(mnumfe)%rhs(2*nodeidx - 1)
          disp(2*i) = assembly(mnumfe)%rhs(2*nodeidx)
      end do
      call getB_2d(posn,r,s,eltypenum,B,J)
      strain = matmul(B,disp)
      
      end function getFEStrainAtPoint
************************************************************************
      subroutine updateFENodalPosn(mnumfe)

C     Subroutine: updateFENodalPosn

C     Inputs: mnumfe --- number of continuum (FE) material

C     Outputs: None

C     Purpose: Update FE nodal positions, displacements using total fields
C     (FE solution from assembly(mnumfe)%(rhs) and DD fields)

C     Notes: Technically the x-fixed and y-fixed checks are not necessary,
C     but using them ensures that the fixed points are truly fixed (otherwise,
C     small errors might accrue with many FE steps due to finite precision).

C     Other note: This subroutine is expected to be quite computationally
C     intensive because DD fields are requested at basically all FE nodes.
C     Also, it is not needed for main procedure. So, it should be invoked
C     only when a dump is requested.

      implicit none

C     input variables
      integer :: mnumfe
      
C     local variables
      integer :: i
      integer :: node
      integer :: nodetype, bcflag
      logical :: xfixed, yfixed
      real(dp) :: udisp, vdisp
      real(dp) :: posnundef(2)
      real(dp) :: dispdd(2)
      
      do i = 1, size(feelements(mnumfe)%nodelist)
          node = feelements(mnumfe)%nodelist(i)
          posnundef = nodes%posn(1:2,node) - nodes%posn(4:5,node)
          nodetype = nodes%types(2,node)
          bcflag = nodes%types(3,node)
          xfixed = ((nodetype==2).or.(bcflag==1).or.(bcflag==3))
          yfixed = ((nodetype==2).or.(bcflag==2).or.(bcflag==3))
          if (.not.(xfixed.and.yfixed)) then ! only call tilde routine if needed
              dispdd = getTildeDispAtPointAll(posnundef,mnumfe)
          end if    
          if (.not.xfixed) then
              udisp = assembly(mnumfe)%rhs(2*i - 1)
              udisp = udisp + dispdd(1)
              nodes%posn(1,node) = udisp + posnundef(1)
              nodes%posn(4,node) = udisp
          end if
          if (.not.yfixed) then
              vdisp = assembly(mnumfe)%rhs(2*i)
              vdisp = vdisp + dispdd(2)
              nodes%posn(2,node) = vdisp + posnundef(2)
              nodes%posn(5,node) = vdisp
          end if
      end do    
      
      end subroutine updateFENodalPosn
************************************************************************
      end module   
      