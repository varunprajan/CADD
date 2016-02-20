      module mod_fe_main_2d_no_disl

C     Purpose: Initializes, assembles, factors, and solves finite element
C     equations, for "bare" FE problem (no dislocations).
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
      implicit none
      
      private
      public :: solveAllNoDisl, assembleAndFactorAllNoDisl, solveOneMat,
     &  factorOneMat, initAssemblyNoDisl, countEqns, assembleKNormal,
     &  assembleLagrange, getDispRHS, getTotalDispAtPointNoDisl,
     &  getFEStrainAtPoint, getFEDispAtPoint, updateFENodalPosn,
     &  assembly, getFEStressAtPoint, updateFENodalPosnAllNoDisl
      
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
      
C     HARD-CODED CONSTANTS
      real(dp), parameter :: lfactfac = 1.25_dp
      real(dp), parameter :: lifactfac = 1.25_dp
      real(dp), parameter :: lkeepfac = 1.25_dp
      real(dp), parameter :: lworkfac = 1.25_dp
      
      contains
************************************************************************ 
      subroutine solveAllNoDisl()

C     Function: solveAllNoDisl

C     Inputs: None

C     Outputs: None

C     Purpose: Loop over fe materials, assembling right-hand-side (f)
C     and solving K*u = f for each. Requires prior assembly of K, using
C     assembleAndFactorAllNoDisl
      
C     local variables
      integer :: i
      
      do i = 1, nfematerials          
          assembly(i)%rhs = 0.0_dp ! initialize/wipe previous solution
          call getDispRHS(i)
          call solveOneMat(i)
      end do
      
      end subroutine solveAllNoDisl
************************************************************************
      subroutine initAssemblyNoDisl()

C     Function: initAssemblyNoDisl

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
      
      end subroutine initAssemblyNoDisl
************************************************************************
      subroutine assembleAndFactorAllNoDisl()

C     Function: assembleAndFactorAllNoDisl

C     Inputs: None

C     Outputs: None

C     Purpose: Loop over fe materials, assembling K ("normal" K and Lagrange)
C     and factoring K for each.
      
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
      
      end subroutine assembleAndFactorAllNoDisl
************************************************************************
      subroutine updateFENodalPosnAllNoDisl()

C     Subroutine: updateFENodalPosnAllNoDisl

C     Inputs: None

C     Outputs: None

C     Purpose: Loop over fe materials, updating positions and displacements
C     of fe nodes for each
      
      implicit none
      
C     local variables
      integer :: i
      
      do i = 1, nfematerials       
          call updateFENodalPosn(i)
      end do
      
      end subroutine updateFENodalPosnAllNoDisl
************************************************************************ 
      subroutine solveOneMat(mnumfe)

C     Function: solveOneMat

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

C     Function: factorOneMat

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

C     Function: countEqns

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

C     Function: assembleKNormal

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

C     Function: assembleLagrange

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
      subroutine getDispRHS(mnumfe)

C     Function: getDispRHS

C     Inputs: mnumfe --- number of continuum (FE) material

C     Outputs: None

C     Purpose: Get right-hand-sides of (Lagrange) displacement boundary conditions

      implicit none
      
C     input variables
      integer :: mnumfe
      
C     local variables
      integer :: i
      integer :: node, nodeidx
      integer :: bcflag, nodetype
      logical :: xfixed, yfixed
      real(dp) :: posnnew(2)
      real(dp) :: dispbc(2)
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
              dispbc = nodes%posn(4:5,node)
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
      function getTotalDispAtPointNoDisl(posn,mnumfe,eltypenum,element,
     &                                   r,s) result(disp)

C     Function: getTotalDispAtPointNoDisl

C     Inputs: posn --- (undeformed) coordinates of point (vector of length 2)
C             mnumfe --- number of continuum (FE) material
C             eltypenum --- number of element type in fe element library              
C             element --- number of FE element
C             r, s --- local coordinates of point of interest w.r.t. element

C     Outputs: displacement --- vector, length 2, of displacements at point of interest

C     Purpose: Get total displacement at point. Since there are no dislocations,
C     this is just the FE displacement (from assembly(mnumfe)%rhs)
      
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
      
      end function getTotalDispAtPointNoDisl                                                     
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

C     Purpose: Update FE nodal positions, displacements
C     (using FE solution from assembly(mnumfe)%(rhs) )

C     Notes: Technically the x-fixed and y-fixed checks are not necessary,
C     but using them ensures that the fixed points are truly fixed (otherwise,
C     small errors might accrue with many FE steps due to finite precision).

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
      
      do i = 1, size(feelements(mnumfe)%nodelist)
          node = feelements(mnumfe)%nodelist(i)
          posnundef = nodes%posn(1:2,node) - nodes%posn(4:5,node)
          nodetype = nodes%types(2,node)
          bcflag = nodes%types(3,node)
          xfixed = ((nodetype==2).or.(bcflag==1).or.(bcflag==3))
          yfixed = ((nodetype==2).or.(bcflag==2).or.(bcflag==3))
          if (.not.xfixed) then
              udisp = assembly(mnumfe)%rhs(2*i - 1)
              nodes%posn(1,node) = udisp + posnundef(1)
              nodes%posn(4,node) = udisp
          end if
          if (.not.yfixed) then
              vdisp = assembly(mnumfe)%rhs(2*i)
              nodes%posn(2,node) = vdisp + posnundef(2)
              nodes%posn(5,node) = vdisp
          end if
      end do    
      
      end subroutine updateFENodalPosn
************************************************************************
      end module   
      