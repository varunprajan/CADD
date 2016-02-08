      module mod_disl_detect_pass

C     Purpose: Reads/writes/stores information related to dislocation
C     detection and passing. Also contains procedures for detecting
C     dislocations (using strains of detection elements) and passing dislocations
C     (both from the atomistic region to the continuum and vice versa).
      
      use mod_types, only: dp
      use mod_utils, only: writeIntMat, prettyPrintRealMat
      use mod_math, only: getIdentityMatrix,
     &        getIntersectionTwoLines, getUniqueInts, piconst
      use mod_nodes, only: nodes
      use mod_fe_el_2d, only: getBAlt_2d
      use mod_fe_elements, only: interfaceedges, fematerials
      use mod_groups, only: groups, tempgroupnum
      use mod_disl_fields2, only: getDispAtPoint, getImageDispAtPoint
      use mod_disl, only: addDislocation, addImageDislocation,
     &                    deleteDislocation
      use mod_mesh_find, only: findInAllWithGuess
      use mod_materials, only: materials, nmaterials
      use mod_integrate, only: loopVerlet
      use mod_neighbors, only: updateNeighbors, getAtomsInBoxGroup
      implicit none
 
      private
      public :: getDetectionElementStrain, getEForElement,
     &          getDetectionElementNorm, detectDislocationSub,
     &          updateDislSlip, findInterfaceIntersection,
     &          passAtomisticToContinuum, detection,
     &          passContinuumToAtomistic, addImageDislCrack,
     &          updateAtomsPassing, imposeDipoleDispOnAtoms,
     &          initDetectionData, readDetectionData,
     &          processDetectionData, writeDetectionData,
     &          processDetectNodeLists, processBAlt, processComparisons,
     &          getEcomp
      
      type detectiondata
C     (read-in)
      integer, allocatable :: connect(:,:)
      character(len=20) :: elname
      integer, allocatable :: interfaceedges(:,:)
      logical :: iscrackproblem
      integer, allocatable :: neighbors(:,:)
      real(dp) :: passdistance
      real(dp) :: timestepdamp
      real(dp) :: ycrack
C     (processed)
C     processDetectionData
      character(len=20) :: lattice
      real(dp) :: burgers
      integer :: neldof
      integer :: nelnodes
      integer :: nelements
C     processNodeLists
      integer, allocatable :: nodelist(:)
      integer, allocatable :: nodeinvlist(:)
      real(dp), allocatable :: prevslip(:,:,:)
C     processBalt
      real(dp), allocatable :: Balt(:,:,:)
C     processComparisons
      real(dp), allocatable :: Ecompmat(:,:,:)
      real(dp), allocatable :: thetacomp(:)
      integer, allocatable :: sgncomp(:)      
      end type
      
      type(detectiondata) :: detection

C     TODO: check if use of deformed/undeformed positions is correct
C     Problem: Dislocation is found in mesh using mesh_find, which uses
C     undeformed positions, but this module uses deformed positions to track crossings, etc....
C     so, dislocation may be "lost"

C     module variables (private)
      real(dp), parameter :: distfudge = 1.5_dp ! this parameter is not too important, see passAtomistictoContinuum
      real(dp), parameter :: distfudge2 = 3.0_dp ! see passAtomistictoContinuum. If dislocation path crosses interface at p,
                                                 !  dislocation is placed at p + distfudge2*b*m, where m is direction of travel
      real(dp), parameter :: boxfudge = 10.0_dp ! see updateAtomsPassing
      real(dp), parameter :: boxfudge2 = 20.0_dp ! see updateAtomsPassing
      integer, parameter :: passMDincr = 100 ! see passAtomistictoContinuum
      real(dp), parameter :: passgamma = 0.5_dp ! see passAtomistictoContinuum
      real(dp), parameter :: distfudge3 = 10.0_dp ! this parameter is not too important, see updateDislSlip
      
      contains
************************************************************************
      subroutine initDetectionData(detectionfile)

C     Subroutine: initDetectionData

C     Inputs: detectionfile --- filename where dislocation detection
C     and passing data is stored
C     (should be something like '[filepref]_detection')

C     Outputs: None

C     Purpose: Read, initialize data in "detection" structure, which holds
C     information about dislocation detection and passing, including
C     the detection band elements and edges, passing distances, etc.
      
      implicit none
      
C     input variables
      character(len=*) :: detectionfile
      
      call readDetectionData(detectionfile)
      call processDetectionData()
      call processBalt()
      call processDetectNodeLists()
      call processComparisons()
      
      end subroutine initDetectionData
************************************************************************
      subroutine readDetectionData(detectionfile)

C     Subroutine: readDetectionData

C     Inputs: detectionfile --- filename where dislocation detection
C     and passing data is stored
C     (should be something like '[filepref]_detection')

C     Outputs: None

C     Purpose: Read, detection data (detection elements, edges, etc.)
C     from file, initialize/allocate structures/arrays
           
      implicit none
      
C     input variables
      character(len=*) :: detectionfile
      
C     local variables
      integer :: iunit, j, k
      integer :: m, nrow
      integer :: temp
      
      open(newunit=iunit,file=detectionfile)
      
      read(iunit,*) m, nrow
      allocate(detection%connect(nrow,m))
      read(iunit,*) ((detection%connect(k,j),k=1,nrow),j=1,m)
      read(iunit,*) detection%elname
      read(iunit,*) m, nrow
      allocate(detection%interfaceedges(nrow,m))
      read(iunit,*) ((detection%interfaceedges(k,j),k=1,nrow),j=1,m)
      read(iunit,*) temp
C     use explicit integer -> logical conversion
      detection%iscrackproblem = (temp /= 0)
      read(iunit,*) m, nrow
      allocate(detection%neighbors(nrow,m))
      read(iunit,*) ((detection%neighbors(k,j),k=1,nrow),j=1,m)
      read(iunit,*) detection%passdistance
      read(iunit,*) detection%timestepdamp
      read(iunit,*) detection%ycrack
      
      close(iunit)
      
      end subroutine readDetectionData
************************************************************************
      subroutine processDetectionData()

C     Subroutine: processDetectionData

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize some basic properties of detection structure:
C     lattice name, burgers vector, nelnodes, neldof

C     Notes/TODO: No support for 3D yet...
      
      implicit none

C     local variables
      integer :: i

      detection%lattice = materials(1)%lattice
      detection%burgers = materials(1)%burgers
C     check to see that these are the same across materials...
C     (is this check strictly a good idea?)
      do i = 2, nmaterials
          if (detection%lattice/=materials(i)%lattice) then
              write(*,*) 'Dissimilar lattices'
              write(*,*) 'Unclear how dislocation passing should work'
              stop
          end if
          if (detection%burgers/=materials(i)%burgers) then
              write(*,*) 'Dissimilar magnitudes of burgers vectors'
              write(*,*) 'Unclear how dislocation passing should work'
              stop
          end if
      end do
      
      detection%nelnodes = size(detection%connect,1)
      if (detection%elname == 'CPE3') then
          if (detection%lattice/='hex') then ! only CPE3 is appropriate for hex lattice, I think
              write(*,*) 'Inappropriate lattice for this element type'
              write(*,*) 'Element name: ', detection%elname
              write(*,*) 'Lattice: ', detection%lattice
              stop
          end if
          detection%neldof = 6
      else
          write(*,*) 'Unknown element type'
          write(*,*) 'Currently defined elements:'
          write(*,*) '1) CPE3'
          stop
      end if
      
      detection%nelements = size(detection%connect,2)
      allocate(detection%prevslip(2,detection%nelnodes,
     &                            detection%nelements))
      detection%prevslip = 0.0_dp
      
      end subroutine processDetectionData
************************************************************************
      subroutine processBalt()

C     Subroutine: processBalt

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize array that stores shape function derivatives for
C     all detection band elements

C     Notes: This may be rather inefficient from the perspective of memory usage,
C     particularly in 3D and for many detection elements. On the other hand,
C     it's faster than calling getBalt every time...

      implicit none

C     local variables
      integer :: i, k
      integer :: ndetectelements
      integer :: node
      real(dp) :: posn(2,detection%nelnodes)
      real(dp) :: J(2,2)
      real(dp) :: r, s
      real(dp) :: Balt(4,detection%neldof)

      ndetectelements = size(detection%connect,2)
C     r, s, are irrelevant for triangles/tetrahedra
C     TODO: may need to fix for non-triangles/tetrahedra
      r = 0.0_dp
      s = 0.0_dp 
      allocate(detection%Balt(4,detection%neldof,ndetectelements))
      do k = 1, size(detection%connect,2)
          do i = 1, detection%nelnodes
              node = detection%connect(i,k)
              posn(:,i) = nodes%posn(1:2,node) - nodes%posn(4:5,node) ! undeformed positions
          end do
          call getBAlt_2d(posn,r,s,detection%neldof,detection%nelnodes,
     &                    detection%elname,Balt,J)
          detection%Balt(:,:,k) = Balt
      end do
      
      end subroutine processBalt
************************************************************************
      subroutine processDetectNodeLists()

C     Subroutine: processDetectNodeLists

C     Inputs: None

C     Outputs: None

C     Purpose: Determine nodes in detection elements, and "inverse" node list.
      
C     local variables
      integer :: counter

      allocate(detection%nodeinvlist(nodes%nnodes))
      call getUniqueInts(detection%connect,nodes%nnodes,
     &       counter,detection%nodelist,detection%nodeinvlist)
      
      end subroutine processDetectNodeLists
************************************************************************
      subroutine processComparisons()

C     Subroutine: processComparisons

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize array that stores strains for all possible
C     slips within element. This will be compared against the actual
C     strain in the element using the L2norm.

C     Notes/TODO: Can we distinguish between dislocation with b, theta
C     and dislocation with -b, theta + pi/2?

C     Notes/TODO: No 3D support...

      implicit none

C     local variables
      integer :: i, j, comp
      integer :: ncomparisons
      integer :: bsgn
      real(dp) :: theta
      real(dp) :: burgers, d
      real(dp) :: Ecomp(2,2)

      burgers = detection%burgers
      if (detection%lattice=='hex') then
          d = burgers*sqrt(3.0_dp)*0.5_dp ! slip plane spacing
          ncomparisons = 7
          allocate(detection%sgncomp(ncomparisons))
          allocate(detection%thetacomp(ncomparisons))
          allocate(detection%Ecompmat(2,2,ncomparisons))
          comp = 1
          detection%sgncomp(comp) = 0 ! fake value (no slip)
          detection%thetacomp(comp) = 0.0_dp ! fake value (no slip)
          detection%Ecompmat(:,:,comp) = 0.0_dp ! (no slip)
          do i = -1, 1
              theta = piconst/3.0_dp*i ! -60, 0, and 60 slip planes
              do j = -1, 1, 2
                  bsgn = j
                  Ecomp = getEcomp(bsgn,theta,burgers,d)
                  comp = comp + 1
                  detection%sgncomp(comp) = bsgn
                  detection%thetacomp(comp) = theta
                  detection%Ecompmat(:,:,comp) = Ecomp
              end do
          end do
      end if    
      
      end subroutine processComparisons
************************************************************************
      subroutine writeDetectionData(detectionfile)

C     Subroutine: writeDetectionData

C     Inputs: detectiontfile --- filename where where dislocation detection
C     and passing data is stored
C     (should be something like '[filepref]_detection')

C     Outputs: None

C     Purpose: Write dislocation detection and passing data to file (essentially
C     inverse of readDetectionData). Useful in creating "restart" file
      
      implicit none
      
C     input variables
      character(len=*) :: detectionfile
      
C     local variables
      integer :: iunit
      integer :: temp
      
      open(newunit=iunit,file=detectionfile)
      
      call writeIntMat(detection%connect,iunit,.true.)
      write(iunit,*) detection%elname
      call writeIntMat(detection%interfaceedges,iunit,.true.)
C     use explicit logical -> integer conversion
      if (detection%iscrackproblem) then
          temp = 1
      else
          temp = 0
      end if
      write(iunit,*) temp
      call writeIntMat(detection%neighbors,iunit,.true.)
      write(iunit,*) detection%passdistance
      write(iunit,*) detection%timestepdamp
      write(iunit,*) detection%ycrack
      
      close(iunit)
      
      end subroutine writeDetectionData
************************************************************************
      function getEcomp(bsgn,theta,burgers,d) result(Ecomp)

C     Function: getEcomp

C     Inputs: bsgn - sign of dislocation (+1 or -1)
C             theta - angle of dislocation
C             burgers - magnitude of burgers vector
C             d - slip plane spacing

C     Outputs: Ecomp - plastic strain matrix (Ep), 2 by 2 (in 2D), symmetric

C     Purpose: Generate plastic strain matrix, Ep (see Shilkrot et al., JMPS, 2004)
C     for dislocation slip, using dislocation attributes (b, slip plane spacing, etc.)
C     Should give result in Equation 19.

C     Notes/TODO: No 3d support...

      implicit none
      
C     input variables
      integer :: bsgn
      real(dp) :: theta
      real(dp) :: burgers
      real(dp) :: d
      
C     output variables
      real(dp) :: Ecomp(2,2)
      
C     local variables
      integer :: i, j
      real(dp) :: c, s
      real(dp) :: bvec(2), mvec(2)
      real(dp) :: dUdxmat(2,2)
      
      c = cos(theta)
      s = sin(theta)
      bvec = bsgn*burgers*[c,s]
      mvec = [-s,c] ! slip plane normal
      do i = 1, 2
          do j = 1, 2
C             see Equation 17, Shilkrot et al., JMPS, 2004
              dUdxmat(i,j) = bvec(i)*mvec(j)/d
          end do
      end do
      Ecomp = getEForElement(dUdxmat)
      
      end function getEcomp
************************************************************************
      function getEForElement(dUdxmat) result(E)

C     Function: getEforElement

C     Inputs: dUdxmat - matrix of displacement derivatives, du_i/dx_j = F - I

C     Outputs: Ecomp - plastic strain matrix (Ep), 2 by 2 (in 2D), symmetric

C     Purpose: Generate plastic strain matrix, Ep (see Shilkrot et al., JMPS, 2004)
C     for dislocation slip, using dislocation attributes (b, slip plane spacing, etc.)
C     Should give result in Equation 19.

C     Notes/TODO: No 3d support...

      implicit none
      
C     input variables
      real(dp) :: dUdxmat(2,2)
      
C     output variables
      real(dp) :: E(2,2)
      
C     local variables
      real(dp) :: idmat(2,2)
      real(dp) :: F(2,2)
      
      idmat = getIdentityMatrix(2)
      F = idmat + dUdxmat
      E = 0.5_dp*(matmul(transpose(F),F) - idmat)
      
      end function getEForElement
************************************************************************
      subroutine getDetectionElementStrain(element,posndef,dUdxmat)

C     Subroutine: getDetectionElementStrain

C     Inputs: element --- number of detection band element 
          
C     Outputs: posndef --- matrix, 2 by nelnodes, of deformed positions of element nodes
C              dUdxmat --- matrix, 2 by 2 of strains in detection band element

C     Purpose: Generates deformed positions and strains within detection band element

C     Notes/TODO: No 3d support...      

      implicit none
      
C     input variables
      integer :: element
      
C     output variables
      real(dp) :: dUdxmat(2,2)
      
C     local variables
      integer :: i
      real(dp) :: posndef(2,detection%nelnodes)
      integer :: node, nodeidx
      real(dp) :: disptot(2), disptilde(2)
      real(dp) :: disp(detection%neldof)
      real(dp) :: dUdx(4) 
      
      do i = 1, detection%nelnodes
          node = detection%connect(i,element)
          disptot = nodes%posn(4:5,node)
          nodeidx = detection%nodeinvlist(node)
          disptilde = detection%prevslip(:,i,element)
          disp(2*i-1:2*i) = disptot - disptilde
          posndef(:,i) = nodes%posn(1:2,node)
      end do
      
      dUdx = matmul(detection%Balt(:,:,element),disp)
C     assemble dUdx into matrix
      dUdxmat(1,:) = dUdx(1:2)
      dUdxmat(2,:) = dUdx(3:4)
      
      end subroutine getDetectionElementStrain
************************************************************************
      function getDetectionElementNorm(E,Ecomp) result(norm)

C     Subroutine: getDetectionElementNorm

C     Inputs: E - matrix of plastic strains (ndim by ndim)
C             Ecomp - comparison matrix of plastic strains (ndim by ndim)
          
C     Outputs: norm --- scalar, measuring the difference between E and Ecomp
C             (in this case, L2 norm *squared*)

C     Purpose: Measure difference between E and Ecomp

      implicit none
      
C     input variables
      real(dp) :: E(:,:)
      real(dp) :: Ecomp(:,:)
      
C     output variables
      real(dp) :: norm
      
      norm = sum((E - Ecomp)**2)    
      
      end function getDetectionElementNorm
************************************************************************
      subroutine detectAndPassDislocations()

      implicit none
     
C     local variables
      integer :: i
      integer :: compbest
      real(dp) :: dislpos(2), disldisp(2)
      real(dp) :: dislposold(2), dislposnew(2), dispfudge(2)
      integer :: bsgn
      real(dp) :: theta
      
      do i = 1, size(detection%connect,2)
          call detectDislocationSub(i,compbest,dislpos)
          if (compbest/=1) then ! we have detected a dislocation
              bsgn = detection%sgncomp(compbest)
              theta = detection%thetacomp(compbest)
C             assume dislocation is traveling in direction of positive b? (TODO: check!)
              disldisp = bsgn*[cos(theta),sin(theta)]
              call passAtomisticToContinuum(dislpos,disldisp,bsgn,theta)
              
C             Update "prevslip" structure, which stores slip
C             from all previous dislocations passing through element.
C             The positions below are "fake", and are designed so that
C             the dipole is far enough way from the detection element.
              dispfudge = disldisp*detection%burgers*distfudge3
              dislposold = dislpos - dispfudge
              dislposnew = dislpos + dispfudge
              call updateDislSlip(i,dislposold,dislposnew,bsgn,theta)
              call updateNeighborDislSlip(i,
     &                              dislposold,dislposnew,bsgn,theta)
          end if
      end do
      
      end subroutine detectAndPassDislocations
************************************************************************
      subroutine detectDislocationSub(element,compbest,dislpos)

      implicit none
      
C     input variables
      integer :: element 
      
C     output variables
      integer :: compbest
      real(dp) :: dislpos(2)
      
C     local variables
      integer :: i
      integer :: ncomparisons
      real(dp) :: dUdxmat(2,2)
      real(dp) :: E(2,2), Ecomp(2,2)
      real(dp) :: norm, normbest
      real(dp) :: posndef(2,detection%nelnodes)
      
      call getDetectionElementStrain(element,posndef,dUdxmat)
      E = getEforElement(dUdxmat)
      ncomparisons = size(detection%Ecompmat,3)
      normbest = huge(normbest)
      call prettyPrintRealMat(E,'E')
      do i = 1, ncomparisons
          Ecomp = detection%Ecompmat(:,:,i)
          call prettyPrintRealMat(Ecomp,'Ecomp')
          norm = getDetectionElementNorm(E,Ecomp)
          if (norm < normbest) then
              normbest = norm
              compbest = i
          end if    
      end do
      dislpos = sum(posndef,2)/size(posndef,2)
      
      end subroutine detectDislocationSub
************************************************************************
      subroutine updateDislSlip(element,dislposold,dislposnew,
     &                          bsgnnew,theta)

      implicit none
      
C     input variables
      integer :: element
      real(dp) :: dislposold(2), dislposnew(2)
      integer :: bsgnnew
      real(dp) :: theta
      
C     local variables
      integer :: i
      integer :: node, nodeidx
      real(dp) :: atompos(2)
      real(dp) :: disp(2)
      integer :: mnum
      
      do i = 1, size(detection%connect,1)
          node = detection%connect(i,element)
          nodeidx = detection%nodeinvlist(node)
          atompos = nodes%posn(1:2,node)
          mnum = nodes%types(1,node)
          disp = getImageDispAtPoint(atompos,dislposold,
     &                               theta,-bsgnnew,mnum)
          disp = disp + getImageDispAtPoint(atompos,dislposnew,
     &                                      theta,bsgnnew,mnum)
          detection%prevslip(:,i,element) =
     &                detection%prevslip(:,i,element) + disp
      end do    
      
      end subroutine updateDislSlip
************************************************************************
      subroutine updateNeighborDislSlip(element,dislposold,dislposnew,
     &                                  bsgnnew,theta)

      implicit none
      
C     input variables
      integer :: element
      real(dp) :: dislposold(2), dislposnew(2)
      integer :: bsgnnew
      real(dp) :: theta
     
C     local variables
      integer :: i, neighbor
     
      do i = 1, size(detection%neighbors,1)
          neighbor = detection%neighbors(i,element)
          if (neighbor == 0) then ! no neighbor
              return
          else
              call updateDislSlip(neighbor,dislposold,dislposnew,
     &                            bsgnnew,theta)
          end if
      end do
      
      end subroutine updateNeighborDislSlip
************************************************************************
      subroutine passAtomisticToContinuum(dislpos,disldisp,bsgn,theta)
      
C     Subroutine: passAtomisticToContinuum

C     Inputs: dislpos --- position (length 2) of detected dislocation
C             disldisp --- vector (length 2) indicating direction of movement of disl.
C             bsgn --- sign of dislocation (+1 or -1)
C             theta --- orientation of dislocation (radians)

C     Outputs: None

C     Purpose: Pass a detected dislocation from the atomistic region to the
C              continuum region

C     Notes: There are a lot of hard-coded constants here...can we fix this?

C     Notes/TODO: I assume that if we are dealing with a crack problem,
C     then the dislocation was emitted from the crack, which resides on
C     the plane y = ycrack. This is somewhat hackish.
      
      implicit none
      
C     input variables
      real(dp) :: dislpos(2)
      real(dp) :: disldisp(2)
      integer :: bsgn
      real(dp) :: theta
      
C     local variables
      real(dp) :: disldispnorm(2)
      real(dp) :: dislposnew(2)
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      integer :: mnumfe
      integer :: element
      real(dp) :: xd, yd
      real(dp) :: r, s
      logical :: badflip
      real(dp) :: burgers

C     place end of dislocation path well within continuum
C     tacit assumption that continuum is sufficiently "thick"
C     (so there aren't multiple intersections)
C     if this is an issue, findInterfaceIntersection should be modified
C     to find *nearest* intersection
      disldispnorm = disldisp/sqrt(disldisp(1)**2 + disldisp(2)**2)
      dislposnew = dislpos +
     &                    disldispnorm*detection%passdistance*distfudge
      call findInterfaceIntersection(interfaceedges%array,.true.,
     &                           dislpos,dislposnew,pint,isint,edgenum)
      mnumfe = interfaceedges%array(3,edgenum)
      element = interfaceedges%array(4,edgenum)
      if (.not.isint) then ! if no intersection, something went wrong
          write(*,*) 'Dislocation intersection with interface'
          write(*,*) ' could not be found.'
          write(*,*) 'Disl. line from', dislpos
          write(*,*) 'to', dislposnew
          stop
      end if
      
C     determine actual new position of dislocation
      burgers = detection%burgers
      dislposnew = pint + disldispnorm*distfudge2*burgers
      
C     search for it in the mesh, using our very good initial guess (determined
C     from the intersection)
      xd = dislposnew(1)
      yd = dislposnew(2)
      call findInAllWithGuess(xd,yd,mnumfe,element,r,s,badflip)
      
      if (badflip) then ! dislocation couldn't be placed in mesh
          write(*,*) 'Dislocation position in mesh'
          write(*,*) ' could not be found.'
          write(*,*) 'xp', xd, 'yp', yd
      end if
      
C     add discrete dislocation to DD array
      call addDislocation(mnumfe,element,xd,yd,theta,bsgn,r,s)
      
C     Add image dislocation with *opposite* sign to image array, if we have
C     a free surface. (Must cancel out dislocation in continuum.)
      if (detection%iscrackproblem) then
          call addImageDislCrack(detection%ycrack,dislpos,dislposnew,
     &                           theta,-bsgn,mnumfe)
      end if
      
C     update atom positions, including imposing dipole displacements
C     and minimizing atomic positions near core
      call updateAtomsPassing(dislpos,dislposnew,
     &                        burgers,bsgn,theta,.true.)
      
      end subroutine passAtomisticToContinuum
************************************************************************
      subroutine passContinuumToAtomistic(dislpos,pint,bsgn,theta,
     &                                    dislnum,mnumfe)
      
C     Subroutine: passContinuumToAtomistic

C     Inputs: dislpos --- position (length 2) of discrete dislocation (before movement)
C             pint --- position (length 2) of intersection of dislocation path
C                      with interface between atomistic/continuum regions
C             bsgn --- sign of dislocation (+1 or -1)
C             theta --- orientation of dislocation (radians)
C             dislnum --- number of dislocation
C             mnumfe --- continuum material that discrete dislocation used to belong to

C     Outputs: None

C     Purpose: Pass a discrete dislocation from the continuum region to the
C              atomistic region

C     Notes: There are a lot of hard-coded constants here...can we fix this?

C     Notes/TODO: I assume that if we are dealing with a crack problem,
C     then the dislocation was emitted from the crack, which resides on
C     the plane y = ycrack, so this is where an image dislocation goes.
C     This is somewhat hackish.
      
      implicit none
      
C     input variables
      real(dp) :: dislpos(2)
      real(dp) :: pint(2)
      integer :: bsgn
      real(dp) :: theta
      integer :: dislnum
      integer :: mnumfe
      
C     local variables
      real(dp) :: disldisp(2), disldispnorm(2)
      real(dp) :: dislposnew(2)
      real(dp) :: pintnew(2)
      real(dp) :: burgers
      logical :: isint
      integer :: edgenum
      integer :: element
    
C     determine intersection of disl. path with detection band
      disldisp = pint - dislpos
      disldispnorm = disldisp/sqrt(disldisp(1)**2 + disldispnorm(2)**2)
      dislposnew = pint + disldispnorm*detection%passdistance*distfudge
      call findInterfaceIntersection(detection%interfaceedges,.false., ! use deformed positions of interfaceedges
     &                   pint,dislposnew,pintnew,isint,edgenum)
      element = detection%interfaceedges(3,edgenum)
      
C     place dislocation just inside detection band
      burgers = detection%burgers
      dislposnew = pintnew + disldispnorm*distfudge2*burgers
      
C     update prev slip array
      call updateDislSlip(element,pint,dislposnew,bsgn,theta)
      call updateNeighborDislSlip(element,pint,dislposnew,bsgn,theta)
      
C     delete discrete dislocation from DD array
      call deleteDislocation(mnumfe,dislnum)
      
C     Add image dislocation with *same* sign to image array, if we have
C     a free surface. (Sign is same because dislocation is leaving continuum,
C     which is equivalent to opposite sign dislocation being added to continuum,
C     and image must cancel this out.)
      if (detection%iscrackproblem) then
          call addImageDislCrack(detection%ycrack,dislpos,dislposnew,
     &                           theta,bsgn,mnumfe)
      end if
      
C     update atom positions, including imposing dipole displacements
C     and minimizing atomic positions near core
      call updateAtomsPassing(dislpos,dislposnew,
     &                        burgers,bsgn,theta,.false.)
      
      end subroutine passContinuumToAtomistic
************************************************************************
      subroutine addImageDislCrack(ycrack,dislpos,dislposnew,
     &                             theta,bsgn,mnumfe)
     
C     input variables
      real(dp) :: ycrack
      real(dp) :: dislpos(2), dislposnew(2)
      real(dp) :: theta
      integer :: bsgn
      integer :: mnumfe
      
C     local variables
      real(dp) :: s
      real(dp) :: imagepos(2)
     
C     parametrize line connecting dislpos and dislposnew using s, solve
C     for intersection with y = ycrack
      s = (ycrack-dislposnew(2))/(dislpos(2)-dislposnew(2))
      imagepos = dislpos*s + dislposnew*(1.0_dp-s)
      call addImageDislocation(imagepos(1),imagepos(2),
     &                         theta,bsgn,mnumfe)
      
      end subroutine addImageDislCrack
************************************************************************
      subroutine updateAtomsPassing(dislposold,dislposnew,
     &                              burgers,bsgn,theta,tocontinuum)

C     Subroutine: updateAtomsPassing

C     Inputs: dislpos --- old position (length 2) of dislocation
C             disldisp --- vector (length 2) indicating direction of movement of disl.
C             bsgn --- sign of dislocation (+1 or -1)
C             theta --- orientation of dislocation (radians)
C             tocontinuum --- flag if we're passing to the continuum

C     Outputs: None

C     Purpose: Update atomic positions as a result of dislocation being passed
C     (in either direction). Helper for passAtomistictoContinuum and
C     passContinuumtoAtomistic

C     Notes: There are a lot of hard-coded constants here...can we fix this?

      implicit none
      
C     input variables
      real(dp) :: dislposold(2)
      real(dp) :: dislposnew(2)
      real(dp) :: burgers
      integer :: bsgn
      real(dp) :: theta
      logical :: tocontinuum
      
C     local variables
      real(dp) :: xmin, xmax, ymin, ymax
      logical :: dampflag
      real(dp) :: corepos(2)

C     create group of atoms for dislocation dipole displacements
C     we only need to get displacements right for atoms sufficiently close to dipole...
C     (since dipole field decays rapidly in space)
C     so, we define a box around dipole with a "fudge" of boxfudge*burgers
C     (this is a bit hackish),
      xmin = min(dislposnew(1),dislposold(1)) - boxfudge*burgers
      xmax = max(dislposnew(1),dislposold(1)) + boxfudge*burgers
      ymin = min(dislposnew(2),dislposold(2)) - boxfudge*burgers
      ymax = max(dislposnew(2),dislposold(2)) + boxfudge*burgers      
      call getAtomsInBoxGroup(xmin,xmax,ymin,ymax,tempgroupnum)
      
C     add dislocation dipole displacements (subtract disp. for dislposold; add for dislposnew)
      call imposeDipoleDispOnAtoms(dislposnew,dislposold,
     &                             bsgn,theta,tempgroupnum)

C     regenerate neighbor list
      call updateNeighbors()

C     figure out core position, based on flag indicating which way dislocation is going
      if (tocontinuum) then
          corepos = dislposold
      else
          corepos = dislposnew
      end if
C     create box around core for atom group
      xmin = corepos(1) - boxfudge2*burgers
      xmax = corepos(1) + boxfudge2*burgers
      ymin = corepos(2) - boxfudge2*burgers
      ymax = corepos(2) + boxfudge2*burgers
      call getAtomsInBoxGroup(xmin,xmax,ymin,ymax,tempgroupnum)

C     run damped dynamics near (former) atomistic core
      dampflag = .true.
      call loopVerlet(passMDincr,detection%timestepdamp,
     &                tempgroupnum,tempgroupnum,dampflag,passgamma)
      
      end subroutine updateAtomsPassing
************************************************************************
      subroutine findInterfaceIntersection(edgesarray,undeformed,
     &                                p0,p1,pint,isint,edgenum)

C     Subroutine: findInterfaceIntersection
      
C     Inputs: edgesarray --- array (2 by N), where each column consists of two nodes defining an edge
C             p0, p1 --- coordinates (vector, length 2) of points defining line
      
C     Outputs: pint --- coordinates of intersection point (if it exists) (vector, length 2)
C              isint --- boolean indicating whether intersection exists
C              mnumfe --- fe material that interface edge belongs to
C              element --- element that interface edge belongs to
C              undeformed --- (logical) flag indicating whether undeformed positions are to be used

C     Purpose: Find where a line between p0 and p1 (i.e. a path of a dislocation) intersects
C     an edge in the array edgesarray. If search is successful,
C     return pint (intersection point) and mnumfe (number of continuum material)
      
C     Notes: Assumes multiple intersections are not possible...

      implicit none
      
C     input variables
      integer :: edgesarray(:,:)
      logical :: undeformed
      real(dp) :: p0(2), p1(2)
      
C     output variables
      real(dp) :: pint(2)
      logical :: isint
      integer :: edgenum
      
C     local variables
      integer :: i
      integer :: node1, node2
      real(dp) :: posn1(2), posn2(2)

C     loop over interface edges
      do i = 1, size(edgesarray,2)
          node1 = edgesarray(1,i)
          node2 = edgesarray(2,i)
          posn1 = nodes%posn(1:2,node1)
          posn2 = nodes%posn(1:2,node2)
          if (undeformed) then
              posn1 = posn1 - nodes%posn(4:5,node1)
              posn2 = posn2 - nodes%posn(4:5,node2)             
          end if
          call getIntersectionTwoLines(p0,p1,posn1,posn2,pint,isint)
          if (isint) then
              edgenum = i
              return
          end if    
      end do
      
C     if we have gotten here, no intersection was found
      isint = .false.

      end subroutine findInterfaceIntersection
************************************************************************
      subroutine imposeDipoleDispOnAtoms(dislpos1,dislpos2,
     &                                            bsgn1,theta,gnum)
     
C     Subroutine: imposeDipoleDispOnAtoms

C     Inputs: dislpos1 --- position of first dislocation (vector, length 2)
C             dislpos2 --- position of second dislocation (vector, length 2)
C             bsgn1 --- sign of burgers vector of first dislocation (+1 or -1)
C                       (sign of second dislocation is -bsgn1, since we assume a dipole)
C             theta --- orientation of dislocation
C             gnum --- number of group for which to apply displacements to atoms (see mod_groups)
      
C     Outputs: None

C     Purpose: Apply displacements of a dislocation dipole to all atoms in a group (see mod_groups)
C              *including* pad atoms
      
C     input variables
      real(dp) :: dislpos1(2), dislpos2(2)
      integer :: bsgn1
      real(dp) :: theta
      integer :: gnum
      
C     local variables
      integer :: i
      integer :: node
      integer :: mnum
      real(dp) :: atompos(2)
      real(dp) :: disp(2)
      
      do i = 1, nodes%natoms
      if (groups(gnum)%maskatoms(i)) then
          node = nodes%atomlist(i)
          atompos = nodes%posn(1:2,node)
          mnum = nodes%types(1,node) ! this tacitly assumes no elastic mismatch in the box
          disp = getDispAtPoint(atompos,dislpos1,theta,bsgn1,mnum)
          disp = disp + getDispAtPoint(atompos,dislpos2,
     &                                 theta,-bsgn1,mnum)
          nodes%posn(1:2,node) = nodes%posn(1:2,node) + disp
          nodes%posn(4:5,node) = nodes%posn(4:5,node) + disp
      end if
      end do
      
      end subroutine imposeDipoleDispOnAtoms
************************************************************************     
      
      end module