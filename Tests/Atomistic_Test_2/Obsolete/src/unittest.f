      program unittest

      use mod_types, only: dp
      use mod_utils, only: prettyPrintIntMat, prettyPrintRealMat,
     &                      writeRealMat,writeIntMat
      use mod_math, only: rotateVec2d, rotateStress2d,
     &    piconst, invertmat2, getUnitNormalRHR, getDuplicates,
     &    getIntersectionTwoLines, getUniqueInts, getIdentityMatrix,
     &    linearInterp, searchSortedBinary, searchSortedBrute,
     &    searchSortedSpecial, isMultiple
      use mod_nodes, only: nodes, initNodeData, readNodeData,
     &                          processNodeData, writeNodeData
      use mod_misc, only: misc, initMiscData, readMiscData,
     &                          writeMiscData
      use mod_neighbors, only: neighbors, updateNeighbors,
     &                         initNeighborData, writeNeighborData,
     &                         genNeighList, readNeighborData,
     &                         processNeighborData, getXYAtomBounds,
     &                         checkDisp, updatePosSinceLastCheck,
     &                         getAtomsInBoxGroup
      use mod_materials, only: initMaterialData, materials,
     &    readMaterialData, processMaterialData, writeMaterialData
      use mod_potentials, only: initPotentialData, potentials,
     &                          readPotentialData, writePotentialData
      use mod_disl, only: initDislData, disl, readDislData,
     &    writeDislData, checkTooManyDisl, assignDislocation,
     &    addDislocation, deleteDislocation, imagedisl,
     &    readImageDislData, writeImageDislData, addImageDislocation
      use mod_interactions, only: initInteractionData, interactions,
     &  readInteractionData, writeInteractionData,
     &  processInteractionData
      use mod_poten_pair, only: getPotForceTable, getPotForcesAll
      use mod_integrate, only: loopVerlet, velocityVerletSubAll,
     &                         velocityVerletSub1, velocityVerletSub2
      use mod_groups, only: groups, initGroupData, processGroupData,
     &            readGroupData, writeGroupData, getMaskSortedIntersect,
     &            ngroups, genGroupMaskAll, genGroupMaskAtoms,
     &            genGroupMaskFENodes
      use mod_kdispfield, only: applyKDispIso
      use mod_disl_fields2, only: getDispAtPointSub,
     & getStressAtPointSub, adjustDxnDyn, getDispAtPoint,
     & getStressAtPoint, getDispAtPointHelper, getImageDispAtPoint,
     & getImageDispAtPointSub
      use mod_fe_el_2d, only: getK_2d, getGaussEdge_2d, getN_2d,
     &   getdN_2d, getJ_2d, getB_2d, getGauss_2d, getBAlt_2d, getBSub_2d
      use mod_fe_elements, only: readFEElementData,
     &   feelements, processFEElementData, processEdges, interfaceedges,
     &   processNodeLists, processMaterialList, fematerials,
     &   writeFEElementData, initFEElementData, processInterfaceEdges
      use mod_damping, only: getDampingForcesAll, initDampingData,
     &          writeDampingData,
     &          getDampingForce, readDampingData, damping
      use mod_sort, only: mergeSub, mergeSort, mergeSortCols,
     &                    mergeSubCols
      use mod_mesh_find, only: getLocalCoords, findinCPE4,
     &   getElGuessBrute, getNodeGuessBrute, checkSameSide, checkCPE3,
     &   checkCPE4, findInOneMat, findinOneMatAlt, findInAllWithGuess,
     &   findinCPE3, findinElement, checkElement, findInAllInitially,
     &   findInAllSub, findInOneMatInitially
      use mod_fe_main_2d, only: solveAll, solveOneMat, countEqns,
     &  assembleAndFactorAll, assembleKNormal, assembleLagrange,
     &  getDislForceRHS, getDislForceSub, getDispRHS, initAssembly,
     &  updateNodalPosn, getFEDispAtPoint, getFEStrainAtPoint, assembly
      use mod_disl_detect_pass, only: getDetectionElementStrain,
     &  getEForElement, getDetectionElementNorm, detectDislocationSub,
     &  updateDislSlip, findInterfaceIntersection,
     &  passAtomisticToContinuum, detection, passContinuumToAtomistic,
     &  addImageDislCrack, updateAtomsPassing, imposeDipoleDispOnAtoms,
     &  initDetectionData, readDetectionData, processDetectionData,
     &  writeDetectionData, processDetectNodeLists, processBAlt,
     &  processComparisons, getEcomp
      use mod_io, only: writeDump, initAll, writeAll, initAtomistic
     
      implicit none
      
C     integer :: element
C     integer :: compbest
C     real(dp) :: theta
C     real(dp) :: disp(2), dislpos(2)
C     
C     call initAtomistic()      
C     call initDetectionData(trim(misc%simname)//'_detection')
C     element = 3
C     theta = piconst
C     disp = [cos(theta),sin(theta)]
C     nodes%posn(1:2,7) = nodes%posn(1:2,7) + disp
C     nodes%posn(4:5,7) = nodes%posn(4:5,7) + disp
C     call detectDislocationSub(element,compbest,dislpos)
C     write(*,*) 'compbest', compbest
C     write(*,*) 'thetacomp', detection%thetacomp(compbest)
C     write(*,*) 'sgncomp', detection%sgncomp(compbest)
C     write(*,*) 'dislpos', dislpos

C     call initAtomistic()
C     call readDetectionData(trim(misc%simname)//'_detection')
C     call processDetectionData()
C     call processDetectNodeLists()
C     call processComparisons()
C     write(*,*) 'Sgn', detection%sgncomp
C     write(*,*) 'Theta', detection%thetacomp
C     call prettyPrintRealMat(detection%Ecompmat(:,:,1),'compnone')
C     call prettyPrintRealMat(detection%Ecompmat(:,:,2),'comp-60neg')
C     call prettyPrintRealMat(detection%Ecompmat(:,:,3),'comp-60')
C     call prettyPrintRealMat(detection%Ecompmat(:,:,4),'comp0neg')
C     call prettyPrintRealMat(detection%Ecompmat(:,:,5),'comp0')
C     call prettyPrintRealMat(detection%Ecompmat(:,:,6),'comp60neg') 
C     call prettyPrintRealMat(detection%Ecompmat(:,:,7),'comp60')

      real(dp) :: disp(2)
      real(dp) :: theta
      integer :: element
      real(dp) :: posndef(2,3)
      real(dp) :: dUdxmat(2,2)
      real(dp) :: E(2,2)

      call initAtomistic()      
      nodes%posn(1:2,5) = [1.0_dp,0.0_dp]
      nodes%posn(1:2,6) = [0.0_dp,0.0_dp]
      nodes%posn(1:2,7) = [0.5_dp,sqrt(3.0_dp)/2.0_dp]
      call initDetectionData(trim(misc%simname)//'_detection')
      nodes%posn(4:5,5:7) = 0.0_dp
      
      theta = -piconst/3.0_dp
      disp = [cos(theta),sin(theta)]
      write(*,*) 'disp', disp
      nodes%posn(1:2,6) = nodes%posn(1:2,6) + disp
      nodes%posn(4:5,6) = nodes%posn(4:5,6) + disp
      element = 3
      call getDetectionElementStrain(element,posndef,dUdxmat)
      call prettyPrintRealMat(dUdxmat,'dUdxmat')
      E = getEforElement(dUdxmat)
      call prettyPrintRealMat(E,'E')

      nodes%posn(1:2,5) = [1.0_dp,0.0_dp]
      nodes%posn(1:2,6) = [0.0_dp,0.0_dp]
      nodes%posn(1:2,7) = [0.5_dp,sqrt(3.0_dp)/2.0_dp]
      nodes%posn(4:5,5:7) = 0.0_dp
      
      theta = 0.0_dp
      disp = [cos(theta),sin(theta)]
      write(*,*) 'disp', disp
      nodes%posn(1:2,7) = nodes%posn(1:2,7) + disp
      nodes%posn(4:5,7) = nodes%posn(4:5,7) + disp
      element = 3
      call getDetectionElementStrain(element,posndef,dUdxmat)
      call prettyPrintRealMat(dUdxmat,'dUdxmat')
      E = getEforElement(dUdxmat)
      call prettyPrintRealMat(E,'E')
      
      end program