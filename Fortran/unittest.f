      program unittest

      use mod_types, only: dp
      use mod_utils, only: prettyPrintMat, writeMat, prettyPrintVec,
     &                     writeMatTranspose, readMatTransposeSize
      use mod_math, only: rotateVec2d, rotateStress2d, linspace,
     &    piconst, invertmat2, getUnitNormalRHR, getDuplicates,
     &    getIntersectionTwoLines, getUniqueInts, normalizeVec,
     &    linearInterp, searchSortedBinary, searchSortedBrute,
     &    searchSortedSpecial, isMultiple, projectVec, getPerpDistance,
     &    sameSign, tolconst, findPointBetween, intToLogical,
     &    getCircumradiusSqForTriangle, logicalToInt, getDeterminant2,
     &    nearestMultiple
      use mod_nodes, only: nodes, initNodeData, readNodeData,
     &  processNodeData, writeNodeData, getXYAtomBounds, getXYBounds,
     &  getXYAtomBoundsUndef, getXYAtomBoundsDef
      use mod_misc, only: misc, initMiscData, readMiscData,
     &                    writeMiscData, updateMiscIncrementCurr
      use mod_neighbors, only: neighbors, updateNeighborsCheck,
     & initNeighborData, writeNeighborData, genNeighList,
     & readNeighborData, processNeighborData, checkDisp, atombox,
     & updatePosSinceLastCheck, getAtomsInBoxGroupTemp, genAtomBins,
     & updateNeighborsNoCheck, updateNeighIncrementCurr, getAtomBin,
     & genAtomBinsUndeformed, genAtomBinsCurrent
      use mod_materials, only: initMaterialData, materials,
     &    readMaterialData, processMaterialData, writeMaterialData,
     &    getNucleationLength, nmaterials, getMuNuApprox
      use mod_potentials, only: initPotentialData, potentials,
     &                          readPotentialData, writePotentialData
      use mod_disl_misc, only: initDislMiscData, readDislMiscData,
     &  writeDislMiscData, dislmisc
      use mod_disl_try, only: initDislData, disl, readDislData,
     &    writeDislData, checkTooManyDisl,
     &    addDislocation, deleteDislocation,
     &    obstacles, sources, readDislObsData, assignDislLocalPos,
     &    readDislSourceData, processDislSourceData,
     &    initDislObsData, initDislSourceData, processDislObsData,
     &    initDislSortedPlanes, initObsSortedPlanes,
     &    sortPlanes, sortPlane, assignDislSortedPlanes,
     &    assignObsSortedPlanes, assignSourcesSortedPlanes,
     &    sortObsPlanes, initSourcesSortedPlanes, countActiveDisl,
     &    sortDislPlanes, sortedplanedata, addDislocationSub,
     &    addObjSub, processDislData, findEmptyDislSlot,
     &    deleteDislocationSub2, checkTooManyObj, zeroDislDisp,
     &    writeDislObsData, writeDislSourceData, setupSources,
     &    assignSourcesLocalPos, assignNucleationLength,
     &    zeroObstacles, activateDislocations,
     &    assignObsLocalPos, deleteDislocationSub, sourcet
      use mod_slip_sys, only: initSlipSysData, readSlipSysData,
     &    processSlipSysData, getSlipPlane, writeSlipSysData,
     &    slipsys, invResolveDisp, resolveStress
      use mod_interactions, only: initInteractionData, interactions,
     &  readInteractionData, writeInteractionData,
     &  processInteractionData
      use mod_poten_pair_table, only: getPotForceTable,
     7   getPotForcesAllTable
      use mod_integrate, only: loopVerlet, velocityVerletSubAll,
     &                         velocityVerletSub1, velocityVerletSub2
      use mod_groups, only: groups, initGroupData, processGroupData,
     &            readGroupData, writeGroupData, getMaskSortedIntersect,
     &            ngroups, genGroupMaskAll, genGroupMaskAtoms,
     &            genGroupMaskFENodes, getGroupNum
      use mod_kdispfield, only: applyKDispIso
      use mod_disl_fields2, only: getDispAtPointSub,
     & getStressAtPointSub, adjustDxnDyn, getDispAtPoint,
     & getStressAtPoint
      use mod_fe_el_2d, only: getK_2d, getElTypeNum, getB_2d, 
     &    initFELibrary, felib, getBSub_2d, getJ_2d, findInCPE4,
     &          findInCPE3, checkCPE4, checkCPE3, getN_CPE3_2d,
     &          getN_CPE4_2d, getdN_CPE3_2d, getdN_CPE4_2d, getBAlt_2d,
     &          getGauss_CPE4_2d, getGauss_CPE3_2d
      use mod_fe_elements, only: readFEElementData,
     &   feelements, processFEElementData, processEdges, interfaceedges,
     &   processNodeLists, processMaterialList, fematerials,
     &   writeFEElementData, initFEElementData, processInterfaceEdges
      use mod_damping, only: getDampingForcesAll, initDampingData,
     &          writeDampingData, writeDampingDataSub, normaldamping,
     &          getDampingForce, readDampingData, readDampingDataSub
      use mod_sort, only: mergeSub, mergeSort, mergeSortCols,
     &  mergeSubCols, insertionSortPlane, mergeSortReal, mergeSubReal,
     &  mergeSortColsReal, mergeSubColsReal
      use mod_mesh_find, only: getLocalCoords, findInOneMatSub,
     &   getElGuessBrute, getNodeGuessBrute,
     &   findInOneMat, findinOneMatSub, findInAllWithGuess,
     &   findInAllSub, findInOneMatInitially
      use mod_fe_main_2d, only: solveAll, solveOneMat, countEqns,
     &  assembleAndFactorAll, assembleKNormal, assembleLagrange,
     &  getDislForceRHS, getDislForceSub, getDispRHS, initAssembly,
     &  updateFENodalPosn, getFEDispAtPoint,
     &  assembly
      use mod_fe_main_2d_no_disl, only: getFEStressAtPoint,
     &  getFEStrainAtPoint
      use mod_disl_detect_pass, only: initDetectionData,
     &  readDetectionData, processDetectionData, BOXWIDTHNORM,
     &  writeDetectionData, findDetectionNodes, placeDetectionSub,
     &  getDislBranchCut, insideAnnulus, getPaddedParamsAnnulus,
     &  insideRectAnnulus, detection, errorInterface,
     &  getDislPropsFromBurgersVec, detectAndPassDislocations,
     &  passContinuumToAtomistic, passAtomistictoContinuum,
     &  updateAtomsPassing, assignDetectionBand, placeInsideDetection,
     &  imposeDipoleDispOnAtoms, findInterfaceIntersectionDeformed,
     &  getPaddedParamsRectAnnulus, placeOutsideDetection
      use mod_io, only: initCADD, writeRestartCADD, initAtomistic,
     & initSimulation, initDD, initFE, writeRestart_ptr, writeDump_ptr,
     & initGeneralChunk, initAtomisticChunk, initFEChunk,
     & getRestartPrefSuff, writeRestartDD, writeRestartCADD,
     & writeRestartCADDNoDisl, writeRestartFE, writeRestartAtomistic
      use mod_delaunay, only: dtris2, genDelaunay, regenDelaunay,
     & identifyLargeTri, getTriCenter, getTriNodes, delaunaydata,
     & genBadTriangles, setDelaunayPos, setDelaunayPosDef,
     & setDelaunayPosUndef
      use mod_disl_ident_simple, only: initDislIdentData,
     & computeCircuits, identsimple, getIdealVector
      use mod_dd_main, only: runDDStep,
     & updateDislocations, enforceObstacles, isObstacleBetween,
     & getResolvedStressOnObstacle, updateDislRelPos,
     & insertionSortPlaneWithCrossing, annihilateDislocations,
     & annihilateDislocationsSub, updateDislPos, updateDislPosSub,
     & updateSources, getResolvedStressOnSource, updateSource,
     & createDipole, updateDislPosSub_ptr
      use mod_dd_integrate, only: dispFromPKOneMat, dispFromPK,
     & dispFromPKOneMatVCorr, velFromTau, velFromTauCorr,
     & getResolvedStressOnDisl, dispFromVel
      use mod_disl_escaped, only: escapeddisl, initEscapedDislData,
     &  readEscapedDislData, box, getMaxLen,
     &  processEscapedDislData, writeEscapedDislData, 
     &  addEscapedDislocation, getEscapedPos, getEscapedRegion,
     &  getEscapedRegionCrack, getEscapedDispAtPointAll,
     &  getEscapedDispAtPoint, getEscapedDispAtPointSub,
     &  getEscapedDispAtPointCrack, getEscapedDispAtPointCrackSub
      use mod_disl_ghost, only: ghostdisl, initGhostDislData,
     &  readGhostDislData, writeGhostDislData, addGhostDislocation
      use mod_pad_atoms, only: padatoms, initPad, assignPad, updatePad
      use mod_fe_main_2d_assign, only: assignFE,
     & getTotalDispAtPoint_ptr, solveAll_ptr, assembleAndFactorAll_ptr,
     & updateFENodalPosnAll_ptr
      use mod_compute, only: readComputeData, compute, initComputeData,
     & writeComputeData, getCentroAtom, getCentroAtoms
      use mod_dump, only: getDumpFilename, writeDumpCADD,
     & writeDumpCADDNoDisl, writeDumpFE, writeDumpDD,
     & writeDumpAtomistic, writeDumpNodesElementsChunkSub,
     & writeDumpDDChunk, writeDumpNodesDefChunk, writeDumpNodesDef,
     & writeDumpNodesUndef, writeDumpNodesDisp, writeDumpNodesTypes,
     & writeDumpFEElements, writeDumpDisl, writeDumpSources,
     & writeDumpObstacles, writeDumpCompute, writeDumpCentro,
     & writeDumpNodesUndefChunk, writeDumpNodesDefElementsChunk,
     & writeDumpNodesUndefElementsChunk
      use mod_find_crack_atomistic, only: initAtomFindCrackData,
     & readAtomFindCrackData, writeAtomFindCrackData, isTriangleOnEdge,
     & processAtomFindCrackData, crackfinding, findCrack
      use mod_moving_mesh_crack_cadd, only: moveMeshCrack, 
     &  interpFromFEPoint, interpFromAtomPoint, readCADDMovingMeshData,
     &  writeCADDMovingMeshData, passDislocationsBefore, shiftPosnDisp,
     &  passDislocationsAfter, passDislocationsSub, shiftPosnDisp,
     &  interpFromFEPoint, interpFromAtomPoint, checkTriangle,
     &  locateClosestAtomUndef, shiftDDPos, shiftSlipPlanes,
     &  updateDislPosMovingMesh, updateDislPosMovingMeshSub,
     &  updateSourcePosMovingMesh, updateObstaclePosMovingMesh,
     &  initCADDMovingMeshData, movingmesh, processCADDMovingMeshData,
     &  getNewDetectionBandAfter, getNewDetectionBandBefore,
     &  atomicDispInterpolation
     
      implicit none
      
      real(dp) :: orig(2), a1(2), a2(2)
      real(dp) :: fac
      integer :: i
      real(dp) :: y, ydisp
      logical :: failure
      real(dp) :: disp(2), posnundef(2)

      nodes%nnodes = 9
      allocate(nodes%posn(7,nodes%nnodes))
      nodes%posn = 0.0_dp
      orig = [0.0_dp,0.0_dp]
      a1 = [1.0_dp,0.0_dp]
      a2 = [-0.5_dp,0.5_dp*sqrt(3.0_dp)]
      nodes%posn(1:2,3) = orig
      nodes%posn(1:2,4) = orig + a1
      nodes%posn(1:2,5) = orig - a1
      nodes%posn(1:2,6) = orig + a2
      nodes%posn(1:2,7) = orig + a1 + a2
      nodes%posn(1:2,8) = orig - a2
      nodes%posn(1:2,9) = orig - a1 - a2
      
      allocate(nodes%types(3,nodes%nnodes))
      nodes%types(:,1) = [1,0,0]
      nodes%types(:,2) = [1,0,0]
      do i = 3, nodes%nnodes
          nodes%types(:,i) = [1,1,0]
      end do
      
      call processNodeData()
      
C     stretch in y-direction
      fac = 0.2_dp
      do i = 1, size(nodes%posn,1)
          y = nodes%posn(2,i)
          ydisp = fac*y
          nodes%posn(5,i) = nodes%posn(5,i) + ydisp
          nodes%posn(2,i) = nodes%posn(2,i) + ydisp
      end do
      
      neighbors%nmaxbin = 30
      neighbors%nmaxneigh = 30
      neighbors%rneigh = 0.5_dp
      neighbors%rneighsq = 0.25_dp
      neighbors%rhomax = 1.0_dp
      allocate(neighbors%binlist(2,nodes%natoms)) 
      
      call genAtomBinsUndeformed()
      movingmesh%delaunay%nodenums = nodes%atomlist
      call setDelaunayPosUndef(movingmesh%delaunay)
      call genDelaunay(movingmesh%delaunay)
      
      call initFELibrary()
      movingmesh%eltypenum = getElTypeNum('CPE3') ! triangle
      
      posnundef = [0.9_dp,0.8_dp]
      call atomicDispInterpolation(posnundef,failure,disp)
      write(*,*) 'Failed?', failure
      
      
            
      
      
      end program