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
     &  processNodeData, writeNodeData, getXYAtomBounds, getXYBounds
      use mod_misc, only: misc, initMiscData, readMiscData,
     &                    writeMiscData, updateMiscIncrementCurr
      use mod_neighbors, only: neighbors, updateNeighborsCheck,
     & initNeighborData, writeNeighborData, genNeighList,
     & readNeighborData, processNeighborData, checkDisp,
     & updatePosSinceLastCheck, getAtomsInBoxGroupTemp, genAtomBins,
     & updateNeighborsNoCheck, updateNeighIncrementCurr
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
     &  writeDetectionData, assignDetectionPoints, placeDetectionSub,
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
      use mod_delaunay, only: dtris2, delaunay, genDelaunay,
     & identifyLargeTri, getTriCenter, getTriNodes
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
     
      implicit none

      integer :: mnumfe, iplane, isys
      real(dp) :: relpos, pt(2), disp(2)

      
      call initSimulation('simple2_dd','dd')

      pt = [1.0_dp,2.0_dp]
      isys = 1
      mnumfe = 1
      disp = [3.0_dp,0.0_dp]
      call getSlipPlane(pt,mnumfe,isys,iplane,relpos)
      write(*,*) 'relpos', relpos
      write(*,*) 'iplane', iplane
      slipsys(mnumfe)%origin(:,isys) = slipsys(mnumfe)%origin(:,isys)
     &                               - disp
      pt = pt - disp
      call getSlipPlane(pt,mnumfe,isys,iplane,relpos)
      write(*,*) 'relpos', relpos
      write(*,*) 'iplane', iplane

      

      
      end program