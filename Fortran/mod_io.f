      module mod_io
      
C     Purpose: Contains routines for reading/initializing and writing data
C     from/to files. Can be used to write "dump" or "restart" files.
      
      use mod_types, only: dp
      use mod_nodes, only: nodes, initNodeData, writeNodeData
      use mod_misc, only: misc, initMiscData, writeMiscData
      use mod_damping, only: initDampingData, writeDampingData
      use mod_fe_elements, only: initFEElementData, writeFEElementData
      use mod_materials, only: initMaterialData, writeMaterialData
      use mod_potentials, only: initPotentialData, writePotentialData
      use mod_poten_assign, only: initPotentialStyle
      use mod_interactions, only: initInteractionData,
     &                            writeInteractionData
      use mod_neighbors, only: initNeighborData, writeNeighborData
      use mod_disl_misc, only: initDislMiscData, writeDislMiscData
      use mod_disl_try, only: initDislData, writeDislData, disl,
     & initDislSourceData, writeDislSourceData, 
     & initDislObsData, writeDislObsData, sources, obstacles
      use mod_disl_ghost, only: initGhostDislData, writeGhostDislData
      use mod_disl_escaped, only: initEscapedDislData,
     &                            writeEscapedDislData
      use mod_slip_sys, only: initSlipSysData, writeSlipSysData
      use mod_disl_detect_pass, only: initDetectionData,
     &                                writeDetectionData
      use mod_groups, only: initGroupData, writeGroupData
      use mod_math, only: isMultiple
      use mod_fe_el_2d, only: initFELibrary
      use mod_disl_ident_simple, only: initDislIdentData
      use mod_dump, only: writeDumpCADD, writeDumpCADDNoDisl,
     &  writeDumpAtomistic, writeDumpFE, writeDumpDD
      use mod_pad_atoms, only: initPad
      use mod_fe_main_2d_assign, only: assignFE, initAssembly_ptr,
     &  assembleAndFactorAll_ptr
      use mod_compute, only: initComputeData, writeComputeData
      use mod_dd_main, only: assignDD
      
      implicit none
      
      private :: suffixnodes, suffixmaterials, suffixpotentials,
     &           suffixmisc, suffixdisl, suffixghostdisl, 
     &           suffixinteractions, suffixneighbors, suffixdamping,
     &           suffixfeelements, suffixescapeddisl
      public :: initCADD, writeRestartCADD, initAtomistic,
     & initGeneralChunk, initAtomisticChunk, initFEChunk, initDDChunk,
     & initSimulation, assignInitAndWrite,
     & writeOutput, initCADDNoDisl, initFE, initDD, initCADDChunk,
     & initCADDNoDislChunk, getRestartPrefSuff,
     & writeRestartCADDNoDisl, writeRestartFE, writeRestartDD,
     & writeRestartAtomistic, writeRestartGeneralChunk,
     & writeRestartFEChunk, writeRestartDDChunk, writeRestartCADDChunk
      
      procedure(Dummy), pointer :: init_ptr => NULL()
      procedure(Dummy2), pointer :: writeRestart_ptr => NULL()
      procedure(Dummy3), pointer :: writeDump_ptr => NULL()
      
C     module variables
      character(len=*), parameter :: suffixnodes = '_nodes'
      character(len=*), parameter :: suffixmaterials = '_materials'
      character(len=*), parameter :: suffixpotentials = '_potentials'
      character(len=*), parameter :: suffixmisc = '_misc'
      character(len=*), parameter :: suffixslipsys = '_slipsys'
      character(len=*), parameter :: suffixdislmisc = '_dislmisc'
      character(len=*), parameter :: suffixdisl = '_disl'
      character(len=*), parameter :: suffixghostdisl = '_ghostdisl'
      character(len=*), parameter :: suffixescapeddisl = '_escapeddisl'
      character(len=*), parameter :: suffixdislsource = '_sources'
      character(len=*), parameter :: suffixdislobs = '_obstacles'
      character(len=*), parameter :: suffixinteractions ='_interactions'
      character(len=*), parameter :: suffixneighbors = '_neighbors'
      character(len=*), parameter :: suffixdamping = '_damping'
      character(len=*), parameter :: suffixfeelements = '_feelements'
      character(len=*), parameter :: suffixdetection = '_detection'
      character(len=*), parameter :: suffixgroups = '_groups'
      character(len=*), parameter :: suffixcompute = '_compute'
      
      contains
************************************************************************
      subroutine initSimulation(simname,simtype)
      
C     Subroutine: initSimulation

C     Inputs: simname --- name of simulation
C             simtype --- simulation type: atomistic, fe, dd, cadd_nodisl, or cadd

C     Outputs: None

C     Purpose: Initialize simulation; assign pointers based on simulation type,
C     and call initialization routine (e.g., initCADD) based on simulation type

      implicit none

C     input variables
      character(len=*) :: simname
      character(len=*) :: simtype

      misc%simname = simname
      misc%simtype = simtype
      call assignInitAndWrite(simtype) ! assign pointers for dump, restart, initialization, based on simulation type
      call init_ptr() ! call the initialization function to initialize the simulation
      
      end subroutine initSimulation
************************************************************************
      subroutine assignInitAndWrite(simtype)
      
C     Subroutine: initSimulation

C     Inputs: simtype --- simulation type: atomistic, fe, dd, cadd_nodisl, or cadd

C     Outputs: None

C     Purpose: Assign pointers based on simulation type (initialization, writing dump and restart files)

      implicit none

C     input variables
      character(len=*) :: simtype
      
      select case (trim(simtype))
          case ('cadd')
              init_ptr => initCADD
              writeDump_ptr => writeDumpCADD
              writeRestart_ptr => writeRestartCADD
          case ('cadd_nodisl')
              init_ptr => initCADDNoDisl
              writeDump_ptr => writeDumpCADDNoDisl
              writeRestart_ptr => writeRestartCADDNoDisl
          case ('atomistic')
              init_ptr => initAtomistic
              writeDump_ptr => writeDumpAtomistic
              writeRestart_ptr => writeRestartAtomistic
          case ('fe')
              init_ptr => initFE
              writeDump_ptr => writeDumpFE
              writeRestart_ptr => writeRestartFE
          case ('dd')
              init_ptr => initDD
              writeDump_ptr => writeDumpDD
              writeRestart_ptr => writeRestartDD
          case default
              write(*,*) 'Simulation type has not yet been defined'
              stop
      end select
      
      end subroutine assignInitAndWrite
************************************************************************
      subroutine Dummy()
      
      end subroutine Dummy
************************************************************************ 
      subroutine Dummy2()
      
      end subroutine Dummy2
************************************************************************
      subroutine Dummy3()
      
      end subroutine Dummy3
************************************************************************
      subroutine writeOutput()
      
C     Subroutine: writeOutput

C     Inputs: None

C     Outputs: None

C     Purpose: Wrapper for writeDump_ptr and writeRestart_ptr
C     Write dump/restart file if current increment is multiple of dumpincrement/restartincrement

      implicit none
      
      if (isMultiple(misc%incrementcurr,misc%dumpincrement)) then
          call writeDump_ptr()
      end if
      if (isMultiple(misc%incrementcurr,misc%restartincrement)) then
          call writeRestart_ptr()
      end if
      
      end subroutine writeOutput
************************************************************************       
      subroutine initCADD()
      
C     Subroutine: initCADD

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize cadd simulation (atoms, fe, dd, detection/passing)

      implicit none

C     local variables      
      character(len=:), allocatable :: pref

      pref = trim(misc%simname)
      call initGeneralChunk(pref)
      call initAtomisticChunk(pref)
      call initFEChunk(pref)
      call initDDChunk(pref)
      call initCADDNoDislChunk()
      call initCADDChunk(pref)
      
      end subroutine initCADD
************************************************************************
      subroutine initCADDNoDisl()
      
C     Subroutine: initCADDNoDisl

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize cadd_nodisl simulation (atoms + fe; no dd)

      implicit none

C     local variables      
      character(len=:), allocatable :: pref

      pref = trim(misc%simname)
      call initGeneralChunk(pref)
      call initAtomisticChunk(pref)
      call initFEChunk(pref)
      call initCADDNoDislChunk()
      
      end subroutine initCADDNoDisl
************************************************************************
      subroutine initFE()
      
C     Subroutine: initFE

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize pure fe simulation (fe nodes and elements; no atoms, discrete dislocations)

      implicit none
      
C     local variables      
      character(len=:), allocatable :: pref

      pref = trim(misc%simname)
      call initGeneralChunk(pref)
      call initFEChunk(pref)

      end subroutine initFE     
************************************************************************
      subroutine initDD()
      
C     Subroutine: initDD

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize dd simulation (fe nodes/elements + dislocations, sources, obstacles; no atoms)

      implicit none

C     local variables      
      character(len=:), allocatable :: pref

      pref = trim(misc%simname)
      call initGeneralChunk(pref)
      call initFEChunk(pref)
      call initDDChunk(pref)
      
      end subroutine initDD    
************************************************************************    
      subroutine initAtomistic()
      
C     Subroutine: initAtomistic

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize atomistic simulation (atoms only, no finite elements/dd)

      implicit none

C     local variables      
      character(len=:), allocatable :: pref

      pref = trim(misc%simname)
      call initGeneralChunk(pref)
      call initAtomisticChunk(pref)
      
      end subroutine initAtomistic
************************************************************************
      subroutine initGeneralChunk(pref)
      
C     Subroutine: initGeneralChunk

C     Inputs: pref --- file prefix (= misc%simname)

C     Outputs: None

C     Purpose: Initialize data needed for all simulations:
C     misc, nodes, materials, groups, and compute

      implicit none

C     input variables      
      character(len=*) :: pref
      
      call initMiscData(pref//suffixmisc)
      call initNodeData(pref//suffixnodes)
      call initMaterialData(pref//suffixmaterials)
      call initGroupData(pref//suffixgroups)
      call initComputeData(pref//suffixcompute)
      
      end subroutine initGeneralChunk
************************************************************************
      subroutine initAtomisticChunk(pref)
      
C     Subroutine: initAtomisticChunk

C     Inputs: pref --- file prefix (= misc%simname)

C     Outputs: None

C     Purpose: Initialize data needed for simulations with atoms:
C     potential style, potentials, interactions, neighbors, damping
      
      implicit none

C     input variables      
      character(len=*) :: pref
      
      call initPotentialStyle()
      call initPotentialData(pref//suffixpotentials)
      call initInteractionData(pref//suffixinteractions)
      call initNeighborData(pref//suffixneighbors)
      call initDampingData(pref//suffixdamping)
      
      end subroutine initAtomisticChunk
************************************************************************
      subroutine initFEChunk(pref)
      
C     Subroutine: initFEChunk

C     Inputs: pref --- file prefix (= misc%simname)

C     Outputs: None

C     Purpose: Initialize data needed for simulations with finite elements:
C     fe library, fe elements, assemble and factor stiffness matrix
      
      implicit none

C     input variables      
      character(len=*) :: pref
            
      call initFELibrary()
      call initFEElementData(pref//suffixfeelements)
      call assignFE(misc%simtype) ! see mod_assign_fe -> need to distinguish between simulations with discrete dislocations
                                  ! and those without
      call initAssembly_ptr()
      call assembleAndFactorAll_ptr()
      
      end subroutine initFEChunk
************************************************************************
      subroutine initDDChunk(pref)
      
C     Subroutine: initDDChunk

C     Inputs: pref --- file prefix (= misc%simname)

C     Outputs: None

C     Purpose: Initialize data needed for simulations with discrete dislocations:
C     slip systems, dislocations, dislocation sources and obstacles, ghost dislocations,
C     escaped dislocations

C     Notes/TODO: Ghost dislocations are not needed for pure DD, but are needed for mod_disl_fields...could be initialized to zero

      implicit none

C     input variables      
      character(len=*) :: pref
      
      call initSlipSysData(pref//suffixslipsys)
      call initDislMiscData(pref//suffixdislmisc)
      call initDislData(pref//suffixdisl)
      call initDislSourceData(pref//suffixdislsource)
      call initDislObsData(pref//suffixdislobs)
      call initGhostDislData(pref//suffixghostdisl)
      call initEscapedDislData(pref//suffixescapeddisl)
      call assignDD(misc%simtype)
      
      end subroutine initDDChunk
************************************************************************
      subroutine initCADDNoDislChunk()
      
C     Subroutine: initCADDNoDislChunk

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize data needed for simulations with pad atoms

      implicit none
      
      call initPad()
      
      end subroutine initCADDNoDislChunk
************************************************************************
      subroutine initCADDChunk(pref)
      
C     Subroutine: initCADDChunk

C     Inputs: None

C     Outputs: None

C     Purpose: Initialize data needed for (cadd) simulations with disl. detection and passing
      
      implicit none
      
C     input variables
      character(len=*) :: pref
      
      call initDetectionData(pref//suffixdetection)
      call initDislIdentData()
      
      end subroutine initCADDChunk
************************************************************************
      subroutine getRestartPrefSuff(pref,suff)
      
C     Subroutine: getRestartPrefSuff

C     Inputs: None

C     Outputs: pref --- prefix for restart file name (= misc%simname)
C              suff --- suffix for restart file name (e.g. at increment 0 -> '.0.restart')

C     Purpose: Generate prefix and suffix for restart file name
            
      character(len=15) :: incrementstr
      character(len=:), allocatable :: suff
      character(len=:), allocatable :: pref
      
      write (incrementstr,'(I0)') misc%incrementcurr
      suff = '.'//trim(incrementstr)//'.restart'
      pref = trim(misc%simname)
      
      end 
************************************************************************
      subroutine writeRestartCADD()
      
C     Subroutine: writeRestartCADD

C     Inputs: None

C     Outputs: None

C     Purpose: Write a set of "restart" files for CADD simulation.
C     These can be used directly as input files
C     (after renaming files) to restart a simulation
C     (slight errors may arise, of course, due to finite precision).

      implicit none

C     local variables
      character(len=:), allocatable :: suff
      character(len=:), allocatable :: pref
      
      call getRestartPrefSuff(pref,suff)
      call writeRestartGeneralChunk(pref,suff)
      call writeRestartAtomisticChunk(pref,suff)
      call writeRestartFEChunk(pref,suff)
      call writeRestartDDChunk(pref,suff)
      call writeRestartCADDChunk(pref,suff)
      
      end subroutine writeRestartCADD
************************************************************************
      subroutine writeRestartCADDNoDisl()
      
C     Subroutine: writeRestartCADDNoDisl

C     Inputs: None

C     Outputs: None

C     Purpose: Write a set of "restart" files for CADD simulation with no dislocations.
C     These can be used directly as input files
C     (after renaming files) to restart a simulation
C     (slight errors may arise, of course, due to finite precision).

      implicit none

C     local variables
      character(len=:), allocatable :: suff
      character(len=:), allocatable :: pref
      
      call getRestartPrefSuff(pref,suff)
      call writeRestartGeneralChunk(pref,suff)
      call writeRestartAtomisticChunk(pref,suff)
      call writeRestartFEChunk(pref,suff)
      
      end subroutine writeRestartCADDNoDisl
************************************************************************
      subroutine writeRestartFE()
      
C     Subroutine: writeRestartFE

C     Inputs: None

C     Outputs: None

C     Purpose: Write a set of "restart" files for pure FE simulation (with no dislocations).
C     These can be used directly as input files
C     (after renaming files) to restart a simulation
C     (slight errors may arise, of course, due to finite precision).

      implicit none
      
C     local variables
      character(len=:), allocatable :: suff
      character(len=:), allocatable :: pref

      call getRestartPrefSuff(pref,suff)
      call writeRestartGeneralChunk(pref,suff)
      call writeRestartFEChunk(pref,suff)
      
      end subroutine writeRestartFE
************************************************************************
      subroutine writeRestartDD()
      
C     Subroutine: writeRestartDD

C     Inputs: None

C     Outputs: None

C     Purpose: Write a set of "restart" files for DD simulation (no atoms).
C     These can be used directly as input files
C     (after renaming files) to restart a simulation
C     (slight errors may arise, of course, due to finite precision).

      implicit none
      
C     local variables
      character(len=:), allocatable :: suff
      character(len=:), allocatable :: pref

      call getRestartPrefSuff(pref,suff)
      call writeRestartGeneralChunk(pref,suff)
      call writeRestartFEChunk(pref,suff)
      call writeRestartDDChunk(pref,suff)
      
      end subroutine writeRestartDD
************************************************************************
      subroutine writeRestartAtomistic()
      
C     Subroutine: writeRestartAtomistic

C     Inputs: None

C     Outputs: None

C     Purpose: Write a set of "restart" files for pure atomistic simulation (no fe/dd).
C     These can be used directly as input files
C     (after renaming files) to restart a simulation
C     (slight errors may arise, of course, due to finite precision).

      implicit none
      
C     local variables
      character(len=:), allocatable :: suff
      character(len=:), allocatable :: pref

      call getRestartPrefSuff(pref,suff)
      call writeRestartGeneralChunk(pref,suff)
      call writeRestartAtomisticChunk(pref,suff)
      
      end subroutine writeRestartAtomistic
************************************************************************
      subroutine writeRestartGeneralChunk(pref,suff)
      
C     Subroutine: writeRestartGeneralChunk

C     Inputs: pref --- prefix for restart file name (= misc%simname)
C             suff --- suffix for restart file name (e.g. at increment 0 -> '.0.restart')

C     Outputs: None

C     Purpose: Write restart files relevant to all simulations:
C     misc, nodes, materials, groups, compute

      implicit none
      
C     input variables      
      character(len=*) :: pref
      character(len=*) :: suff
      
      call writeMiscData(pref//suffixmisc//suff)
      call writeNodeData(pref//suffixnodes//suff)
      call writeMaterialData(pref//suffixmaterials//suff)
      call writeGroupData(pref//suffixgroups//suff)
      call writeComputeData(pref//suffixcompute//suff)
      
      end subroutine writeRestartGeneralChunk
************************************************************************
      subroutine writeRestartAtomisticChunk(pref,suff)
      
C     Subroutine: writeRestartAtomisticChunk

C     Inputs: pref --- prefix for restart file name (= misc%simname)
C             suff --- suffix for restart file name (e.g. at increment 0 -> '.0.restart')

C     Outputs: None

C     Purpose: Write restart files relevant to simulations with atoms:
C     potentials, interactions, neighbors, damping

      implicit none

C     input variables      
      character(len=*) :: pref
      character(len=*) :: suff
      
      call writePotentialData(pref//suffixpotentials//suff)
      call writeInteractionData(pref//suffixinteractions//suff)
      call writeNeighborData(pref//suffixneighbors//suff)
      call writeDampingData(pref//suffixdamping//suff)
      
      end subroutine writeRestartAtomisticChunk
************************************************************************
      subroutine writeRestartFEChunk(pref,suff)
      
C     Subroutine: writeRestartFEChunk

C     Inputs: pref --- prefix for restart file name (= misc%simname)
C             suff --- suffix for restart file name (e.g. at increment 0 -> '.0.restart')

C     Outputs: None

C     Purpose: Write restart file for fe elements

      implicit none

C     input variables      
      character(len=*) :: pref
      character(len=*) :: suff
      
      call writeFEElementData(pref//suffixfeelements//suff)
      
      end subroutine writeRestartFEChunk
************************************************************************
      subroutine writeRestartDDChunk(pref,suff)
      
C     Subroutine: writeRestartDDChunk

C     Inputs: pref --- prefix for restart file name (= misc%simname)
C             suff --- suffix for restart file name (e.g. at increment 0 -> '.0.restart')

C     Outputs: None

C     Purpose: Write restart files relevant to simulations with discrete dislocations:
C     slip systems, dislocations, obstacles, sources, ghost dislocations, escaped dislocations

C     Notes/TODO: Ghost dislocations are not needed for pure DD, but are needed for mod_disl_fields...could be initialized to zero

      implicit none

C     input variables      
      character(len=*) :: pref
      character(len=*) :: suff
      
      call writeSlipSysData(pref//suffixslipsys//suff)
      call writeDislMiscData(pref//suffixdislmisc//suff)
      call writeDislData(pref//suffixdisl//suff)
      call writeDislSourceData(pref//suffixdislsource//suff)
      call writeDislObsData(pref//suffixdislobs//suff)
      call writeGhostDislData(pref//suffixghostdisl//suff)
      call writeEscapedDislData(pref//suffixescapeddisl//suff)
      
      end subroutine writeRestartDDChunk
************************************************************************
      subroutine writeRestartCADDChunk(pref,suff)
      
C     Subroutine: writeRestartCADDChunk

C     Inputs: pref --- prefix for restart file name (= misc%simname)
C             suff --- suffix for restart file name (e.g. at increment 0 -> '.0.restart')

C     Outputs: None

C     Purpose: Write restart files relevant to simulations with dislocation detection and passing

      implicit none

C     input variables      
      character(len=*) :: pref
      character(len=*) :: suff
      
      call writeDetectionData(pref//suffixdetection//suff)
      
      end subroutine writeRestartCADDChunk
************************************************************************
      end module