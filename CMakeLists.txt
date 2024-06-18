if(photospline_FOUND)

i3_project(LeptonInjector_HNL
    PYTHON_DIR python
    PYTHON_DEST icecube/LeptonInjector_HNL)

set(LIB_${PROJECT_NAME}_SOURCEFILES private/LeptonInjector_HNL/LeptonInjector_HNL.cxx private/LeptonInjector_HNL/serialization.cxx private/LeptonInjector_HNL/HNLDecayEnergy.cxx private/LeptonInjector_HNL/DecayWidths.cxx private/LeptonInjector_HNL/HNLCrossSection.cxx)
set(LIB_${PROJECT_NAME}_PROJECTS serialization icetray dataclasses phys-services sim-services earthmodel-service)
set(LIB_${PROJECT_NAME}_TOOLS boost python photospline)

if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)
    list(APPEND LIB_${PROJECT_NAME}_PROJECTS tableio)
    list(APPEND LIB_${PROJECT_NAME}_SOURCEFILES private/LeptonInjector_HNL/converter/LeptonInjectionConfigurationConverter.cxx)
    add_definitions(-DUSE_TABLEIO)
    message(STATUS "     + tableio support")
endif()

i3_add_library(${PROJECT_NAME}
    ${LIB_${PROJECT_NAME}_SOURCEFILES}
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)

i3_add_pybindings(LeptonInjector_HNL
    private/pybindings/LeptonInjector_HNL.cxx
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS} LeptonInjector_HNL
)

if(BUILD_PYBINDINGS)
  i3_add_pybindings(LeptonInjector_HNL
    private/pybindings/HNLDecayEnergy.cxx
    private/pybindings/LeptonInjector_HNL.cxx
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS} 
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS} LeptonInjector_HNL
  )
endif(BUILD_PYBINDINGS)

i3_test_executable(test
    private/test/main.cxx
    private/test/RangedInjection.cxx
    private/test/VolumeInjection.cxx
    private/test/MultiInjection.cxx
    private/test/Utils.cxx
    private/test/tools.cxx
    USE_TOOLS gsl
    USE_PROJECTS icetray dataclasses dataio phys-services earthmodel-service LeptonInjector_HNL
)

else()
    colormsg(RED "+-- Omitting LeptonInjector, photospline *not* found")
endif()
