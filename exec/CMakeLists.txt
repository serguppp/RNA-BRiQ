add_executable(BRiQ_InitPDB GenerateInitModel.cpp)
add_executable(BRiQ_AssignSS TestAssignment.cpp)
add_executable(BRiQ_Predict RunPredict.cpp)
add_executable(BRiQ_Refinement RunRefinement.cpp)
add_executable(BRiQ_Energy RunEnergy.cpp)

add_dependencies(BRiQ_InitPDB combinelib)
add_dependencies(BRiQ_AssignSS combinelib)
add_dependencies(BRiQ_Predict combinelib)
add_dependencies(BRiQ_Refinement combinelib)
add_dependencies(BRiQ_Energy combinelib)

target_link_libraries (BRiQ_InitPDB briq)
target_link_libraries (BRiQ_AssignSS briq)
target_link_libraries (BRiQ_Predict briq)
target_link_libraries (BRiQ_Refinement briq)
target_link_libraries (BRiQ_Energy briq)

# Embed source identity in run.json; callers may provide an explicit build value.
if(NOT BRIQ_SOURCE_VERSION)
  execute_process(COMMAND git describe --always --dirty
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    OUTPUT_VARIABLE BRIQ_SOURCE_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if(NOT BRIQ_SOURCE_VERSION)
  set(BRIQ_SOURCE_VERSION "unknown")
endif()
set_property(TARGET BRiQ_Energy APPEND PROPERTY COMPILE_DEFINITIONS
  BRIQ_SOURCE_VERSION="${BRIQ_SOURCE_VERSION}")
