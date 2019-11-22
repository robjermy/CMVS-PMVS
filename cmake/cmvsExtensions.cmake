# Adapted from https://github.com/kracejic/cleanCppProject

# Helper macro for creating convenient targets
find_program(GDB_PATH gdb)

# Adds -run and -dbg targets
macro(addRunAndDebugTargets TARGET)
  add_custom_target(${TARGET}-run
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    USES_TERMINAL
    DEPENDS ${TARGET}
    COMMAND ./${TARGET})

  # convenience run gdb target
  if(GDB_PATH)
    add_custom_target(${TARGET}-gdb
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      USES_TERMINAL
      DEPENDS ${TARGET}
      COMMAND ${GDB_PATH} ./${TARGET})
  endif()
endmacro()

macro(ExternalGitProject LIBNAME REPOSITORY GIT_TAG)
  set(${LIBNAME}_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${LIBNAME})

  # clone repository if not done
  if(IS_DIRECTORY ${${LIBNAME}_SOURCE_DIR})
    message(STATUS "Already downloaded: ${REPOSITORY}")
  else()
    message(STATUS "Clonning: ${REPOSITORY}")
    execute_process(
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND ${GIT_EXECUTABLE} clone -b ${GIT_TAG} --single-branch --depth 1 ${REPOSITORY} ${CMAKE_CURRENT_SOURCE_DIR}/${LIBNAME}
      )
  endif()
endmacro()