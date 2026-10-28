find_library(ARIADNE_CONSTRAINT_SOLVING_LIBRARY NAMES ariadne-constraint-solving)

find_package(Threads REQUIRED)

find_path(ARIADNE_CONSTRAINT_SOLVING_INCLUDE_DIR ariadne-constraint-solving.hpp PATH_SUFFIXES ariadne-constraint-solving)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneConstraintSolving DEFAULT_MSG
  ARIADNE_CONSTRAINT_SOLVING_LIBRARY
  ARIADNE_CONSTRAINT_SOLVING_INCLUDE_DIR
)

if(AriadneConstraintSolving_FOUND)
  get_filename_component(ARIADNE_CONSTRAINT_SOLVING_INCLUDE_PARENT_DIR ${ARIADNE_CONSTRAINT_SOLVING_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_CONSTRAINT_SOLVING_INCLUDE_DIRS
    ${ARIADNE_CONSTRAINT_SOLVING_INCLUDE_PARENT_DIR}
    ${ARIADNE_CONSTRAINT_SOLVING_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneConstraintSolving::ariadne-constraint-solving)
    add_library(AriadneConstraintSolving::ariadne-constraint-solving UNKNOWN IMPORTED)
    set_target_properties(AriadneConstraintSolving::ariadne-constraint-solving PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_CONSTRAINT_SOLVING_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_CONSTRAINT_SOLVING_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "Threads::Threads"
    )
  endif()
endif()
