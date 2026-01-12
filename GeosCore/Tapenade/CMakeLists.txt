# Tapenade/CMakeLists.txt

# Add a library containing the TAPENADE runtime source files.
# Make sure to list all C and Fortran files from TAPENADE here.
add_library(tapenade_runtime
    adStack.c
    # List any other C or Fortran source files from TAPENADE here, e.g.:
    # adFortranSupport.f90
)

# This makes the 'Tapenade' directory and its headers available
# to the compiler when it processes the source files in this library.
target_include_directories(tapenade_runtime PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# This is an important line. The GEOS-Chem model often uses specific
# compile flags. You can inherit them from the parent project.
target_compile_options(tapenade_runtime PUBLIC
    -fPIE
)
