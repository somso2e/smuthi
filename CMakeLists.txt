cmake_minimum_required(VERSION 3.15)
project(smuthi LANGUAGES C CXX Fortran)

# Define the directory containing the Fortran code
set(FORTRAN_SRC_DIR ${CMAKE_SOURCE_DIR}/smuthi/linearsystem/tmatrix/nfmds/NFM-DS/TMATSOURCES)

# Run the existing Makefile inside TMATSOURCES
add_custom_target(nfmds ALL
    COMMAND make -C ${FORTRAN_SRC_DIR}
    WORKING_DIRECTORY ${FORTRAN_SRC_DIR}
)

# Install the compiled library
install(DIRECTORY ${FORTRAN_SRC_DIR}
    DESTINATION smuthi/linearsystem/tmatrix/nfmds
)

