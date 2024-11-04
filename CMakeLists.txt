cmake_minimum_required(VERSION 3.10)  # Adjust based on your version of CMake
project(gevolution)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find packages
find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(FFTW3 REQUIRED)  # This assumes that FFTW3 has FindFFTW3.cmake or similar
find_package(GSL REQUIRED)    # This assumes GSL has a FindGSL.cmake

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/LATfield2)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${FFTW3_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})

# Source files
set(SOURCE_FILES main.cpp)

# Executable definitions
add_executable(gevolution ${SOURCE_FILES})

# Compiler options
add_compile_options(-O3)

# Definitions for LATfield2 and gevolution
add_compile_definitions(FFT3D HDF5 PHINONLINEAR BENCHMARK EXACT_OUTPUT_REDSHIFTS)
# Uncomment the following lines as needed
# add_compile_definitions(H5_HAVE_PARALLEL)
# add_compile_definitions(EXTERNAL_IO)
# add_compile_definitions(SINGLE)
# add_compile_definitions(VELOCITY)
# add_compile_definitions(COLORTERMINAL)
# add_compile_definitions(CHECK_B)
# add_compile_definitions(HAVE_CLASS)
# add_compile_definitions(HAVE_HEALPIX)

# Link libraries
target_link_libraries(gevolution ${HDF5_LIBRARIES} ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

# Additional executables if needed
add_executable(lccat lccat.cpp)
target_link_libraries(lccat ${GSL_LIBRARIES})

add_executable(lcmap lcmap.cpp)
target_link_libraries(lcmap ${FFTW3_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES})
target_compile_options(lcmap PRIVATE -fopenmp)
