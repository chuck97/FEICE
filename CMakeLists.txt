cmake_minimum_required(VERSION 3.17.0)

project(DYNAMICS_TESTS)

enable_language (CXX C)
set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(Service)
add_subdirectory(Mesh)
add_subdirectory(Configuration)
add_subdirectory(DynamicsTest)
add_subdirectory(Advection)
add_subdirectory(Momentum)
add_subdirectory(Forcing)
add_subdirectory(Headers)

include_directories("Headers")

set(INMOST_BUILD_ARCH "INMOST-intel-build")

set(INMOST_SRC "/data90t/geosci/spetrov/INMOST_BUILD_PVTU/INMOST")
set(INMOST_BUILD "/data90t/geosci/spetrov/INMOST_BUILD_PVTU/${INMOST_BUILD_ARCH}")
set(PETSC_DIR "/data90t/geosci/spetrov/INMOST_BUILD_PVTU/petsc-3.5.2")
set(PETSC_ARCH "linux-intel-debug")
set(JSON_DIR "/data90t/geosci/spetrov/json_lib/nlohmann_json_cmake_fetchcontent/include")
set(PROJ4 "/data90t/geosci/spetrov/proj-4/proj-4-build")

include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})
include_directories(${JSON_DIR}/)


include_directories(${PETSC_DIR}/include)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)

include_directories(${PROJ4}/include)

set(
	SOURCE
	${SOURCE}
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
	${HEADER}
)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE})

target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2clapack.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2cblas.a)
target_link_libraries(${EXEC_NAME} ${PROJ4}/lib/libproj.a)
target_link_libraries(${EXEC_NAME} -lnetcdf)
target_link_libraries(${EXEC_NAME} -lhdf5)
target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
