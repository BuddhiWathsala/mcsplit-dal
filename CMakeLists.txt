cmake_minimum_required(VERSION 3.26)
project(main)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_FLAGS "-largp")

set(BOOST_INCLUDE_DIR /opt/homebrew/Cellar/boost/1.85.0/include)
set(BOOST_LIBRARY_DIR /opt/homebrew/Cellar/boost/1.85.0/lib)
set(BOOST_INCLUDEDIR /opt/homebrew/Cellar/boost/1.85.0/include/boost)
set(ARGP_INCLUDE_DIR /opt/homebrew/Cellar/argp-standalone/1.3/include)
set(ARGP_LIBRARY_DIR /opt/homebrew/Cellar/argp-standalone/1.3/lib)
set(ARGP_INCLUDEDIR /opt/homebrew/Cellar/argp-standalone/1.3/lib)

add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

# Add -O3 optimization flag for Release builds
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

include_directories(${BOOST_INCLUDE_DIR})
link_directories(${BOOST_LIBRARY_DIR})
find_package(Boost COMPONENTS system filesystem REQUIRED)

include_directories(${ARGP_INCLUDE_DIR})
link_directories(${ARGP_LIBRARY_DIR})
#find_package()


#add_executable(Boost main.cpp graph.cpp graph.h modular.cpp modular.h primes.h search_ops.cpp search_ops.h structs.h utils.cpp utils.h)
#add_executable(Argp main.cpp graph.cpp graph.h modular.cpp modular.h primes.h search_ops.cpp search_ops.h structs.h utils.cpp utils.h)
add_executable(main mcsplit_plus.cpp graph.cpp graph.h)