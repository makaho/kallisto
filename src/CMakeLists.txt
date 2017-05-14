file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

include_directories(../ext/)

add_library(kallisto_core ${sources} ${headers})
target_include_directories(kallisto_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#add_executable(kallisto main.cpp)

find_package( Threads REQUIRED )
target_link_libraries(kallisto_core pthread)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)
    set(HDF5_USE_STATIC_LIBRARIES 1)

    if (UNIX AND NOT APPLE)
        #set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libstdc++")
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

    SET_TARGET_PROPERTIES(kallisto kallisto_core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINK MATCHES static)

if(LINK MATCHES static)
    if (UNIX AND NOT APPLE)
       # target_link_libraries(kallisto librt.a)
    endif()
else()
    if (UNIX AND NOT APPLE)
        #target_link_libraries(kallisto rt)
    endif()
endif(LINK MATCHES static)
