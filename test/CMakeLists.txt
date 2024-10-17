if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_MACOSX_RPATH 1)
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    
    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endif()

add_executable(tests_streambuffer tests_streambuffer.cc)
add_executable(tests_feature_extraction tests_feature_extraction.cc)
add_executable(tests_streamcompressor tests_streamcompressor.cc)
add_executable(tests_gbdlib tests_gbdlib.cc)

target_link_libraries(tests_streambuffer PRIVATE util ${LibArchive_LIBRARIES})
target_link_libraries(tests_feature_extraction PRIVATE util solver extract ${LibArchive_LIBRARIES})
target_link_libraries(tests_streamcompressor PRIVATE util ${LibArchive_LIBRARIES})
target_link_libraries(tests_gbdlib PRIVATE util ${LIBS})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
