add_executable(tests_streambuffer tests_streambuffer.cc)
add_executable(tests_cnfbasefeatures tests_cnfbasefeatures.cc)
target_link_libraries(tests_streambuffer PUBLIC util ${LibArchive_LIBRARIES})
target_link_libraries(tests_cnfbasefeatures PUBLIC util ${LibArchive_LIBRARIES})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
