add_executable(tests tests.cc)
target_link_libraries(tests PUBLIC util ${LibArchive_LIBRARIES})