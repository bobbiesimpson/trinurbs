add_executable(trinurbstest trinurbstest.cpp)

target_link_libraries(trinurbstest ${THIS_LINK_LIBS} trinurbs)

add_test(trinurbstest trinurbstest)
