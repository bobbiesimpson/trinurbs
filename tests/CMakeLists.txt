add_executable(trinurbstest trinurbstest.cpp)
add_executable(multiscaletest multiscaleforesttest.cpp)
add_executable(periodic_connectivity_test periodictest.cpp)

target_link_libraries(trinurbstest ${THIS_LINK_LIBS} trinurbs)
target_link_libraries(multiscaletest ${THIS_LINK_LIBS} trinurbs)
target_link_libraries(periodic_connectivity_test ${THIS_LINK_LIBS} trinurbs)

add_test(trinurbstest trinurbstest)
add_test(multiscaletest multiscaletest)
add_test(periodic_connectivity_test periodic_connectivity_test)
