add_executable(trinurbstest trinurbstest.cpp)
add_executable(multiscaletest multiscaleforesttest.cpp)
add_executable(periodic_connectivity_test periodictest.cpp)
add_executable(mpi_fe_test mpiFEassemblytest.cpp)

target_link_libraries(trinurbstest ${THIS_LINK_LIBS} trinurbs)
target_link_libraries(multiscaletest ${THIS_LINK_LIBS} trinurbs)
target_link_libraries(periodic_connectivity_test ${THIS_LINK_LIBS} trinurbs)
target_link_libraries(mpi_fe_test ${THIS_LINK_LIBS} ${MPI_CXX_LIBRARIES} ${Trilinos_LIBRARIES} trinurbs)

#add_test(trinurbstest trinurbstest)
#add_test(multiscaletest multiscaletest)
#add_test(periodic_connectivity_test periodic_connectivity_test)
