# create the executables for the apps

cmessage( STATUS "Creating daq_tps_clustering app..." )
add_executable( daq_tps_clustering ${CMAKE_CURRENT_SOURCE_DIR}/daq_tps_clustering.cpp )
target_link_libraries( daq_tps_clustering utils-libraries )
install( TARGETS daq_tps_clustering DESTINATION bin )


cmessage( STATUS "Creating extract_neutrino_candidates app..." )
add_executable( extract_neutrino_candidates ${CMAKE_CURRENT_SOURCE_DIR}/extract_neutrino_candidates.cpp )
target_link_libraries( extract_neutrino_candidates utils-libraries )
install( TARGETS extract_neutrino_candidates DESTINATION bin )
