# create the executables for the apps

cmessage( STATUS "Creating daq_tps_clustering app..." )
add_executable( daq_tps_clustering ${CMAKE_CURRENT_SOURCE_DIR}/daq_tps_clustering.cpp )
target_link_libraries( daq_tps_clustering utils-libraries )
install( TARGETS daq_tps_clustering DESTINATION bin )
