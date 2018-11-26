


set(_orbweaver_sources
	allvars.cxx
	calcorbitprops.cxx
	cleanorbits.cxx
	cosmofuncs.cxx
	interpfuncs.cxx
	io.cxx
	main.cxx
	ui.cxx)


add_executable(orbweaver ${_orbweaver_sources})
target_link_libraries(orbweaver nbodylib ${ORB_LIBS})
target_compile_definitions(orbweaver PRIVATE ${ORB_DEFINES})
if (ORB_CXX_FLAGS)
	set_target_properties(orbweaver PROPERTIES COMPILE_FLAGS ${ORB_CXX_FLAGS})
endif()
if (ORB_LINK_FLAGS)
	set_target_properties(orbweaver PROPERTIES LINK_FLAGS ${ORB_LINK_FLAGS})
endif()
set_target_properties(orbweaver PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})