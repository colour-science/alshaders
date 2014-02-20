set(SHADERS alFractal alCellNoise alFlowNoise alGaborNoise)

foreach(SHADER ${SHADERS})
	set(SRC ${SHADER}.cpp ../common/alUtil.cpp)
	set(UI ${CMAKE_CURRENT_SOURCE_DIR}/${SHADER}.ui.py)
	set(MTD ${CMAKE_CURRENT_BINARY_DIR}/${SHADER}.mtd)
	set(AE ${CMAKE_CURRENT_BINARY_DIR}/${SHADER}Template.py)
	set(SPDL ${CMAKE_CURRENT_BINARY_DIR}/${SHADER}.spdl)
	
	add_library(${SHADER} SHARED ${SRC})

	target_link_libraries(${SHADER} ai)
	set_target_properties(${SHADER} PROPERTIES PREFIX "")

	add_custom_command(OUTPUT ${MTD} COMMAND python ARGS ${CMAKE_SOURCE_DIR}/uigen.py ${UI} ${MTD} ${AE} ${SPDL} DEPENDS ${UI})
	add_custom_target(${SHADER}UI ALL DEPENDS ${MTD})

	if(DEFINED MTOA_SHADERS)
	install(TARGETS ${SHADER} DESTINATION ${MTOA_SHADERS})
	install(FILES ${MTD} DESTINATION ${MTOA_SHADERS})
	install(FILES ${AE} DESTINATION ${MTOA_UI})
	endif()
	if(DEFINED HTOA_SHADERS)
	install(TARGETS ${SHADER} DESTINATION ${HTOA_SHADERS})
	install(FILES ${MTD} DESTINATION ${HTOA_SHADERS})
	endif()
endforeach()
