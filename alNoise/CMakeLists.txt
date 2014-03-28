set(SHADERS alFractal alCellNoise alFlowNoise alGaborNoise)

foreach(SHADER ${SHADERS})
	set(SRC ${SHADER}.cpp ../common/alUtil.cpp)
	set(UI ${CMAKE_CURRENT_SOURCE_DIR}/${SHADER}.ui)
	set(MTD ${CMAKE_CURRENT_BINARY_DIR}/${SHADER}.mtd)
	set(AE ${CMAKE_CURRENT_BINARY_DIR}/${SHADER}Template.py)
	set(SPDL ${CMAKE_CURRENT_BINARY_DIR}/${SHADER}.spdl)
	set(KARGS ${CMAKE_CURRENT_BINARY_DIR}/${SHADER}.args)

	add_library(${SHADER} SHARED ${SRC})

	target_link_libraries(${SHADER} ai)
	set_target_properties(${SHADER} PROPERTIES PREFIX "")

	add_custom_command(OUTPUT ${MTD} COMMAND python ARGS ${CMAKE_SOURCE_DIR}/uigen.py ${UI} ${MTD} ${AE} ${SPDL} ${KARGS} DEPENDS ${UI})
	add_custom_target(${SHADER}UI ALL DEPENDS ${MTD})

	install(TARGETS ${SHADER} DESTINATION ${MTOA_SHADERS})
	install(FILES ${MTD} DESTINATION ${MTOA_SHADERS})
	install(FILES ${AE} DESTINATION ${MTOA_UI})
endforeach()