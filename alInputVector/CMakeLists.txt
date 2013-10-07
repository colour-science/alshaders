set(SHADERS alInputVector)

foreach(SHADER ${SHADERS})
	set(SRC ${SHADER}.cpp ../common/alUtil.cpp)
	set(UI ${CMAKE_CURRENT_SOURCE_DIR}/${SHADER}.ui.py)
	set(MTD ${CMAKE_CURRENT_SOURCE_DIR}/${SHADER}.mtd)
	set(AE ${CMAKE_CURRENT_SOURCE_DIR}/${SHADER}Template.py)
	set(SPDL ${CMAKE_CURRENT_SOURCE_DIR}/${SHADER}.spdl)

	add_library(${SHADER} SHARED ${SRC})

	target_link_libraries(${SHADER} ai)
	set_target_properties(${SHADER} PROPERTIES PREFIX "")

	add_custom_command(OUTPUT ${MTD} COMMAND ${CMAKE_SOURCE_DIR}/uigen.py ${UI} ${MTD} ${AE} ${SPDL} DEPENDS ${UI})
	add_custom_target(${SHADER}UI ALL DEPENDS ${MTD})

	install(TARGETS ${SHADER} DESTINATION ${MTOA_SHADERS})
	install(FILES ${MTD} DESTINATION ${MTOA_SHADERS})
	install(FILES ${UI} DESTINATION ${MTOA_UI})
endforeach()

