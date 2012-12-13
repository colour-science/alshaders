set(SHADERS alIES)
include_directories(/Users/anders/vfx/boost/1.43.0/include)
foreach(SHADER ${SHADERS})
	set(SRC ${SHADER}.cpp ../common/alUtil.cpp ../common/IES.cpp)
	set(MTD ${SHADER}.mtd)
	set(MEL AE${SHADER}Template.mel)

	add_library(${SHADER} SHARED ${SRC})

	target_link_libraries(${SHADER} ai)
	set_target_properties(${SHADER} PROPERTIES PREFIX "")

	install(TARGETS ${SHADER} DESTINATION ${MTOA_SHADERS})
	install(FILES ${MTD} DESTINATION ${MTOA_SHADERS})
	install(FILES ${MEL} DESTINATION ${MTOA_SCRIPTS})
endforeach()

