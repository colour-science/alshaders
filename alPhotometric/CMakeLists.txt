set(SHADERS alPhotometric)

foreach(SHADER ${SHADERS})
    set(SRC ${SHADER}.cpp ../common/alUtil.cpp Photometric.cpp)
    set(MTD ${SHADER}.mtd)
    set(UI ${SHADER}Template.py)
    
    add_library(${SHADER} SHARED ${SRC})
    
    target_link_libraries(${SHADER} ai)
    set_target_properties(${SHADER} PROPERTIES PREFIX "")
    
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

