set(SHADERS alRemapFloat alRemapColor)

foreach(SHADER ${SHADERS})
    set(SRC ${SHADER}.cpp ../common/alUtil.cpp)
    set(MTD ${SHADER}.mtd)
    set(UI ${SHADER}Template.py)
    
    add_library(${SHADER} SHARED ${SRC})
    
    target_link_libraries(${SHADER} ai)
    set_target_properties(${SHADER} PROPERTIES PREFIX "")
    
    install(TARGETS ${SHADER} DESTINATION ${MTOA_SHADERS})
    install(FILES ${MTD} DESTINATION ${MTOA_SHADERS})
    install(FILES ${UI} DESTINATION ${MTOA_UI})
endforeach()

