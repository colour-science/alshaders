set(SHADER alSpectral)        
set(SRC ${SHADER}.cpp  SPD.cpp ../common/alUtil.cpp)
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

install(TARGETS ${SHADER} DESTINATION ${DSO_INSTALL_DIR})
install(FILES ${MTD} DESTINATION ${MTD_INSTALL_DIR})
install(FILES ${AE} DESTINATION ${AE_INSTALL_DIR})
install(FILES ${SPDL} DESTINATION ${SPDL_INSTALL_DIR})
install(FILES ${KARGS} DESTINATION ${ARGS_INSTALL_DIR})


set(SHADER alSpectralColor)        
set(SRC ${SHADER}.cpp  SPD.cpp ../common/alUtil.cpp)
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

install(TARGETS ${SHADER} DESTINATION ${DSO_INSTALL_DIR})
install(FILES ${MTD} DESTINATION ${MTD_INSTALL_DIR})
install(FILES ${AE} DESTINATION ${AE_INSTALL_DIR})
install(FILES ${SPDL} DESTINATION ${SPDL_INSTALL_DIR})
install(FILES ${KARGS} DESTINATION ${ARGS_INSTALL_DIR})
