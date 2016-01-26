set(SHADER alWriteIds)

set(SRC ${SHADER}.cpp)
set(MTD ${SHADER}.mtd)

add_library(${SHADER} SHARED ${SRC})

target_link_libraries(${SHADER} ai)
set_target_properties(${SHADER} PROPERTIES PREFIX "")

install(TARGETS ${SHADER} DESTINATION ${DSO_INSTALL_DIR})
install(FILES ${MTD} DESTINATION ${MTD_INSTALL_DIR})

