set(SHADER alWriteIdsDriver)

set(SRC ${SHADER}.cpp)

add_library(${SHADER} SHARED ${SRC})

target_link_libraries(${SHADER} ai)
set_target_properties(${SHADER} PROPERTIES PREFIX "")

install(TARGETS ${SHADER} DESTINATION ${CMAKE_INSTALL_EXEC_PREFIX}/plugins)
