set(DRIVER alWriteIdsDriver)

set(SRC ${DRIVER}.cpp)

add_library(${DRIVER} SHARED ${SRC})

target_link_libraries(${DRIVER} ai)
set_target_properties(${DRIVER} PROPERTIES PREFIX "")

install(TARGETS ${DRIVER} DESTINATION ${DSO_INSTALL_DIR})
