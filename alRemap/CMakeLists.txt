set(SHADER_F alRemapFloat)
set(SHADER_C alRemapColor)

set(SRC_F ${SHADER_F}.cpp ../common/alUtil.cpp)
set(SRC_C ${SHADER_C}.cpp ../common/alUtil.cpp)
set(MTD_F ${SHADER_F}.mtd)
set(MTD_C ${SHADER_C}.mtd)
set(MEL_F AE${SHADER_F}Template.mel)
set(MEL_C AE${SHADER_C}Template.mel)

add_library(${SHADER_F} SHARED ${SRC_F})
add_library(${SHADER_C} SHARED ${SRC_C})

target_link_libraries(${SHADER_F} ai)
target_link_libraries(${SHADER_C} ai)
set_target_properties(${SHADER_F} PROPERTIES PREFIX "")
set_target_properties(${SHADER_C} PROPERTIES PREFIX "")

install(TARGETS ${SHADER_F} ${SHADER_C} DESTINATION ${MTOA_SHADERS})
install(FILES ${MTD_F} ${MTD_C} DESTINATION ${MTOA_SHADERS})
install(FILES ${MEL_F} ${MTD_C} DESTINATION ${MTOA_SCRIPTS})
