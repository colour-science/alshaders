set(UI alShaders.py)

set(MEL alsUtil.mel shelf_alShaders.mel alsIdWindow.mel)
set(ICONS als_lg_0.png als_lg_1.png als_lg_2.png als_lg_3.png als_lg_4.png als_lg_5.png als_lg_6.png als_lg_7.png als_lg_8.png als_id.png)

install(FILES ${UI} DESTINATION ${MTOA_UI})
install(FILES ${MEL} DESTINATION ${MTOA_SCRIPTS})
install(FILES ${ICONS} DESTINATION ${MTOA_ICONS})