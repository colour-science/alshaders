add_library(cryptomatte_filter SHARED cryptomatte_filter.cpp MurmurHash3.cpp)
target_link_libraries(cryptomatte_filter ai)
set_target_properties(cryptomatte_filter PROPERTIES PREFIX "")
install(TARGETS cryptomatte_filter DESTINATION ${DSO_INSTALL_DIR})
install(DIRECTORY nuke DESTINATION ${INSTALL_DIR})

