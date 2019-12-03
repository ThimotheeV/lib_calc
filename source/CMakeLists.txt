project(libsource)

aux_source_directory(. SOURCES_LIST)
if (WIN32)
    add_library(${PROJECT_NAME} STATIC ${SOURCES_LIST})
else()
    add_library(${PROJECT_NAME} SHARED ${SOURCES_LIST})
endif()

SET(MY_DEBUG_OPTIONS "-Wall" "-Werror" "-Wfatal-errors" "-g")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:${MY_DEBUG_OPTIONS}>")