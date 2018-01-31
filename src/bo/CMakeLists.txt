set(SOURCES
bo_mod.f90
)

add_library(${BO_LIB} STATIC ${SOURCES})
add_dependencies(${BO_LIB} ${GLOBAL_LIB})

target_link_libraries(${BO_LIB} ${GLOBAL_LIB})
