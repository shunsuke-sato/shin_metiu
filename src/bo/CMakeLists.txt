set(SOURCES
bo_mod.f90
)

add_library(${BO_LIB} STATIC ${SOURCES})
add_dependencies(${BO_LIB} ${GLOBAL_LIB})
list(APPEND LIBRARIES ${BO_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)
