set(SOURCES
tdse_mod.f90
variables_for_tdse.f90
)

add_library(${TDSE_LIB} STATIC ${SOURCES})
add_dependencies(${TDSE_LIB} ${GLOBAL_LIB})

