set(SOURCES
math_parameters.f90
model_parameters.f90
global_variables.f90
io_global.f90
)

add_library(${GLOBAL_LIB} STATIC ${SOURCES})
add_dependencies(${GLOBAL_LIB} ${IO_LIB})

