set(SOURCES
model_parameters.f90
math_parameters.f90
control_parameters.f90
io_global.f90
)

add_library(${GLOBAL_LIB} STATIC ${SOURCES})
add_dependencies(${GLOBAL_LIB} ${IO_LIB})
list(APPEND LIBRARIES ${GLOBAL_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)
