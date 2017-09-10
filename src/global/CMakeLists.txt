set(SOURCES
model_parameters.f90
math_parameters.f90
)

add_library(${GLOBAL_LIB} STATIC ${SOURCES})
list(APPEND LIBRARIES ${GLOBAL_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)
