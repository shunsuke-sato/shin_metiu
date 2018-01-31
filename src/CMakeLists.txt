add_subdirectory("./global")
add_subdirectory("./tdse")
add_subdirectory("./bo")
add_subdirectory("./io")

list(APPEND LIBRARIES 
  ${GLOBAL_LIB}
  ${IO_LIB}
  ${TDSE_LIB}
  ${BO_LIB}
)

set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)

# Executable binary
add_executable(${TARGET_NAME} main.f90)


target_link_libraries(${TARGET_NAME} ${LIBRARIES})
add_dependencies(${TARGET_NAME} ${LIBRARIES})
