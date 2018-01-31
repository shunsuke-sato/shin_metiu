add_subdirectory("./global")
add_subdirectory("./tdse")
add_subdirectory("./bo")
add_subdirectory("./io")

# Executable binary
add_executable(${TARGET_NAME} main.f90)
target_link_libraries(${TARGET_NAME} ${LIBRARIES})
add_dependencies(${TARGET_NAME} ${LIBRARIES})
