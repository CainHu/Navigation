# include_directories("E:\\eigen3")
# include_directories("${PROJECT_SOURCE_DIR}\\eskf\\include")

# aux_source_directory(. SOURCES)
# add_executable(${PROJECT_NAME} ${SOURCES})

# set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}\\bin)


include_directories("E:\\eigen3")
# include_directories("${PROJECT_SOURCE_DIR}\\eskf\\include")


aux_source_directory(. SOURCES)
add_library(ESKF ${SOURCES})
target_include_directories(ESKF PUBLIC "${PROJECT_SOURCE_DIR}\\eskf\\include")