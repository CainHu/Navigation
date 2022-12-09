include_directories("E:\\eigen3")


# aux_source_directory(. SOURCES)
# add_library(GESKF ${SOURCES})

add_library(GESKF state.cpp covariance.cpp fuse_position.cpp fuse_velocity.cpp fuse_magnet.cpp fuse_declination.cpp)
target_include_directories(GESKF PUBLIC "${PROJECT_SOURCE_DIR}\\global_eskf\\include")

target_link_libraries(GESKF ESKF)