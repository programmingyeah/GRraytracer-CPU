cmake_minimum_required(VERSION 3.10)
project(RelativisticRaymarcher)

set(CMAKE_CXX_STANDARD 11)

# Include directories
include_directories(include)

# Add GLM to the project (assuming glm folder is in the project directory)
set(GLM_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include/glm")
include_directories(${GLM_INCLUDE_DIRS})

# Define source files
set(SOURCES
    src/main.cpp
    src/raymarcher.cpp
)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# No need to link GLM as it's header-only

# FFmpeg command to generate video from frames after the build
#add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND ffmpeg -r 30 -f image2 -s 800x600 -i output/frame%d.ppm -vcodec libx264 -crf 25 -pix_fmt yuv420p output/output.mp4
#    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#    COMMENT "Generating video from frames using FFmpeg..."
#)
