#include "raymarcher.hpp"
#include <iostream>
#include <glm/glm.hpp>

int main() {
    const int numFrames = 1;
    double c = 299792458.0;

    glm::vec3 cameraVelocity(0.0f, 0.0f, 0.0f);

    for (int frame = 0; frame < numFrames; ++frame) {
        std::cout << "Rendering frame: " << frame << std::endl;
        renderFrame(frame, cameraVelocity, c);
    }

    std::cout << "Rendering complete!" << std::endl;

    return 0;
}
