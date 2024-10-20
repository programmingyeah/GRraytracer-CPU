#ifndef RAYMARCHER_HPP
#define RAYMARCHER_HPP

#include <iostream>

#include <glm/glm.hpp>
#include <vector>

struct Ray {
    glm::vec3 origin;
    glm::vec3 direction;
};

glm::vec3 raymarch(const Ray& ray, const glm::vec3& cameraVelocity, double speedOfLight);

void renderFrame(int frameNumber, const glm::vec3& cameraVelocity, double speedOfLight);

#endif
