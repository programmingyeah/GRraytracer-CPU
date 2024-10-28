#ifndef RAYMARCHER_HPP
#define RAYMARCHER_HPP

#include <iostream>

#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <vector>

struct Ray {
    glm::vec3 origin;
    glm::vec3 direction;
};

bool loadHDRSkybox(const char* hdrPath);

glm::vec3 raymarch(const Ray& ray, const glm::vec3& cameraVelocity, double speedOfLight);

void renderFrame(int frameNumber, const glm::vec3& cameraVelocity, double speedOfLight);

#endif
