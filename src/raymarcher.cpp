#include "raymarcher.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>


const double MAX_DIST = 100.0f;
const int MAX_STEPS = 100;
const double EPSILON = 0.001f;

struct State {
    glm::dvec4 pos;
    glm::dvec4 vel;
    int index;
};

class Data {   
public:
    float ***gamma;   
    float **g;        
    State *s;   

    Data() {
        s = new State();

        gamma = new float**[4];
        for (int i = 0; i < 4; ++i) {
            gamma[i] = new float*[4];
            for (int j = 0; j < 4; ++j) {
                gamma[i][j] = new float[4];
            }
        }

        g = new float*[4];
        for (int i = 0; i < 4; ++i) {
            g[i] = new float[4];
        }
    }

    ~Data() {
        deallocateMemory();
    }

    void deallocateMemory() {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                delete[] gamma[i][j];
            }
            delete[] gamma[i];
        }
        delete[] gamma;
        for (int i = 0; i < 4; ++i) {
            delete[] g[i];
        }
        delete[] g;

        delete s;
    }
};

double W_0(double x) //laaaambert
{
    if (x == 0) return 0;

    double e = exp(1);
    unsigned int N = 3;

    double W = 1;

    if (x > e) {
        double l = log(x);
        W = l - log(l);
    } else if (x > 0) {
        W = x/e;
    } else {
        W = e*x*log(1+sqrt(1+e*x))/(2+e*x+sqrt(1+e*x));
    }

    for (unsigned int n = 1; n <= N; n++) 
    {
        W = W*(1+log(x/W))/(1+W);
    }

    return W;
}

double DW_0(double x) { //derivative of lambert
    double W = W_0(x);
    return W/(x*(1+W));
}

double sdfSphere(glm::dvec3*& p, glm::dvec3 center, double radius) {
    return glm::length((*p) - center) - radius;
}

double sceneSDF(glm::dvec3*& p) {
    double dist1 = sdfSphere(p, glm::dvec3(0, 0, -3), 2.0f);
    double dist2 = sdfSphere(p, glm::dvec3(2, 0, -3), 2.0f);
    return dist1;
}

glm::dmat4 Metric(State*& s) {

    glm::dmat4 metric = glm::dmat4(1.0);
    metric[0][0] = -1.0;

    double r = (*s).pos.y;
    double theta = (*s).pos.z;
    double phi = (*s).pos.w;

   /* metric[0][0] = -(1-rs/r);
    metric[1][1] = 1/(1-rs/r);
    metric[2][2] = r * r;
    metric[3][3] = metric[2][2] * sin(theta) * sin(theta);*/

    return metric;
}

double partialg(State*& state, int mu, int nu, int direction)
{
    glm::dvec4 pos = (*state).pos;

    double delta = 0.01;
    
    State forward = *state;
    forward.pos = pos + glm::dvec4(((direction == 0) ? delta : 0.0),
    ((direction == 1) ? delta : 0.0),
    ((direction == 2) ? delta : 0.0),
    ((direction == 3) ? delta : 0.0));
    
    State backward = *state;
    backward.pos = pos - glm::dvec4(((direction == 0) ? delta : 0.0),
    ((direction == 1) ? delta : 0.0),
    ((direction == 2) ? delta : 0.0),
    ((direction == 3) ? delta : 0.0));

    State* pF = &forward;
    State* pB = &backward;

    glm::dmat4 gFor = Metric(pF);
    glm::dmat4 gBac = Metric(pB);

    return (gFor[mu][nu] - gBac[mu][nu]) / (2 * delta);
}

void Gamma(State*& state, double*** gamma)
{
    glm::dmat4 g = Metric(state);
    glm::dmat4 ginv = inverse(g);

    for (int lambda = 0; lambda < 4; lambda++)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            for (int nu = 0; nu < 4; nu++)
            {
                double gamma_value = 0.0;

                for (int sigma = 0; sigma < 4; sigma++)
                {
                    gamma_value += ginv[lambda][sigma] * (
                        partialg(state,nu, sigma, mu) +
                        partialg(state,mu, sigma, nu) -
                        partialg(state,nu, mu, sigma)
                    );
                }
                gamma[lambda][mu][nu] = 0.5 * gamma_value;
            }
        }
    }
}

State GeodesicEquation(State*& s)
{   
    State out;

    //alloc mem for gamma
    double*** gamma = new double**[4]; 

    for (int i = 0; i < 4; ++i) {
        gamma[i] = new double*[4]; 
        for (int j = 0; j < 4; ++j) {
            gamma[i][j] = new double[4];  
        }
    }

    Gamma(s, gamma);
    
    out.pos = (*s).vel;
    
    for (int lambda = 0; lambda < 4; lambda++)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            for (int nu = 0; nu < 4; nu++)
            {
                out.vel[lambda] -= gamma[lambda][mu][nu] * (*s).vel[mu] * (*s).vel[nu];
            }
        }
    }

    return out;
}

void Integrate(State*& s, double dt) {
    /*double3 org = { 0.0, 0.0, 8.0 };
    double3 force = 0.3 * normalize(state.pos - org) / pow(length(org - state.pos), 2);
    state.vel += dt * force;
    state.pos += dt * state.vel;*/ //Old 
    
    State k1 = GeodesicEquation(s);
    /*State h1;
    h1.pos = s.pos + 0.5 * dt * k1.pos;
    h1.vel = s.vel + 0.5 * dt * k1.vel;
    State k2 = GeodesicEquation(h1);
    State h2;
    h2.pos = s.pos + 0.5 * dt * k2.pos;
    h2.vel = s.vel + 0.5 * dt * k2.vel;
    State k3 = GeodesicEquation(h2);
    State h3;
    h3.pos = s.pos + dt * k3.pos;
    h3.vel = s.vel + dt * k3.vel;
    State k4 = GeodesicEquation(h3);*/
    //to switch to rk4, uncomment previous comment and comment out the euler updating in exchange for rk4 updating
    
    (*s).vel += dt * k1.vel;
    (*s).pos += dt * k1.pos; //note that k1.pos is the velocity as k1 is a derivative
    //s.pos = s.pos + (dt / 6.0) * (k1.pos + 2.0 * k2.pos + 2.0 * k3.pos + k4.pos);
    //s.vel = s.vel + (dt / 6.0) * (k1.vel + 2.0 * k2.vel + 2.0 * k3.vel + k4.vel);
}

glm::vec3 raymarch(const Ray& ray, const glm::vec3& cameraVelocity, double speedOfLight, Data*& data) {
    double dt = 0.1;

    (*(*data).s).pos = (glm::dvec4) glm::vec4(0,ray.origin);
    (*(*data).s).vel = (glm::dvec4) glm::vec4(1,ray.direction);

    for (int i = 0; i < MAX_STEPS; i++) {
        glm::dvec3* p = &glm::dvec3((*(*data).s).pos.x,(*(*data).s).pos.y,(*(*data).s).pos.z);
        double dist = sceneSDF(p);
        if (dist < EPSILON) {
            return glm::vec3(1.0f, 1.0f, 1.0f);
        }

        (*(*data).s).pos += (*(*data).s).vel * dt;
    }

    return glm::vec3(1.0f, 0.0f, 0.0f);
}

void renderFrame(int frameNumber, const glm::vec3& cameraVelocity, double speedOfLight) {
    int width = 800;
    int height = 600;

    std::vector<unsigned char> image(width * height * 3);

    glm::vec3 cameraPos(0, 0, 0);  

    Data* data = new Data(); //reusable memory chunk to hold important data

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double u = (x / (double)width) * 2.0f - 1.0f;
            double v = (y / (double)height) * 2.0f - 1.0f;

            glm::vec3 rayDir = glm::normalize(glm::vec3(u, v, 1.0f));

            Ray ray;
            ray.origin = cameraPos;
            ray.direction = rayDir;


            glm::vec3 color = raymarch(ray, cameraVelocity, speedOfLight, data);

            int index = (y * width + x) * 3;
            image[index] = static_cast<unsigned char>(color.r * 255);
            image[index + 1] = static_cast<unsigned char>(color.g * 255);
            image[index + 2] = static_cast<unsigned char>(color.b * 255);
        }
    }

    std::ofstream file("../../output/frame" + std::to_string(frameNumber) + ".ppm");
    file << "P6\n" << width << " " << height << "\n255\n";
    file.write(reinterpret_cast<char*>(&image[0]), image.size());
    file.close();
}
