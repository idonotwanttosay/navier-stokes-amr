#include "physics.hpp"
#include <cmath>

// Initialize accretion disk conditions
void initialize_accretion_disk(FlowField& flow, double G, double M) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < flow.u.nx; ++i) {
        for (int j = 0; j < flow.u.ny; ++j) {
            double x = flow.u.x0 + i * flow.u.dx - 0.5;
            double y = flow.u.y0 + j * flow.u.dy - 0.5;
            double r = std::sqrt(x * x + y * y);
            double v_theta = std::sqrt(G * M / std::max(r, 0.01)); // Avoid division by zero
            flow.u.data[i][j] = -v_theta * y / std::max(r, 0.01);
            flow.v.data[i][j] = v_theta * x / std::max(r, 0.01);
            flow.p.data[i][j] = 0.01 / (r * r); // Reduced pressure magnitude
        }
    }
}

// Apply gravitational force
void apply_gravity(FlowField& flow, double G, double M, double dt) {
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < flow.u.nx - 1; ++i) {
        for (int j = 1; j < flow.u.ny - 1; ++j) {
            double x = flow.u.x0 + i * flow.u.dx - 0.5;
            double y = flow.u.y0 + j * flow.u.dy - 0.5;
            double r = std::sqrt(x * x + y * y);
            double ax = -G * M * x / (std::max(r, 0.01) * std::max(r, 0.01) * std::max(r, 0.01));
            double ay = -G * M * y / (std::max(r, 0.01) * std::max(r, 0.01) * std::max(r, 0.01));
            flow.u.data[i][j] += dt * ax;
            flow.v.data[i][j] += dt * ay;
        }
    }
}