#include "solver.hpp"
#include "physics.hpp"
#include "io.hpp"
#include <iostream>
#include <chrono>

// Main simulation loop
int main() {
    // Simulation parameters
    const int nx = 40, ny = 40;
    const double Lx = 1.0, Ly = 1.0;
    const double dx = Lx / (nx - 1), dy = Ly / (ny - 1);
    const double dt = 0.001, nu = 0.01;
    const int max_steps = 50;
    const int max_iter = 200;
    const double tolerance = 1e-3;
    const double G = 1.0, M = 1.0;

    // Initialize grid, flow, and AMR
    AMRGrid amr(nx, ny, Lx, Ly, 2);
    FlowField flow(nx, ny, dx, dy, 0.0, 0.0);
    initialize_accretion_disk(flow, G, M);

    // Timing
    auto start = std::chrono::high_resolution_clock::now();

    // Main loop
    for (int step = 0; step < max_steps; ++step) {
        std::cout << "Step: " << step << std::endl;
        solve_navier_stokes(amr, flow, dt, nu, max_iter, tolerance);
        if (step % 2 == 0) { // Save every 2 steps
            save_flow_field(flow, "output", step);
        }
    }

    // Output simulation time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Simulation time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}