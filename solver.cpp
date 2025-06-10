#include "solver.hpp"
#include <cmath>
#include <iostream>
#include <limits>

// Flow field initialization
FlowField::FlowField(int nx, int ny, double dx, double dy, double x0, double y0)
    : u(nx, ny, dx, dy, x0, y0), v(nx, ny, dx, dy, x0, y0), p(nx, ny, dx, dy, x0, y0) {}

// Compute source term for pressure Poisson equation
void compute_source(Grid& b, const FlowField& flow, double dt) {
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < b.nx - 1; ++i) {
        for (int j = 1; j < b.ny - 1; ++j) {
            double dudx = (flow.u.data[i + 1][j] - flow.u.data[i - 1][j]) / (2 * b.dx);
            double dvdy = (flow.v.data[i][j + 1] - flow.v.data[i][j - 1]) / (2 * b.dy);
            b.data[i][j] = (dudx + dvdy) / dt;
        }
    }
}

// Solve pressure Poisson equation
void solve_poisson(Grid& p, const Grid& b, int max_iter, double tolerance) {
    Grid temp(p.nx, p.ny, p.dx, p.dy, p.x0, p.y0);
    double error = 1.0;
    int iter = 0;
    while (error > tolerance && iter < max_iter) {
        error = 0.0;
        #pragma omp parallel for collapse(2) reduction(max:error)
        for (int i = 1; i < p.nx - 1; ++i) {
            for (int j = 1; j < p.ny - 1; ++j) {
                temp.data[i][j] = 0.25 * (p.data[i + 1][j] + p.data[i - 1][j] +
                                          p.data[i][j + 1] + p.data[i][j - 1] -
                                          p.dx * p.dx * b.data[i][j]);
                error = std::max(error, std::abs(temp.data[i][j] - p.data[i][j]));
            }
        }
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < p.nx - 1; ++i) {
            for (int j = 1; j < p.ny - 1; ++j) {
                p.data[i][j] = temp.data[i][j];
                // Check for NaN or infinity
                if (!std::isfinite(p.data[i][j])) {
                    std::cerr << "Error: Non-finite pressure detected at i=" << i << ", j=" << j << std::endl;
                    exit(1);
                }
            }
        }
        iter++;
        std::cout << "Poisson iteration " << iter << ": error=" << error << std::endl;
    }
}

// Navier-Stokes solver with AMR
void solve_navier_stokes(AMRGrid& amr, FlowField& flow, double dt, double nu, int max_iter, double tolerance) {
    Grid b(flow.p.nx, flow.p.ny, flow.p.dx, flow.p.dy, flow.p.x0, flow.p.y0);
    
    // Update velocity with viscous and nonlinear terms
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < flow.u.nx - 1; ++i) {
        for (int j = 1; j < flow.u.ny - 1; ++j) {
            double lap_u = (flow.u.data[i + 1][j] + flow.u.data[i - 1][j] +
                           flow.u.data[i][j + 1] + flow.u.data[i][j - 1] -
                           4 * flow.u.data[i][j]) / (flow.u.dx * flow.u.dx);
            double dpdx = (flow.p.data[i + 1][j] - flow.p.data[i - 1][j]) / (2 * flow.p.dx);
            double dudx = (flow.u.data[i + 1][j] - flow.u.data[i - 1][j]) / (2 * flow.u.dx);
            double dudy = (flow.u.data[i][j + 1] - flow.u.data[i][j - 1]) / (2 * flow.u.dy);
            double conv_u = flow.u.data[i][j] * dudx + flow.v.data[i][j] * dudy;
            flow.u.data[i][j] += dt * (nu * lap_u - dpdx - conv_u);
        }
    }
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < flow.v.nx - 1; ++i) {
        for (int j = 1; j < flow.v.ny - 1; ++j) {
            double lap_v = (flow.v.data[i + 1][j] + flow.v.data[i - 1][j] +
                           flow.v.data[i][j + 1] + flow.v.data[i][j - 1] -
                           4 * flow.v.data[i][j]) / (flow.v.dx * flow.v.dx);
            double dpdy = (flow.p.data[i][j + 1] - flow.p.data[i][j - 1]) / (2 * flow.p.dy);
            double dvdx = (flow.v.data[i + 1][j] - flow.v.data[i - 1][j]) / (2 * flow.v.dx);
            double dvdy = (flow.v.data[i][j + 1] - flow.v.data[i][j - 1]) / (2 * flow.v.dy);
            double conv_v = flow.u.data[i][j] * dvdx + flow.v.data[i][j] * dvdy;
            flow.v.data[i][j] += dt * (nu * lap_v - dpdy - conv_v);
        }
    }

    // Compute pressure source and solve
    compute_source(b, flow, dt);
    solve_poisson(flow.p, b, max_iter, tolerance);

    // Dynamic grid refinement with single-threaded loop
    for (int i = 5; i < flow.p.nx - 5; ++i) {
        for (int j = 5; j < flow.p.ny - 5; ++j) {
            if (amr.needs_refinement(flow.p, i, j, 5.0)) { // Increased threshold
                std::cout << "Refining at i=" << i << ", j=" << j << std::endl;
                amr.refine(0, i, j, 12, 12); // Smaller fine grid
            }
        }
    }
}