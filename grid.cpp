#include "grid.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

// Grid initialization and data management
Grid::Grid(int nx_, int ny_, double dx_, double dy_, double x0_, double y0_)
    : nx(nx_), ny(ny_), dx(dx_), dy(dy_), x0(x0_), y0(y0_) {
    data.resize(nx, std::vector<double>(ny, 0.0));
}

void Grid::fill(double value) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            data[i][j] = value;
        }
    }
}

// AMR grid refinement
AMRGrid::AMRGrid(int base_nx, int base_ny, double Lx, double Ly, int max_level_)
    : max_level(max_level_) {
    levels.emplace_back(base_nx, base_ny, Lx / (base_nx - 1), Ly / (base_ny - 1), 0.0, 0.0);
}

void AMRGrid::refine(int level, int start_x, int start_y, int fine_nx, int fine_ny) {
    if (level >= max_level || refinement_count >= MAX_REFINEMENTS) {
        std::cout << "Refinement skipped: max_level=" << max_level 
                  << ", refinement_count=" << refinement_count << std::endl;
        return;
    }
    
    // Ensure valid refinement region
    start_x = std::max(0, std::min(start_x, levels[level].nx - fine_nx / 2));
    start_y = std::max(0, std::min(start_y, levels[level].ny - fine_ny / 2));
    double fine_dx = levels[level].dx / 2.0;
    double fine_dy = levels[level].dy / 2.0;
    double x0 = levels[level].x0 + start_x * levels[level].dx;
    double y0 = levels[level].y0 + start_y * levels[level].dy;

    // Diagnostic output
    std::cout << "Creating fine grid: start_x=" << start_x << ", start_y=" << start_y 
              << ", fine_nx=" << fine_nx << ", fine_ny=" << fine_ny 
              << ", refinement_count=" << refinement_count + 1 << std::endl;

    // Create new fine grid
    levels.emplace_back(fine_nx, fine_ny, fine_dx, fine_dy, x0, y0);
    refinement_count++; // Increment count

    Grid& fine = levels.back();
    Grid& coarse = levels[level];

    // Interpolate coarse to fine grid
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < fine_nx; ++i) {
        for (int j = 0; j < fine_ny; ++j) {
            int coarse_i = std::min(start_x + i / 2, coarse.nx - 1);
            int coarse_j = std::min(start_y + j / 2, coarse.ny - 1);
            fine.data[i][j] = coarse.data[coarse_i][coarse_j];
        }
    }
}

// Refinement criterion based on gradient
bool AMRGrid::needs_refinement(const Grid& grid, int i, int j, double threshold) {
    if (i <= 0 || i >= grid.nx - 1 || j <= 0 || j >= grid.ny - 1) return false;
    double grad_x = std::abs(grid.data[i + 1][j] - grid.data[i - 1][j]) / (2 * grid.dx);
    double grad_y = std::abs(grid.data[i][j + 1] - grid.data[i][j - 1]) / (2 * grid.dy);
    double gradient = grad_x + grad_y;
    std::cout << "Gradient at i=" << i << ", j=" << j << ": " << gradient << std::endl;
    return gradient > threshold;
}