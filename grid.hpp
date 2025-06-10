#ifndef GRID_HPP
#define GRID_HPP

#include <vector>

// Grid and AMR structure definitions
struct Grid {
    int nx, ny;
    double dx, dy, x0, y0;
    std::vector<std::vector<double>> data;
    Grid(int nx_, int ny_, double dx_, double dy_, double x0_, double y0_);
    void fill(double value);
};

struct AMRGrid {
    std::vector<Grid> levels;
    int max_level;
    int refinement_count = 0; // Track refinement count
    static const int MAX_REFINEMENTS = 5; // Limit refinements
    AMRGrid(int base_nx, int base_ny, double Lx, double Ly, int max_level_);
    void refine(int level, int start_x, int start_y, int fine_nx, int fine_ny);
    bool needs_refinement(const Grid& grid, int i, int j, double threshold);
};

#endif