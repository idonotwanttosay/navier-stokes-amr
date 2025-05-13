#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "grid.hpp"

// Flow field and Navier-Stokes solver
struct FlowField {
    Grid u, v, p;
    FlowField(int nx, int ny, double dx, double dy, double x0, double y0);
};

void solve_navier_stokes(AMRGrid& amr, FlowField& flow, double dt, double nu, int max_iter, double tolerance);
void compute_source(Grid& b, const FlowField& flow, double dt);

#endif