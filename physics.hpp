#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "grid.hpp"
#include "solver.hpp"

// Astrophysical scenario setup
void initialize_accretion_disk(FlowField& flow, double G, double M);
void apply_gravity(FlowField& flow, double G, double M, double dt);

#endif