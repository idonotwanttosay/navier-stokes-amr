#pragma once
#include "grid.hpp"
void initialize_MHD_disk(FlowField& flow, int seed = 12345);
void add_divergence_error(FlowField& flow, double amplitude = 0.1); 
