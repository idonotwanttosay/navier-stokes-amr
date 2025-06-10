#ifndef IO_HPP
#define IO_HPP

#include "solver.hpp"
#include <string>

// Data input/output
void save_flow_field(const FlowField& flow, const std::string& prefix, int step);

#endif