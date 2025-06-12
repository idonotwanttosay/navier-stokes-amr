#pragma once
#include "grid.hpp"
#include <vector>
void solve_MHD(AMRGrid& amr, std::vector<FlowField>& flows,double dt,double nu,int max_iter,double tol);
