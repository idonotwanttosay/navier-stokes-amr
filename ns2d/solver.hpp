#pragma once
#include "grid.hpp"
void solve_MHD(AMRGrid& amr, FlowField& flow,double dt,double nu,int max_iter,double tol);
