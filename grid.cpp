#include "grid.hpp"
#include <algorithm>
#include <omp.h>
#include <iostream>

int refinement_count = 0;

Grid::Grid(int nx_,int ny_,double dx_,double dy_,double x0_,double y0_)
    : nx(nx_),ny(ny_),dx(dx_),dy(dy_),x0(x0_),y0(y0_),
      data(nx_, std::vector<double>(ny_,0.0)) {}

void Grid::fill(double v){
#pragma omp parallel for collapse(2)
    for(int i=0;i<nx;++i)
        for(int j=0;j<ny;++j)
            data[i][j]=v;
}

AMRGrid::AMRGrid(int base_nx,int base_ny,double Lx,double Ly,int max_level_)
    : max_level(max_level_)
{
    levels.emplace_back(base_nx,base_ny,Lx/(base_nx-1),Ly/(base_ny-1),0.0,0.0);
}

void AMRGrid::refine(int level,int start_x,int start_y,int fine_nx,int fine_ny){
    if(level>=max_level || refinement_count>=MAX_REFINEMENTS) return;
    start_x=std::max(0,std::min(start_x,levels[level].nx-fine_nx/2));
    start_y=std::max(0,std::min(start_y,levels[level].ny-fine_ny/2));
    double fine_dx = levels[level].dx/2.0;
    double fine_dy = levels[level].dy/2.0;
    double x0=levels[level].x0 + start_x*levels[level].dx;
    double y0=levels[level].y0 + start_y*levels[level].dy;
    std::cout<<"[AMR] refine @"<<level<<" ("<<start_x<<","<<start_y<<")\n";
    levels.emplace_back(fine_nx,fine_ny,fine_dx,fine_dy,x0,y0);
    ++refinement_count;
    // injection
    Grid& fine = levels.back();
    Grid& coarse = levels[level];
#pragma omp parallel for collapse(2)
    for(int i=0;i<fine_nx;++i)
        for(int j=0;j<fine_ny;++j){
            int ci = std::min(start_x+i/2,coarse.nx-1);
            int cj = std::min(start_y+j/2,coarse.ny-1);
            fine.data[i][j]=coarse.data[ci][cj];
        }
}

bool AMRGrid::needs_refinement(const Grid& g,int i,int j,double thres){
    if(i<=0||i>=g.nx-1||j<=0||j>=g.ny-1) return false;
    double gx=std::fabs(g.data[i+1][j]-g.data[i-1][j])/(2*g.dx);
    double gy=std::fabs(g.data[i][j+1]-g.data[i][j-1])/(2*g.dy);
    return (gx+gy)>thres;
}

FlowField::FlowField(int nx,int ny,double dx,double dy)
    : rho(nx,ny,dx,dy), u(nx,ny,dx,dy), v(nx,ny,dx,dy),
      p(nx,ny,dx,dy), e(nx,ny,dx,dy),
      bx(nx,ny,dx,dy), by(nx,ny,dx,dy), psi(nx,ny,dx,dy) {}
