#include "solver.hpp"
#include <omp.h>
#include <cmath>
#include <vector>

static constexpr double ETA=0.01;
static constexpr double CH =1.0;
static constexpr double CR =0.1;
static constexpr double gamma_gas=1.4;

static inline double laplacian(const Grid& g,int i,int j){
    return (g.data[i+1][j]-2*g.data[i][j]+g.data[i-1][j])/(g.dx*g.dx)
         + (g.data[i][j+1]-2*g.data[i][j]+g.data[i][j-1])/(g.dy*g.dy);
}

void solve_MHD(AMRGrid& amr, FlowField& flow,double dt,double nu,int max_iter, double tol){
    Grid& grid=flow.rho;
    auto rho_new=flow.rho.data;
    auto u_new=flow.u.data;
    auto v_new=flow.v.data;
    auto e_new=flow.e.data;
    auto p_new=flow.p.data;
    auto bx_new=flow.bx.data;
    auto by_new=flow.by.data;
    auto psi_new=flow.psi.data;

#pragma omp parallel for collapse(2)
    for(int i=1;i<grid.nx-1;++i){
        for(int j=1;j<grid.ny-1;++j){
            double rho=flow.rho.data[i][j];
            double u=flow.u.data[i][j];
            double v=flow.v.data[i][j];
            double bx=flow.bx.data[i][j];
            double by=flow.by.data[i][j];
            double psi=flow.psi.data[i][j];

            auto ddx=[&](const Grid& g){ return (g.data[i+1][j]-g.data[i-1][j])/(2*grid.dx); };
            auto ddy=[&](const Grid& g){ return (g.data[i][j+1]-g.data[i][j-1])/(2*grid.dy); };

            // density update
            double div_rho_u = ddx(flow.rho)*u + ddy(flow.rho)*v + rho*(ddx(flow.u)+ddy(flow.v));
            rho_new[i][j]=rho - dt*div_rho_u;

            // momentum primitive update for simplicity
            double adv_u = -u*ddx(flow.u)-v*ddy(flow.u);
            double adv_v = -u*ddx(flow.v)-v*ddy(flow.v);
            double diff_u=nu*laplacian(flow.u,i,j);
            double diff_v=nu*laplacian(flow.v,i,j);
            double p_gradx=-(flow.p.data[i+1][j]-flow.p.data[i-1][j])/(2*grid.dx*rho);
            double p_grady=-(flow.p.data[i][j+1]-flow.p.data[i][j-1])/(2*grid.dy*rho);
            double Jz = (flow.by.data[i+1][j]-flow.by.data[i-1][j])/(2*grid.dx)
                       -(flow.bx.data[i][j+1]-flow.bx.data[i][j-1])/(2*grid.dy);
            double fx=Jz*by/rho, fy=-Jz*bx/rho;

            u_new[i][j]=u+dt*(adv_u+diff_u+p_gradx+fx);
            v_new[i][j]=v+dt*(adv_v+diff_v+p_grady+fy);

            double ke_new=0.5*rho_new[i][j]*(u_new[i][j]*u_new[i][j]+v_new[i][j]*v_new[i][j]);
            p_new[i][j]=flow.p.data[i][j]; // still isothermal
            e_new[i][j]=p_new[i][j]/(gamma_gas-1.0)+ke_new;

            // induction + cleaning
            bx_new[i][j] = bx + dt*( ETA*laplacian(flow.bx,i,j) - (psi_new[i+1][j]-psi_new[i-1][j])/(2*grid.dx) );
            by_new[i][j] = by + dt*( ETA*laplacian(flow.by,i,j) - (psi_new[i][j+1]-psi_new[i][j-1])/(2*grid.dy) );
            double divB=ddx(flow.bx)+ddy(flow.by);
            psi_new[i][j]=psi - dt*(CH*CH*divB + CR*CR*psi);
        }
    }

    flow.rho.data=std::move(rho_new);
    flow.u.data=std::move(u_new);
    flow.v.data=std::move(v_new);
    flow.e.data=std::move(e_new);
    flow.p.data=std::move(p_new);
    flow.bx.data=std::move(bx_new);
    flow.by.data=std::move(by_new);
    flow.psi.data=std::move(psi_new);
}
