#include "solver.hpp"
#include "physics.hpp"
#include "io.hpp"

#include <filesystem>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>

static std::string prepare_output_dir(){
    namespace fs = std::filesystem;
    fs::path base("Result");
    if(fs::exists(base) && !fs::is_empty(base)){
        auto ts = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        fs::rename(base, "Result_"+std::to_string(ts));
    }
    fs::create_directory(base);
    return "Result";
}

int main(){
    const int nx=64, ny=64;
    const double Lx=1.0,Ly=1.0, dx=Lx/(nx-1), dy=Ly/(ny-1);
    const double nu=0.01;
    const int max_steps=300;
    const int output_every=20;

    std::string out_dir = prepare_output_dir();

    AMRGrid amr(nx,ny,Lx,Ly,1);
    std::vector<FlowField> flows;
    flows.emplace_back(nx,ny,dx,dy);
    initialize_MHD_disk(flows[0]); // deterministic seed default

    auto t0=std::chrono::high_resolution_clock::now();
    for(int step=0; step<=max_steps; ++step){
        double max_speed=0.0;
#pragma omp parallel for reduction(max:max_speed)
        for(int i=0;i<nx;++i)
            for(int j=0;j<ny;++j){
                double cs=std::sqrt(1.4*flows[0].p.data[i][j]/flows[0].rho.data[i][j]);
                double spd=std::abs(flows[0].u.data[i][j])+std::abs(flows[0].v.data[i][j])+cs;
                max_speed=std::max(max_speed,spd);
            }
        double dt=0.4*std::min(dx,dy)/(max_speed+1e-12);

        solve_MHD(amr,flows,dt,nu,0,0.0);

        if(step%output_every==0){
            std::cout << "step "<< std::setw(4) << step
                      << " dt="<<dt<<"\n";
            save_flow_MHD(flows[0],out_dir,step);
        }
    }
    auto t1=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;
    std::cout<<"Total time "<<elapsed.count()<<" s\n";
    return 0;
}
