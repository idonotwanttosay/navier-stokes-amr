#include "io.hpp"
#include <fstream>
#include <iomanip>

// Save flow field to CSV
void save_flow_field(const FlowField& flow, const std::string& prefix, int step) {
    std::ofstream out_p(prefix + "_p_" + std::to_string(step) + ".csv");
    std::ofstream out_u(prefix + "_u_" + std::to_string(step) + ".csv");
    std::ofstream out_v(prefix + "_v_" + std::to_string(step) + ".csv");
    for (int i = 0; i < flow.p.nx; ++i) {
        for (int j = 0; j < flow.p.ny; ++j) {
            double x = flow.p.x0 + i * flow.p.dx;
            double y = flow.p.y0 + j * flow.p.dy;
            out_p << x << "," << y << "," << flow.p.data[i][j] << "\n";
            out_u << x << "," << y << "," << flow.u.data[i][j] << "\n";
            out_v << x << "," << y << "," << flow.v.data[i][j] << "\n";
        }
    }
}