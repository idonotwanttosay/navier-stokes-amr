# Navier-Stokes AMR Solver

This project contains a simple 2-D magneto-hydrodynamics solver with adaptive mesh refinement (AMR).

## Requirements

- **C++ compiler** supporting C++17 and OpenMP (e.g. `g++` 10 or later). The provided `run.sh` expects a compiler named `g++-14`; adjust the script if your compiler has a different name.
- **Python 3** with the packages:
  - `numpy`
  - `pandas`
  - `matplotlib`
  FFMPEG is recommended for creating the animation.

## Building and Running

Execute the helper script `run.sh` from the repository root:

```bash
bash run.sh
```

The script compiles all `*.cpp` sources into the executable `ns2d` and runs the solver.  The environment variable `OMP_NUM_THREADS` controls the number of OpenMP threads (default is 4).
After the simulation finishes, several Python scripts are executed to post-process the results and produce plots and an optional animation.

## Output

All results are written into a directory named `Result`.  For each output step CSV files are created:

- `out_rho_<step>.csv`, `out_u_<step>.csv`, `out_v_<step>.csv`, `out_e_<step>.csv`, `out_bx_<step>.csv`, `out_by_<step>.csv` – fluid variables on the grid.
- `grid_<step>.csv` – AMR grid metadata.

The post-processing scripts generate additional files inside the same directory:

- `density_<step>.png` – contour plot of the density field.
- `conservation.png`, `divB_error.png`, `energy_spectrum.png`, `summary_mass_energy.png`, `summary_divB.png` – analysis figures.
- `animation.mp4` – optional animation of the density and grid evolution.

Previous `Result` directories are automatically renamed with a timestamp when a new run starts.

