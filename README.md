g++-14 -fopenmp poisson_solver.cpp -o poisson_solver
export OMP_NUM_THREADS=4                                  
./navier_stokes_amr > output.txt
python -u "/Users/liuyuen/project/plot_flow.py"    
