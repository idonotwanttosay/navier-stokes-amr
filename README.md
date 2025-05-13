g++-14 -fopenmp *.cpp -o navier_stokes_amr   


export OMP_NUM_THREADS=4                                  
./navier_stokes_amr > output.txt 

python -u "/Users/liuyuen/project/plot_flow.py"           
