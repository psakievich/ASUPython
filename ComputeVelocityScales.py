
import MrVtkImaginay as MIV
import numpy as np
from mpi4py import MPI
import os


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

num_tsteps = 1024
num_waves = 1019

mean_squared = np.zeros(num_waves)
volume = np.pi*3.15**2

file_template = os.environ["SCRATCH"]."FFT_Results/Snaps/{wave}/SymS1_{wave}_TSTEP_{tstep}.vts"

nFiles = num_waves * num_tsteps

def ComputeMeanSquared(fileName):
    grid = MIV.MrVtkVecHandle(fileName).get()
    mean_squared = MIV.point_product(grid,grid.complex_conjugate())
    values = mean_squared.integrated_values()
    return values[9]

for index in range(rank,nFiles,size):
    my_wave = index%num_tsteps
    my_tstep = index//num_tsteps
    f = file_template.format(wave=my_wave, tstep = my_tstep)
    val = ComputeMeanSquared(f)
    if my_wave is not 0:
        val *= 2.0
    mean_squared[my_wave] += val

mean_squared /= volume
mean_squared /= num_tsteps

dummy = mean_squared.copy()

if rank == 0:
    comm.Reduce(dummy, mean_squared, op = MPI.SUM, root = 0)
else:
    comm.Reduce(mean_squared, dummy, op = MPI.SUM, root = 0)

if rank == 0:
    rms = np.sqrt(mean_squared)
    np.savetxt("rmsVerticalVelocityPerWave.txt",rms,delimiter=',')