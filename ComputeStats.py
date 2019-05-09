'''
Compute statistics using a mean value in a single
plane and decomposed instantaneous fileds
'''

from mpi4py import MPI
import MrImaginaryVtk as MIV
import MrRealVtk as MRV
import ParallelMapping as pm
import os

'''
Gather all values of interest
flucs^2
flucs^3
flucs^4
number of values at each rz
compute final QOI
write files
'''

mean_file_path=os.environ["SCRATCH"]+'/FFT_Results/Snaps/0/SymWave_0_TAVG_0-1017.vts'
num_tsteps = 1
num_writing_procs = 3
file_template = "{tstep}/RSnap_{wrank}_TSTEP_{tstep}.vts"

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    global_blocks = (1,20,64)
    
    # all ranks read the mean value
    mean = MIV.MrVtkVecHandle(mean_file_path).get()
    
    mdims = mean.data.GetDimensions()
    
    # switch index to match real data
    mean.data.SetDimensions((mdims[1],mdims[0],mdims[2]))
    
    variance = mean*0.0
    skewness = mean*0.0
    kurtosis = mean*0.0
    
    # each time step divide files amongst processors
    for t in range(num_tsteps):
        for writeRank in range(rank, num_writing_procs, size):
            inst = MRV.MrVtkVecHandle(file_template.format(wrank=writeRank,
                                                           tstep=t)).get()
            v = inst.data.GetPointData().GetArray(2)
            p = inst.data.GetPoints()
            pm.Cart2Cyl(p,v)
            
            meanL,c0,c1 = pm.CreateLocalMeanMr(mean,inst,writeRank,
                                                  global_blocks)
            fluc = inst - meanL
            varL = pm.AzimuthalAverage(fluc.power(2.0),writeRank,global_blocks)
            skeL = pm.AzimuthalAverage(fluc.power(3.0),writeRank,global_blocks)
            kurL = pm.AzimuthalAverage(fluc.power(4.0),writeRank,global_blocks)
            
            pm.AddSubset(variance,varL,c0,c1)
            pm.AddSubset(skewness,skeL,c0,c1)
            pm.AddSubset(kurtosis,kurL,c0,c1)
    
    # time average               
    variance*=1.0/num_tsteps
    skewness*=1.0/num_tsteps
    kurtosis*=1.0/num_tsteps
    
    # parallel assembly of fields
    pm.ReduceGridToRank0(variance,comm)
    pm.ReduceGridToRank0(skewness,comm)
    pm.ReduceGridToRank0(kurtosis,comm)
    
    if rank == 0:
        # normalize skewness and kurtosis
        skewness = MIV.point_division(skewness,variance.power(1.5))
        kurtosis = MIV.point_division(kurtosis,variance.power(2.0))
        
        variance.data.SetDimensions(mdims)
        skewness.data.SetDimensions(mdims)
        kurtosis.data.SetDimensions(mdims)
        
        # write result files
        MIV.MrVtkVecHandle("variance.vts").put(variance)
        MIV.MrVtkVecHandle("skewness.vts").put(skewness)
        MIV.MrVtkVecHandle("kurtosis.vts").put(kurtosis)
    
    
            
            
            
    