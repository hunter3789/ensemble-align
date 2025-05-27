# Description

https://doi.org/10.1175/WAF-D-23-0229.1

# 2023 FFaIR Summer Experiment results
https://caps.ou.edu/clee/ens/ens_view.php

# 1. Required packages/libraries  
- MPI  
- Gfortran  
- ecCodes (Version 2.27.0), for GRIB I/O
- NetCDF

# 2. Build program  
- Set the paths of ecCodes and NetCDF in 'buildit_mpi' (Already set)
- Execute 'buildit_mpi' in the 'build' directory in the project
  
# 3. Run program  
- Command line : ${path_of_program}/ensalign ${YYYYMMDDHH} ${fhr} ${input_file_path}
- Command line (with MPI): mpirun -np ${number_of_processors} ${path_of_program}/ensalign ${YYYYMMDDHH} ${fhr} ${input_file_path}

# 4. Output file
**4.1) command line to read the header of output file**  
- ncdump -h ${filename}  

**4.2) variables**  
- ensfcst(nmembers, y, x) : Each member's original preicipitation field from ensemble system  
- ensshffcst(nshfpass, nmembers, y, x) : Each member's aligned precipitation field    
- shiftvec_u(nshfpass, nmembers, y, x), shiftvec_v(nshfpass, nmembers, y, x) : Each member's shift vector fields  
- ensmean(y, x) : Ensemble mean of precipitation  
- enspm(y, x) : Ensemble PM mean of precipitation  
- enslpm(y, x) : Ensemble LPM mean of precipitation  
- ensshfmean(nshfpass, y, x) : Spatially Aligned Mean (SAM) of precipitation  
- ensshflpm(nshfpass, y, x) : Spatially Aligned Mean with LPM method applied (SAM-LPM) of precipitation  
