# 1. Environment setting
- MPI  
conda install conda-forge::mpich  
- Gfortran  
conda install conda-forge::gfortran
  
- Do not install libraries below. (Version discrepancy) Use the paths of packages installed in Nimbus  
 #Set specific paths to install both packages below:    
 #- ecCodes (Version 2.27.0)  
 #https://confluence.ecmwf.int/display/ECC/Releases  
 #conda install -p path/to/myenv conda-forge::eccodes  
 #- NetCDF  
 #conda install -p path/to/myenv conda-forge::netcdf4  

# 2. Build program  
- Set the paths of ecCodes and NetCDF in 'buildit_mpi' (Already set)
- Execute 'buildit_mpi' in the 'build' directory in the project
  
# 3. Run program  
- Command line : ${path_of_program}/ensalign ${YYYYMMDDHH} ${fhr} ${input_file_path}
- Command line (with MPI): mpirun -np ${number_of_processors} ${path_of_program}/ensalign ${YYYYMMDDHH} ${fhr} ${input_file_path}

# 4. Verification
