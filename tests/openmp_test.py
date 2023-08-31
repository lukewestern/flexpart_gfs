import netCDF4 as nc
import numpy as np
import os

ndata_1 = nc.Dataset("./output_omp1/partoutput_20090101000000_init.nc", 'r')
ndata_32 = nc.Dataset("./output_omp32/partoutput_20090101000000_init.nc", 'r')

dvars = list(ndata_1.variables)
for idvars in dvars:
  temp = ndata_1[idvars][:]-ndata_32[idvars][:]
  print(idvars, ':')
  mn=np.mean(np.abs(temp))
  mx=np.max(np.abs(temp))
  mean1=(np.mean(np.abs(ndata_1[idvars])))
  print("mean difference:",mn," max abs difference:",mx, "mean: ", mean1)

  print("small mean", smallmean)
  if (mn>(mean1/1.e8)):
    print('Mean difference exceeds allowed value of ',mean1/1.e8)
    raise ValueError('Mean allowed difference exceeded.')
ndata_1.close()
ndata_32.close()

# clean up
os.rmdir('./current') 
os.rmdir('./output_omp1')
os.rmdir('./output_omp32')
os.remove('pathnames_omp1')
os.remove('pathnames_omp32')