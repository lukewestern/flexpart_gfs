import netCDF4 as nc
import numpy as np
import os

dirs = os.listdir('./')

paths=[]
for dd in dirs:
  print(dd)
  if dd.startswith('output_'):
    files_lv1 = os.listdir('./'+dd)
    for f1 in files_lv1:
      print(f1)
      if f1.endswith('_init.nc'):
        paths.append('./'+dd+f1)

print(paths)

ndata_1 = nc.Dataset("./output_omp1/partoutput_20090101000000_init.nc", 'r')
ndata_32 = nc.Dataset("./output_omp32/partoutput_20090101000000_init.nc", 'r')

dvars = list(ndata_1.variables)
for idvars in dvars:
  temp = ndata_1[idvars][:]-ndata_32[idvars][:]
  mn=np.mean(np.abs(temp))
  mx=np.max(np.abs(temp))
  mean1=(np.mean(np.abs(ndata_1[idvars])))
  print(idvars,": ","mean difference:",mn," max abs difference:",mx, "mean: ", mean1)
  if (mn>(mean1/1.e8)):
    print('Mean difference exceeds allowed value of ',mean1/1.e8)
    raise ValueError('Mean allowed difference exceeded.')
ndata_1.close()
ndata_32.close()

ndata_1 = nc.Dataset("./output_omp1/grid_conc_20090101000000.nc", 'r')
ndata_32 = nc.Dataset("./output_omp32/grid_conc_20090101000000.nc", 'r')

print("Checking gridded output...")
#Concentrations
temp=ndata_1['spec001_mr'][:]-ndata_32['spec001_mr'][:]
mn_c=np.mean(temp)
mx_c=np.max(temp)
mn_c1=np.mean(ndata_1['spec001_mr'][:])
if (mn_c>(mn_c1/1.e8)):
  print('Gridded concentration output: mean difference exceeds allowed value of ',mn_c1/1.e8)
  raise ValueError('Gridded concentration output deviates between serial and OMP version')

#Wet deposition
temp=ndata_1['WD_spec001_mr'][:]-ndata_32['WD_spec001_mr'][:]
mn_c=np.mean(temp)
mx_c=np.max(temp)
mn_c1=np.mean(ndata_1['spec001_mr'][:])
if (mn_c>(mn_c1/1.e8)):
  print('Gridded WET deposition output: mean difference exceeds allowed value of ',mn_c1/1.e8)
  raise ValueError('Gridded WET deposition output deviates between serial and OMP version')

#Dry deposition
temp=ndata_1['DD_spec001_mr'][:]-ndata_32['DD_spec001_mr'][:]
mn_c=np.mean(temp)
mx_c=np.max(temp)
mn_c1=np.mean(ndata_1['spec001_mr'][:])
if (mn_c>(mn_c1/1.e8)):
  print('Gridded DRY deposition output: mean difference exceeds allowed value of ',mn_c1/1.e8)
  raise ValueError('Gridded DRY deposition output deviates between serial and OMP version')

ndata_1.close()
ndata_32.close()