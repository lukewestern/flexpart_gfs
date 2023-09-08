import netCDF4 as nc
import numpy as np
import os
import re

return_flag = False

resolution=1.e7

compare = False
# Comparing to previous version, if it exists
if (os.path.isfile("./bkw_master.txt")):
  compare = True
  fr = open("./bkw_master.txt", "r")

output_name = 'bkw_test.txt'
with open(output_name, 'a') as f:

  ndata = nc.Dataset("./output_bkw/grid_drydep_20090101030000.nc", 'r')

  print("Checking gridded output...")
  #Concentrations
  temp=np.abs(ndata['spec001_mr'][:,:,1,:,:,:])
  mn_c=np.mean(temp)
  mx_c=np.max(temp)
  f.write(" Dry:\t\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn_c,mx_c))
  if compare:
    tmp=fr.readline()
    tmp2=tmp.split()
    mn_prev=float(tmp2[3])
    mx_prev=float(tmp2[6])
    if (np.abs(mn_c-mn_prev)>(np.abs(mn_c)/resolution)): 
      print("WARNING: CONCENTRATIONS give a different mean value compared to the previous version.")
      print(np.abs(mn_c-mn_prev)/np.abs(mn_c)*100, "percent")
      return_flag = True
    if (np.abs(mx_c-mx_prev)>(np.abs(mx_c)/resolution)):
      print("WARNING: CONCENTRATIONS give a different max value compared to the previous version.")
      print(np.abs(mx_c-mx_prev)/np.abs(mx_c)*100, "percent")
      return_flag = True
    #f.write(tmp)
    #f.write("\n")
  ndata.close()

  ndata = nc.Dataset("./output_bkw_eta/grid_drydep_20090101030000.nc", 'r')
    
  print("Checking gridded output...")
  #Concentrations
  temp=np.abs(ndata['spec001_mr'][:,:,1,:,:,:])
  mn_c=np.mean(temp)
  mx_c=np.max(temp)
  f.write(" Dry ETA:\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn_c,mx_c))
  if compare:
    tmp=fr.readline()
    tmp2=tmp.split()
    mn_prev=float(tmp2[4])
    mx_prev=float(tmp2[7])
    if (np.abs(mn_c-mn_prev)>(np.abs(mn_c)/resolution)): 
      print("WARNING: CONCENTRATIONS give a different mean value compared to the previous version.")
      print(np.abs(mn_c-mn_prev)/np.abs(mn_c)*100, "percent")
      return_flag = True
    if (np.abs(mx_c-mx_prev)>(np.abs(mx_c)/resolution)):
      print("WARNING: CONCENTRATIONS give a different max value compared to the previous version.")
      print(np.abs(mx_c-mx_prev)/np.abs(mx_c)*100, "percent")
      return_flag = True
    #f.write(tmp)
    #f.write("\n")
  ndata.close()

  ndata = nc.Dataset("./output_bkw/grid_wetdep_20090101030000.nc", 'r')
  
  print("Checking gridded output...")
  #Concentrations
  temp=np.abs(ndata['spec001_mr'][:,:,1,:,:,:])
  mn_c=np.mean(temp)
  mx_c=np.max(temp)
  f.write(" Wet:\t\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn_c,mx_c))
  if compare:
    tmp=fr.readline()
    tmp2=tmp.split()
    mn_prev=float(tmp2[3])
    mx_prev=float(tmp2[6])
    if (np.abs(mn_c-mn_prev)>(np.abs(mn_c)/resolution)): 
      print("WARNING: CONCENTRATIONS give a different mean value compared to the previous version.")
      print(np.abs(mn_c-mn_prev)/np.abs(mn_c)*100, "percent")
      return_flag = True
    if (np.abs(mx_c-mx_prev)>(np.abs(mx_c)/resolution)):
      print("WARNING: CONCENTRATIONS give a different max value compared to the previous version.")
      print(np.abs(mx_c-mx_prev)/np.abs(mx_c)*100, "percent")
      return_flag = True
    #f.write(tmp)
    #f.write("\n")
  ndata.close()

  ndata = nc.Dataset("./output_bkw_eta/grid_wetdep_20090101030000.nc", 'r')
  
  print("Checking gridded output...")
  #Concentrations
  temp=np.abs(ndata['spec001_mr'][:,:,1,:,:,:])
  mn_c=np.mean(temp)
  mx_c=np.max(temp)
  f.write(" Wet ETA:\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn_c,mx_c))
  if compare:
    tmp=fr.readline()
    tmp2=tmp.split()
    mn_prev=float(tmp2[4])
    mx_prev=float(tmp2[7])
    if (np.abs(mn_c-mn_prev)>(np.abs(mn_c)/resolution)): 
      print("WARNING: CONCENTRATIONS give a different mean value compared to the previous version.")
      print(np.abs(mn_c-mn_prev)/np.abs(mn_c)*100, "percent")
      return_flag = True
    if (np.abs(mx_c-mx_prev)>(np.abs(mx_c)/resolution)):
      print("WARNING: CONCENTRATIONS give a different max value compared to the previous version.")
      print(np.abs(mx_c-mx_prev)/np.abs(mx_c)*100, "percent")
      return_flag = True
    #f.write(tmp)
    #f.write("\n")
  ndata.close()

if (return_flag):
  raise ValueError('Old and New version deviate')
elif (compare):
  print("Output from previous version is not significantly different from this version.")
