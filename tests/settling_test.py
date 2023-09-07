import netCDF4 as nc
import numpy as np
import os
import re

return_flag = False

resolution=1.e7

compare = False
# Comparing to previous version, if it exists
if (os.path.isfile("./settling_old.txt")):
  compare = True
  fr = open("./settling_old.txt", "r")

output_name = 'settling_test.txt'
with open(output_name, 'a') as f:
  ndata = nc.Dataset("./output_settling/partoutput_20090101000000_init.nc", 'r')
  ndata2 = nc.Dataset("./output_settling_eta/partoutput_20090101000000_init.nc", 'r')

  dvars = list(ndata.variables)
  for idvars in dvars:
    if (idvars in ['time','particle']):
      continue
    mn=np.mean(np.abs(ndata[idvars][:,-1]))
    mx=np.max(np.abs(ndata[idvars][:,-1]))
    mn2=np.mean(np.abs(ndata2[idvars][:,-1]))
    mx2=np.max(np.abs(ndata2[idvars][:,-1]))
    f.write("%19s:\tmean abs: %2.8e\tmax abs: %2.8e\n" %(idvars,mn,mx))
    f.write("%15s ETA:\tmean abs: %2.8e\tmax abs: %2.8e\n" %(idvars,mn2,mx2))
    if compare:
      tmp=fr.readline()
      tmp2=tmp.split()
      mn_prev=float(tmp2[3])
      mx_prev=float(tmp2[6])
      tmp=fr.readline()
      tmp2=tmp.split()
      mn2_prev=float(tmp2[4])
      mx2_prev=float(tmp2[7])
      if (np.abs(mn-mn_prev)>(np.abs(mn)/resolution)): 
        print("WARNING: ",idvars,"gives a different mean value compared to the previous version.")
        print(np.abs(mn-mn_prev)/np.abs(mn)*100, "percent")
        return_flag = True
      if (np.abs(mn2-mn2_prev)>(np.abs(mn2)/resolution)): 
        print("WARNING: ",idvars,"ETA gives a different mean value compared to the previous version.")
        print(np.abs(mn2-mn2_prev)/np.abs(mn2)*100, "percent")
        return_flag = True
      if (np.abs(mx-mx_prev)>(np.abs(mx)/resolution)):
        print("WARNING: ",idvars,"gives a different max value compared to the previous version.")
        print(np.abs(mx-mx_prev)/np.abs(mx)*100, "percent")
        return_flag = True
      if (np.abs(mx2-mx2_prev)>(np.abs(mx2)/resolution)):
        print("WARNING: ",idvars,"ETA gives a different max value compared to the previous version.")
        print(np.abs(mx2-mx2_prev)/np.abs(mx2)*100, "percent")
        return_flag = True
      #f.write(tmp)
      #f.write("\n")

  ndata.close()
  ndata2.close()

  ndata = nc.Dataset("./output_settling/grid_conc_20090101000000.nc", 'r')
  ndata2 = nc.Dataset("./output_settling_eta/grid_conc_20090101000000.nc", 'r')
  
  print("Checking gridded output...")
  #Concentrations
  temp=np.abs(ndata['spec001_mr'][:,:,1,:,:,:])
  mn_c=np.mean(temp)
  mx_c=np.max(temp)
  temp2=np.abs(ndata2['spec001_mr'][:,:,1,:,:,:])
  mn2_c=np.mean(temp2)
  mx2_c=np.max(temp2)
  f.write(" Concentrations:\t\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn_c,mx_c))
  f.write(" Concentrations ETA:\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn2_c,mx2_c))
  if compare:
    tmp=fr.readline()
    tmp2=tmp.split()
    mn_prev=float(tmp2[3])
    mx_prev=float(tmp2[6])
    tmp=fr.readline()
    tmp2=tmp.split()
    mn2_prev=float(tmp2[4])
    mx2_prev=float(tmp2[7])
    if (np.abs(mn_c-mn_prev)>(np.abs(mn_c)/resolution)): 
      print("WARNING: CONCENTRATIONS give a different mean value compared to the previous version.")
      print(np.abs(mn_c-mn_prev)/np.abs(mn_c)*100, "percent")
      return_flag = True
    if (np.abs(mn2_c-mn2_prev)>(np.abs(mn2_c)/resolution)): 
      print("WARNING: ETA CONCENTRATIONS give a different mean value compared to the previous version.")
      print(np.abs(mn2_c-mn2_prev)/np.abs(mn2_c)*100, "percent")
      return_flag = True
    if (np.abs(mx_c-mx_prev)>(np.abs(mx_c)/resolution)):
      print("WARNING: CONCENTRATIONS give a different max value compared to the previous version.")
      print(np.abs(mx_c-mx_prev)/np.abs(mx_c)*100, "percent")
      return_flag = True
    if (np.abs(mx2_c-mx2_prev)>(np.abs(mx2_c)/resolution)):
      print("WARNING: ETA CONCENTRATIONS give a different max value compared to the previous version.")
      print(np.abs(mx2_c-mx2_prev)/np.abs(mx2_c)*100, "percent")
      return_flag = True
    #f.write(tmp)
    #f.write("\n")

  #Wet deposition
  temp=np.abs(ndata['WD_spec001'][:,:,1,:,:])
  mn_c=np.mean(temp)
  mx_c=np.max(temp)
  temp2=np.abs(ndata2['WD_spec001'][:,:,1,:,:])
  mn2_c=np.mean(temp2)
  mx2_c=np.max(temp2)
  f.write(" WET deposition:\t\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn_c,mx_c))
  f.write(" WET deposition ETA:\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn2_c,mx2_c))
  if compare:
    tmp=fr.readline()
    tmp2=tmp.split()
    mn_prev=float(tmp2[4])
    mx_prev=float(tmp2[7])
    tmp=fr.readline()
    tmp2=tmp.split()
    mn2_prev=float(tmp2[5])
    mx2_prev=float(tmp2[8])
    if (np.abs(mn_c-mn_prev)>(np.abs(mn_c)/resolution)): 
      print("WARNING: Wet deposition gives a different mean value compared to the previous version.")
      print(np.abs(mn_c-mn_prev)/np.abs(mn_c)*100, "percent")
      return_flag = True
    if (np.abs(mn2_c-mn2_prev)>(np.abs(mn2_c)/resolution)): 
      print("WARNING: Wet ETA deposition gives a different mean value compared to the previous version.")
      print(np.abs(mn2_c-mn2_prev)/np.abs(mn2_c)*100, "percent")
      return_flag = True
    if (np.abs(mx_c-mx_prev)>(np.abs(mx_c)/resolution)):
      print("WARNING: Wet deposition gives  a different max value compared to the previous version.")
      print(np.abs(mx_c-mx_prev)/np.abs(mx_c)*100, "percent")
      return_flag = True
    if (np.abs(mx2_c-mx2_prev)>(np.abs(mx2_c)/resolution)):
      print("WARNING: Wet ETA deposition gives  a different max value compared to the previous version.")
      print(np.abs(mx2_c-mx2_prev)/np.abs(mx2_c)*100, "percent")
      return_flag = True
    #f.write(tmp)
    #f.write("\n")

#Dry deposition
  temp=np.abs(ndata['DD_spec001'][:,:,1,:,:])
  mn_c=np.mean(temp)
  mx_c=np.max(temp)
  temp2=np.abs(ndata2['DD_spec001'][:,:,1,:,:])
  mn2_c=np.mean(temp2)
  mx2_c=np.max(temp2)
  f.write(" DRY deposition:\t\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn_c,mx_c))
  f.write(" DRY deposition ETA:\tmean abs: %2.8e\tmax abs: %2.8e\n" %(mn2_c,mx2_c))
  if compare:
    tmp=fr.readline()
    tmp2=tmp.split()
    mn_prev=float(tmp2[4])
    mx_prev=float(tmp2[7])
    tmp=fr.readline()
    tmp2=tmp.split()
    mn2_prev=float(tmp2[5])
    mx2_prev=float(tmp2[8])
    if (np.abs(mn_c-mn_prev)>(np.abs(mn_c)/resolution)): 
      print("WARNING: Dry deposition gives a different mean value compared to the previous version.")
      print(np.abs(mn_c-mn_prev)/np.abs(mn_c)*100, "percent")
      return_flag = True
    if (np.abs(mn2_c-mn2_prev)>(np.abs(mn2_c)/resolution)): 
      print("WARNING: Dry ETA deposition gives a different mean value compared to the previous version.")
      print(np.abs(mn2_c-mn2_prev)/np.abs(mn2_c)*100, "percent")
      return_flag = True
    if (np.abs(mx_c-mx_prev)>(np.abs(mx_c)/resolution)):
      print("WARNING: Dry deposition gives a different max value compared to the previous version.")
      print(np.abs(mx_c-mx_prev)/np.abs(mx_c)*100, "percent")
      return_flag = True
    if (np.abs(mx2_c-mx2_prev)>(np.abs(mx2_c)/resolution)):
      print("WARNING: Dry ETA deposition gives a different max value compared to the previous version.")
      print(np.abs(mx2_c-mx2_prev)/np.abs(mx2_c)*100, "percent")
      return_flag = True
    #f.write(tmp)
    #f.write("\n")

  ndata.close()
  ndata2.close()

if (return_flag):
  raise ValueError('Old and New version deviate')
elif (compare):
  print("Output from previous version is not significantly different from this version.")
