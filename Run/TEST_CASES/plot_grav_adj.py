#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from netCDF4 import Dataset

rootgrp = Dataset("gravadj_his.nc", "r")
Temp=rootgrp.variables["temp"][:]
Hz=rootgrp.variables["Hz"][:,:,2,:]
xrho=rootgrp.variables["x_rho"][2,:]
h=rootgrp.variables["h"][:]

xu=np.zeros(xrho.shape[0]-1)
for i in range(xu.shape[0]):
	xu[i]=0.5*(xrho[i]+xrho[i+1])

z_w_u=np.zeros((Hz.shape[1]+1,Hz.shape[2]-1))

x=np.zeros((Hz.shape[1]+1,xu.shape[0]))
for k in range(Hz.shape[1]+1):
	x[k,:]=xu

for i in range(xu.shape[0]):
	z_w_u[0,i]=-(0.5*(h[2,i]+h[2,i+1]))
	
fig=plt.figure()

for num in range(Temp.shape[3]):
	print 'frame = ',num
	for k in range(1,Hz.shape[1]+1):
		z_w_u[k,:]=z_w_u[k-1,:]+0.5*(Hz[num,k-1,:-1]+Hz[num,k-1,1:])
	plt.pcolor(x,z_w_u,Temp[num,:,2,1:-1],vmin=17.8571,vmax=35.7143)
	filename="grav_adj_{:0>4d}.png".format(num)
	plt.savefig(filename)
	plt.clf()
	
#for a movie use
# ffmpeg -i ./grav_adj_%4d.png XXX.mpg
