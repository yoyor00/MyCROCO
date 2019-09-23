#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap, cm

points=np.loadtxt('cpoints.dat')
lonp = points[:, 0]
latp = points[:, 1]

nprocs=8

xlon=list(range(nprocs))
xlat=list(range(nprocs))
xz0b=list(range(nprocs))

lon = np.array([])
lat = np.array([])
dat = np.array([])
points = [None]*nprocs
for proc in range(nprocs):
#    data1=np.loadtxt('z0b.0{0}-011'.format(proc))
    data=np.loadtxt('z0b.{0:02d}-000'.format(proc))
    ni=int(data[0, 0])
    nj=int(data[0, 1])
    lon=np.append(lon, data[1:, 0]) #np.reshape(data[1:, 0], (nj, ni))
    lat=np.append(lat, data[1:, 1])
    dat=np.append(dat, (data[1:, 2]-0.008)/0.008)#np.reshape(data[1:, 1], (nj, ni))
    points[proc] = np.loadtxt('cpoints-{0:02d}.dat'.format(proc))
    print(len(points[proc]))
#    xlon[proc]=np.reshape(xlon[proc],(1,len(xlon[proc])))
#    xlat[proc]=np.reshape(xlat[proc],(1,len(xlat[proc])))
#    xz0b[proc]=np.reshape(xz0b[proc],(1,len(xz0b[proc])))

    
m = Basemap(projection='merc', llcrnrlon=-9, urcrnrlon=1,llcrnrlat=43,urcrnrlat=51, resolution='h')
 

#m.scatter(x, y, 3, marker='o', latlon=True)
#nx = int((m.xmax-m.xmin)/5000.)+1; ny = int((m.ymax-m.ymin)/5000.)+1
#    dat = m.transform_scalar(xz0b[proc], xlon[proc], xlat[proc], nx, ny)
#    im = m.imshow(dat,cm.GMT_haxby)
pc = m.pcolor(lon, lat, dat, latlon=True, tri=True)

for proc in range(nprocs):
    if (len(points[proc]) > 0):
        if (len(points[proc].shape)>1):
            x=points[proc][:, 0]
            y=points[proc][:, 1]
        else:
            x=points[proc][0]
            y=points[proc][1]
            
        m.plot(x, y, 'o', color=(1, float(proc)/float(nprocs), float(proc)/float(nprocs)), latlon=True)


#m.plot(lonp, latp, 'o', color=(0, 1, 0), latlon=True)        
# draw coastlines and borders
m.drawcoastlines()
m.drawcountries()
 
# draw meridians and parallels
m.drawmeridians(np.arange(0, 360, 2), labels=[False,True,True,False])
m.drawparallels(np.arange(-90, 90, 2), labels=[True,False,False,True])
cb = m.colorbar(pc,"right", size="5%", pad='2%')
plt.show()
