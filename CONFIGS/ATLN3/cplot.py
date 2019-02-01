#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

#points=np.loadtxt('cpoints.dat')
points=np.loadtxt('fort.778')

m = Basemap(projection='merc', llcrnrlon=-9, urcrnrlon=1,llcrnrlat=43,urcrnrlat=51, resolution='h')
 
# draw coastlines and borders
m.drawcoastlines()
m.drawcountries()
 
# draw meridians and parallels
m.drawmeridians(np.arange(0, 360, 30))
m.drawparallels(np.arange(-90, 90, 30))
for x,y in points:
    print (x,y)
    m.plot(x,y,'o',latlon=True)
plt.show()
