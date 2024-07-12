c Haversine function
      function distance(lona,lata,lonb,latb) result (dist)
# include "adinc.h"
# include "adparam.h"     
      double precision lata,lona,latb,lonb
      double precision dlat,dlon,rlata,rlatb
      double precision tmp1,tmp2,dist
      
      dlat = (latb-lata)*deg2rad
      dlon = (lonb-lona)*deg2rad
      rlata = lata*deg2rad
      rlatb = latb*deg2rad
      tmp1 = (sin(dlat/2))**2 + cos(rlata)*cos(rlatb)*(sin(dlon/2))**2
      tmp2 = 2*asin(sqrt(tmp1))
      dist = Eradius*tmp2
      end function distance
