function add_arrow(lonmin,lonmax,latmin,latmax,cunit,cscale,width,height)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add an arrow in a vector plot (under m_map)
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x(1),y(1)]=m_ll2xy(lonmin,latmin); %coin bas gauche
[x(2),y(2)]=m_ll2xy(lonmax,latmin); %coin bas droite
[x(3),y(3)]=m_ll2xy(lonmin,latmax); %coin haut gauche
[x(4),y(4)]=m_ll2xy(lonmax,latmax); %coin haut droite
xmin=min(x);
xmax=max(x);
ymin=min(y);
ymax=max(y);
fac=0.01;
long=lonmin+fac*(lonmax-lonmin);
lat=latmax-fac*(latmax-latmin);
U=cunit*cscale;
V=0;
[X,Y]=m_ll2xy(long,lat,'clip','point');
[XN,YN]=m_ll2xy(long,lat+.01,'clip','point');
[XE,YE]=m_ll2xy(long+(.01)./cos(lat*pi/180),lat,'clip','point');
mU=U.*(XE-X)*100 + V.*(XN-X)*100;
mV=U.*(YE-Y)*100 + V.*(YN-Y)*100;
larrow=sqrt(mU.^2+mV.^2);
ratio1=(xmax-xmin)/(ymax-ymin);
ratio2=width/height;
if ratio2>=ratio1
 width=ratio1*width/ratio2;
else
 height=ratio2*width/ratio1;
end 
subplot('position',[0.6 0.06-height width height])
quiver(X,Y,larrow,0,0,'k')
text(X+larrow/2,Y,[num2str(cunit),' m.s^{-1}'],...
'HorizontalAlignment','center',...
'VerticalAlignment','bottom')
axis image                 
axis([xmin xmax ymin ymax])
set(gca,'visible','off')
