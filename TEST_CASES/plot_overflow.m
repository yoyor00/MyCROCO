%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make 1 plot from the results of the Gravitational OVERFLOW test case
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
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
clear all
close all

makepdf=0;
tndx=6;
%
% Read data
%
nc=netcdf('over_his.nc');
h=nc{'h'}(:);
time=nc{'scrum_time'}(tndx)/86400;
y=squeeze(nc{'y_rho'}(:,2));
zeta=squeeze(nc{'zeta'}(tndx,:,:));
t0=squeeze(nc{'temp'}(1,:,:,2));
t=squeeze(nc{'temp'}(tndx,:,:,2));
[N,M]=size(t);
theta_s=nc.theta_s(:);
theta_b=nc.theta_b(:);
hc=nc.hc(:);
close(nc);

zr = zlevs(h,zeta,theta_s,theta_b,hc,N,'r');
zr=squeeze(zr(:,:,1));
yr=reshape(y,1,M);
yr=repmat(yr,[N 1])/1000;

subplot(2,1,1)
contourf(yr,zr,t0,(-0.1:0.1:1),'linestyle','none')
caxis([0 1])
colorbar
title('OVERFLOW - \rho anomaly vertical section')
xlabel('Y [km]')
ylabel('Z [m]')
text(10,-30,'t=0','fontsize',14)
set(gca,'fontsize',15)

subplot(2,1,2)
contourf(yr,zr,t,(-0.1:0.1:1),'linestyle','none')
caxis([0 1])
colorbar
xlabel('Y [km]')
ylabel('Z [m]')
text(10,-30,['t=',num2str(time,'%4.1f'),' days'],'fontsize',14)
set(gca,'fontsize',15)

if makepdf
 export_fig -transparent -pdf overflow.pdf
end



