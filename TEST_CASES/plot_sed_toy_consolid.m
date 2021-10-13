%======================================================================
%
%     ---               Sed toy consolid Test Case                ---     
%
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
%  Patrick Marchesiello - 2012
%======================================================================
clear all
close all
%================== User defined parameters ===========================

fname='sed_toy_consolid_his.nc';

%======================================================================
%
% Read data
%

idy=1;
idz=1;
depth=20;

nc=netcdf(fname);

h=squeeze(nc{'h'}(idy,idz));
T=squeeze(nc{'scrum_time'}(:,:))/86400;

bostr=squeeze(nc{'bostr'}(:,idy,idz));
ALT=squeeze(nc{'act_thick'}(:,idy,idz));

s1=squeeze(nc{'sand_01'}(:,idy,idz));
s2=squeeze(nc{'sand_02'}(:,idy,idz));

m1=squeeze(nc{'mud_01'}(:,idy,idz));
m2=squeeze(nc{'mud_02'}(:,idy,idz));

sand1=s1*depth;
sand2=s2*depth;
mud1=m1*depth;
mud2=m2*depth;

close(nc);


%----------------------------------------------------------
%  Plot 
%----------------------------------------------------------
figure('Position',[1 1 800 800])

subplot(2,1,1)
% bostr
plot(T,bostr,'linewidth',1)
ylabel('Bottom stress (Pa)','fontsize',12);
title(['Double resuspension experiment'],'fontsize',14)
hold off
grid on

subplot(2,1,2)
% AL thickness
plot(T,ALT*100,'linewidth',1)
xlabel('Time (days)','fontsize',12);
ylabel('Active Layer thickness (cm)','fontsize',12);
hold off
grid on

set(gcf,'PaperPositionMode','auto');

export_fig -transparent sed_toy_consolid.pdf




return

