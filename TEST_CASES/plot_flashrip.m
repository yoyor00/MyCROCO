%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the FLASH_RIP test case
% 
%  Further Information:  
%  http://www.crocoagrif.org/
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
%  Ref: Penven, P., L. Debreu, P. Marchesiello and J.C. McWilliams,
%       Application of the ROMS embedding procedure for the Central 
%      California Upwelling System,  Ocean Modelling, 2006.
%
%  Patrick Marchesiello, IRD 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
hisname     = 'rip_his.nc';
avgname     = 'rip_avg.nc';

makepng   = 0;             % make png file
%
%======================================================================

% ---------------------------------------------------------------------
% --- get model grid and variables ---
% ---------------------------------------------------------------------

nc=netcdf(hisname);
tindex=length(nc{'scrum_time'}(:)); % reads last record
time=nc{'scrum_time'}(tindex)/60;

h=nc{'h'}(:);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:)-xl;
y=nc{'y_rho'}(:);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
N=length(nc('s_rho'));
Dcrit=nc{'Dcrit'}(:);
zeta=nc{'zeta'}(tindex,:,:);
close(nc)

nc=netcdf(avgname);
u=nc{'u'}(tindex-1,N,:,:);
v=nc{'v'}(tindex-1,N,:,:);
close(nc);

% Vorticity
vort=psi2rho(vorticity(u,v,pm,pn));

% Zeta
zeta(h<Dcrit)=zeta(h<Dcrit)-Dcrit;

% ---------------------------------------------------------------------
% --- plot ---
% ---------------------------------------------------------------------

hf = figure('position',[500 500 800 500]);
set(gca,'FontSize',15)

subplot(1,2,1)
cmin=-1; cmax=-cmin; cint=0.05;
nbcol=(cmax-cmin)/cint;
map=colormap(jet(nbcol));
colormap(map);
zeta=max(min(zeta,cmax),cmin);
zeta(zeta==0.)=NaN;
[C,hh]=contourf(x,y,zeta,[cmin:cint:cmax],'linestyle','none'); 
colorbar;
axis([-250 -10 0 300])
caxis([cmin cmax])
title(['Sea Level'])
set(gca,'FontSize',15)

subplot(1,2,2)
cmin=-0.05; cmax=-cmin; cint=0.005;
vort=max(min(vort,cmax),cmin);
vort(abs(vort)<cint)=NaN;
[C,hh]=contourf(x,y,vort,[cmin:cint:cmax],'linestyle','none'); 
colorbar;
axis([-250 -10 0 300])
caxis([cmin cmax])
title(['Wave-mean Vorticity'])
set(gca,'FontSize',15)

%----------------------------------

if makepng, 
    set(gcf,'PaperPositionMode','auto');
    export_fig -transparent flashrip.png
end

return





