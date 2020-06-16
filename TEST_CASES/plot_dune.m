%======================================================================
%
%     ---               Dune Test Case                ---     
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

fname='dune_his.nc';

%======================================================================
%
% Read data
%
j=2;
tndx=17;
%
nc=netcdf(fname);
time=nc{'scrum_time'}(:)./86400;
tndx=min(tndx,length(time));
disp([ 'tndx = ',num2str(tndx), ...
    ' - Time = ',num2str(time(tndx)),' days' ])
h=squeeze(nc{'h'}(j,:));
x=squeeze(nc{'x_rho'}(j,:));
zeta=squeeze(nc{'zeta'}(tndx,j,:));
L=length(zeta);
%
hmorph  =squeeze(nc{'hmorph'}(tndx,j,:));
bedthick=squeeze(nc{'bed_thick'}(tndx,:,j,:));
bedfrac =squeeze(nc{'bed_frac_sand_2'}(tndx,:,j,:));
NL=size(bedfrac,1);
%
u=squeeze(nc{'u'}(tndx,:,j,:));
N=size(u,1);
%
theta_s=nc.theta_s(:);
theta_b=nc.theta_b(:);
hc=nc.hc(:);
Vtrans=nc{'Vtransform'}(:);
close(nc);
%
zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',Vtrans);
zr=squeeze(zr);
x =reshape(x,1,L);
xr=repmat(x,[N 1]);
xu=0.5*(xr(:,1:end-1)+xr(:,2:end));
zu=0.5*(zr(:,1:end-1)+zr(:,2:end));
%
NL=NL+1;
bedthick=padarray(bedthick,1,'replicate','pre');
bedfrac =padarray(bedfrac, 1,'replicate','pre');
bedfrac(2:end-1,:)=bedfrac(3:end,:);
xbed=repmat(x,[NL 1]);
h2=repmat(hmorph,[NL 1]);
zbed=-h2;
zbed(2:end,:)=-h2(2:end,:)-cumsum(bedthick(2:end,:),1);
%
%  Plot u and bed stratigraphy
%
figure('position',[500 500 600 600])
subplot(2,1,1)
hold on
pcolor(xbed,zbed,100*bedfrac);
line(xr,-hmorph,'color','k','Linewidth',2)
line(xr,-h,     'color','r','Linewidth',2,'linestyle','--')
hold off
caxis([30 70])
colorbar
axis([-Inf Inf -8 -1])
title(['DUNE Test Case - Fine Sand Fraction - Day ', ...
        num2str(time(tndx))])
set(gca,'fontsize',15);
set(gcf,'PaperPositionMode','auto');
hold off

nc=netcdf(fname);
hmorph=squeeze(nc{'hmorph'}(:,j,:));
close(nc)
subplot(2,1,2)
line(xr,-hmorph( 1,:),'color','k','Linewidth',2)
hold on
for n=1:5
 it=2^(n-1)+1;
 if it<=tndx
  line(xr,-hmorph(it,:),'color','k','Linewidth',2)
 end
end
hold off
grid on
axis([-1 101 -4.5 -1.5])
title(['DUNE Test Case - Bed evolution: 0 1 2 4 8 16 days'])
set(gca,'fontsize',15);
set(gcf,'PaperPositionMode','auto');

export_fig -transparent dune.pdf

return

