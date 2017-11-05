%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SHOREFACE test case
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
%  ------------------------------------
%  Patrick Marchesiello, IRD 2013,2017
%
%  Ref: Uchiyama, Y., McWilliams, J., Shchepetkin, A., 2010. Wave-current 
%  interaction in an oceanic circulation model with a vortex-force formalism: 
%  application to the surf zone. Ocean Modell. 34, 16–35.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = 'shoreface_his.nc';  % croco file name
makepdf   = 0;                   % make pdf file

alphar    = 0.0;                 % roller fraction  (see croco.in)
Cd        = 0.0015;              % drag coefficient (see croco.in

wname     = 'TEST_CASES/shoreface_JW_frc.nc';
%
%======================================================================
%
g = 9.8;
yindex = 2;      
%
% plots axis
xmin0=-700;
xmax0=-30;
zmin0=-8;
zmax0=1;
%
% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname,'r');
tindex=length(nc{'scrum_time'}(:)); % reads last record
%
% horizontal grid
hr=squeeze(nc{'h'}(yindex,:));
xindex=1;
hr=hr(xindex:end);
L=length(hr);
xr=squeeze(nc{'x_rho'}(yindex,xindex:end));
yr=squeeze(nc{'y_rho'}(yindex,xindex:end));
dx=xr(2)-xr(1);
%
% vertical grid
N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(tindex,yindex,xindex:end));
zr=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'r',2));
dzr=zr(2:end,:)-zr(1:end-1,:);               % ---> zw(2:N,:)
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
dzru=zru(2:end,:)-zru(1:end-1,:);            % ---> zwu(2:N,:)
zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
dzw=zw(2:end,:)-zw(1:end-1,:);               % ---> zr
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
dzwu=zwu(2:end,:)-zwu(1:end-1,:);            % ---> zru
%
xr2d=repmat(xr,[N 1]);
xw2d=repmat(xr,[N+1 1]);
D=hr+zeta;
D2d=repmat(D,[N 1]);
h2d=repmat(hr,[N 1]);

% ---------------------------------------------------------------------
% --- read/compute numerical model fields ---
% --------------------------------------------------------------------
time=nc{'scrum_time'}(tindex)/86400;

Dcrit=nc{'Dcrit'}(:);

% ... zonal barotropic velocity ...              ---> xu,zru
ubar=squeeze(nc{'ubar'}(tindex,yindex,xindex:end));

% ... meridional barotropic velocity ...         ---> xr,zr
vbar=squeeze(nc{'vbar'}(tindex,yindex,xindex:end));

% ... zonal velocity ...                         ---> xu,zru
u=squeeze(nc{'u'}(tindex,:,yindex,xindex:end));

% ... meridional velocity ...                    ---> xr,zr
v=squeeze(nc{'v'}(tindex,:,yindex,xindex:end));

% ... vertical velocity ...                      ---> xr,zw
w=squeeze(nc{'w'}(tindex,:,yindex,xindex:end));

% ... vertical viscosity/diffusivity
Akv=squeeze(nc{'AKv'}(tindex,:,yindex,xindex:end));
Akt=squeeze(nc{'AKt'}(tindex,:,yindex,xindex:end));

% -------------- wave fields ------------------

if ~isempty(nc{'hrm'}(:)), 
 wkb_wwave=1; 
else
 wkb_wwave=0;
end

if wkb_wwave,
 % ... wave height ...  
 hrm=squeeze(nc{'hrm'}(tindex,yindex,xindex:end));

 % ... wave breaking dissipation ...  
 eb=squeeze(nc{'epb'}(tindex,yindex,xindex:end));

 % ... wave roller dissipation ...  
 er=squeeze(nc{'epr'}(tindex,yindex,xindex:end));

 % ... wave number x ...
 kx=squeeze(nc{'wkx'}(tindex,yindex,xindex:end));

 % ... wave number y ...
 ky=squeeze(nc{'wke'}(tindex,yindex,xindex:end));

 % ... wave frequency ...
 sig=squeeze(nc{'frq'}(tindex,yindex,xindex:end));

else

 % ... wave height ...  
 hrm=squeeze(nc{'whrm'}(tindex,yindex,xindex:end));

 % ... wave breaking dissipation ...  
 eb=squeeze(nc{'wdsp'}(tindex,yindex,xindex:end));

 % ... wave roller dissipation ...  
 er=0.;

 % ... wave number x ...
 wdrx=squeeze(nc{'wdrx'}(tindex,yindex,xindex:end));

 % ... wave number y ...
 wdry=squeeze(nc{'wdre'}(tindex,yindex,xindex:end));

 % ... wave frequency ...
 sig=squeeze(nc{'wfrq'}(tindex,yindex,xindex:end));
end

% ... wave setup ...  
sup=squeeze(nc{'sup'}(tindex,yindex,xindex:end));

% ... zonal Stokes dritf                         ---> xu,zru
ust2d=squeeze(nc{'ust2d'}(tindex,yindex,xindex:end));

% ... meridional Stokes drift ...                ---> xr,zr
vst2d=squeeze(nc{'vst2d'}(tindex,yindex,xindex:end));

% ... zonal Stokes dritf                         ---> xu,zru
ust=squeeze(nc{'ust'}(tindex,:,yindex,xindex:end));

% ... meridional Stokes drift ...                ---> xr,zr
vst=squeeze(nc{'vst'}(tindex,:,yindex,xindex:end));

% ... vertical Stokes drift ...                  ---> xr,zw
wst=squeeze(nc{'wst'}(tindex,:,yindex,xindex:end));

% eddy viscosity due to depth-induced wave breaking
Akb=squeeze(nc{'Akb'}(tindex,:,yindex,xindex:end));

% eddy diffusivity due to primary waves
Akw=squeeze(nc{'Akw'}(tindex,:,yindex,xindex:end));

close(nc)

%
% Wave height from Haas and Warner 2009 
% (SWAN propagating offshore JONSWAP wave spectrum)
%
nc=netcdf(wname,'r');
hrm_swan=2.*squeeze(nc{'Awave'}(1,yindex,xindex:end));
close(nc)

% ---------------------------------------------------------------------
% --- compute 2D analytical solutions ---
% ---------------------------------------------------------------------

if wkb_wwave==0        % Compute wkx,wke
 dd = D;
 cfrq=sig;
 khd = dd.*cfrq.*cfrq/g;
 kh = sqrt( khd.*khd + khd./(1.0 + khd.*(0.6666666666 ...
           +khd.*(0.3555555555   + khd.*(0.1608465608 ... 
           +khd.*(0.0632098765   + khd.*(0.0217540484 ...
           +khd.*0.0065407983)))))) );
 kw=kh./dd;
 cosw=wdrx;
 sinw=wdry;
 kx=kw.*cosw;
 ky=kw.*sinw;
end

eb=(1-alphar)*eb+er;
zeta_a=cumsum(dx*eb.*kx./(g*D.*sig))+zeta(1);
rhs=eb.*ky./(Cd*sig);
ub=u2rho_2d(ust2d);
vbar_a=-sqrt( 0.5*(sqrt(ub.^4+4*rhs.^2)-ub.^2) );
%vbar_a=-sqrt(abs(rhs)); %if ubar<<vbar
ubar=u2rho_2d(ubar);
ust2d=u2rho_2d(ust2d);
%============================================================
% --- plot ---
%=============================================================
%
xr  =xr  -1000;
xr2d=xr2d-1000;
xw2d=xw2d-1000;
ust(:,2:end)=0.5*(ust(:,1:end-1)+ust(:,2:end));
ust(:,1)=ust(:,2); ust(:,L)=ust(:,L-1);
u(:,2:end)=0.5*(u(:,1:end-1)+u(:,2:end));
u(:,1)=ust(:,2); u(:,L)=u(:,L-1);
zeta(D<Dcrit)  =NaN;
zeta_a(D<Dcrit)=NaN;
u(D2d<Dcrit)   =NaN;
v(D2d<Dcrit)   =NaN;
zeta(hr<0)     =zeta(hr<0)+hr(hr<0);
zeta_a(hr<0)   =zeta_a(hr<0)+hr(hr<0);
Akv(:,1)=NaN;
Akt(:,1)=NaN;
Akb(:,1)=NaN;
Akw(:,1)=NaN;
%
%--------------------------------------------------------
% Wave height, surface elevation and breaking dissipation
%--------------------------------------------------------
% 
figure('Position',[100 100 500 700])
subplot(5,1,1)
xmin=-1000; xmax=-50; zmin=0.; zmax=2.;
plot(xr,hrm,'color','b','LineWidth',3); hold on
plot(xr,hrm_swan,'color','r','LineWidth',2,'linestyle','--');
legend('CROCO','SWAN')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
title(['SHOREFACE: Wave Height at ',num2str(thour),' h'])
hold off
%
subplot(5,1,2)
zmin=-0.1; zmax=0.3;
plot(xr,zeta,'color','b','LineWidth',3); hold on;
plot(xr,zeta_a,'color','r','LineWidth',2,'linestyle','--');
legend('CROCO 3D','Analytical 2D','location','northwest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
title(['SHOREFACE: Wave Setup at ',num2str(thour),' h'])
hold off

subplot(5,1,3)
zmin=-0.15; zmax=0.1;
plot(xr,ubar,'color','b','LineWidth',3); hold on;
plot(xr,-ust2d,'color','r','LineWidth',2,'linestyle','--');
legend('CROCO 3D','Analytical 2D','location','northwest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
title(['SHOREFACE: Crosshore current [m/s] at ',num2str(thour),' h'])

subplot(5,1,4)
zmin=-1.2; zmax=0.5;
plot(xr,vbar,'color','b','LineWidth',3); hold on;
plot(xr,vbar_a,'color','r','LineWidth',2,'linestyle','--');
legend('CROCO 3D','Analytical 2D','location','northwest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
title(['SHOREFACE: Longshore current [m/s] at ',num2str(thour),' h'])

subplot(5,1,5)
zmin=0; zmax=0.1;
plot(xr,eb,'color','b','LineWidth',3);
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
title(['SHOREFACE: Wave breaking diss. [m^3/s^3] at ',num2str(thour),' h'])

if makepdf
 print -dpdf shoreface_z.pdf
 eval('!pdfcrop shoreface_z.pdf shoreface_z.pdf')
end
%
%-----------------------------------
% Eulerian velocities u,v,w
%-----------------------------------
%
figure('Position',[100 150 600 600])

subplot(3,1,1)
cmin=-0.5; cmax=0.5; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,u,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: U at ',num2str(thour),' h'])
hold off

subplot(3,1,2)
cmin=-1.2; cmax=1.2; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,v,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: V at ',num2str(thour),' h'])
hold off

subplot(3,1,3)
cmin=-6.e-3; cmax=6.e-3; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,w,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: W at ',num2str(thour),' h'])
hold off

if makepdf
 print -dpdf shoreface_u.pdf
 eval('!pdfcrop shoreface_u.pdf shoreface_u.pdf')
end

%
%-----------------------------------
% Stokes drift ust,vst,wst
%-----------------------------------
%
figure('Position',[100 150 600 600])

subplot(3,1,1)
cmin=-0.2; cmax=0.2; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,ust,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: UST at ',num2str(thour),' h'])
hold off

subplot(3,1,2)
cmin=-2.e-2; cmax=2.e-2; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,vst,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: VST at ',num2str(thour),' h'])
hold off

subplot(3,1,3)
cmin=-2.e-3; cmax=2.e-3; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,wst,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: WST at ',num2str(thour),' h'])
hold off

if makepdf
 print -dpdf shoreface_ust.pdf
 eval('!pdfcrop shoreface_ust.pdf shoreface_ust.pdf')
end
%
%-------------------------------------------------
% Turbulent and wave-induced Viscosity/diffusivity
%-------------------------------------------------
%
figure('Position',[100 150 600 600])

subplot(3,1,1)
cmin=0; cmax=0.05; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
colormap(map);
contourf(xw2d,zw,Akv,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: total viscosity AKV at ',num2str(thour),' h'])
hold off

subplot(3,1,2)
cmin=0; cmax=0.05; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
colormap(map);
contourf(xw2d,zw,Akb,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: wave-breaking viscosity AKB at ',num2str(thour),' h'])
hold off

subplot(3,1,3)
cmin=0; cmax=0.005; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
colormap(map);
contourf(xw2d,zw,Akw,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([xmin0 xmax0 zmin0 zmax0])
caxis([cmin cmax])
thour=floor(time*24);
title(['SHOREFACE: non-breaking wave diffusivity AKW at ',num2str(thour),' h'])
hold off

if makepdf
 print -dpdf shoreface_Ak.pdf
 eval('!pdfcrop shoreface_Ak.pdf shoreface_Ak.pdf')
end



