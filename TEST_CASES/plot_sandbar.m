%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SANDBAR test case
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
%  Patrick Marchesiello, IRD 2017,2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname      = 'sandbar_his.nc'; % croco file name

Expname    = '1B';     % LIP experiment : 1B: offshore migration
                       %                  1C: onshore  migration
if Expname == '1B',
 morph_fac = 18;       % morphological factor (from sediment.in)
else
 morph_fac = 13; 
end

morph_cpl  = 1;        % feedback to currents
makepdf    = 0;        % make pdf file
%
%======================================================================
%
% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------
%
yindex = 2; % Mm=1 with NS no-slip conditions
%
nc=netcdf(fname,'r');
tindex  =length(nc{'scrum_time'}(:)); % reads last record
tindex0 =min(tindex,5);               %  2 h
%
% horizontal grid
hr=squeeze(nc{'h'}(yindex,:));
xr=squeeze(nc{'x_rho'}(yindex,:));
hu=0.5*(hr(1:end-1)+hr(2:end));
xu=0.5*(xr(1:end-1)+xr(2:end));
L=length(hr);
%
% new bathy from last record
if morph_cpl,
 hnew=squeeze(nc{'hmorph'}(tindex,yindex,:));
 h=hnew;
else
 h=hr;
end
%
% vertical grid
N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(tindex,yindex,:));
Dcrit=1.1*nc{'Dcrit'}(:);
zeta(h<Dcrit)=zeta(h<Dcrit)-h(h<Dcrit); %add land topo
zr=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'r',2));
zw=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'w',2));
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
dz1=zr(1,:)-zw(1,:);
%
xr2d=repmat(xr,[N 1]);
xu2d=repmat(xu,[N 1]);
xw2d=repmat(xr,[N+1 1]);
D   =zw(N+1,:)-zw(1,:);
D2d =repmat(D,[N 1]);
Du  =zwu(N+1,:)-zwu(1,:);
Du2d=repmat(Du,[N 1]);

% ---------------------------------------------------------------------
% --- read/compute numerical model fields (index 1) ---
% --------------------------------------------------------------------
time=morph_fac/86400*(nc{'scrum_time'}(tindex)- ...
                      nc{'scrum_time'}(1));

% ... zonal velocity ...                         ---> xu,zu
u=squeeze(nc{'u'}(tindex,:,yindex,:));

% ... vertical velocity ...                      ---> xr,zw
w=squeeze(nc{'w'}(tindex,:,yindex,:));

% ... sediment concentration ...                 ---> xr,zr
Cbot=squeeze(nc{'sand_1'}(tindex0,1,yindex,:));
Cbot=Cbot.*(dz1./0.01).^1.2; % fitting to Rouse profile
Cbot=max(0.,Cbot);

% ... total viscosity
Akv=squeeze(nc{'AKv'}(tindex,:,yindex,:));

% ... viscosity due to wave breaking ...
Akb=squeeze(nc{'Akb'}(tindex,:,yindex,:));

% ... wave setup ...  
sup=squeeze(nc{'zeta'}(tindex0,yindex,:)); % init time
sup(hr<0)=sup(hr<0)+hr(hr<0)-Dcrit;

% ... u undertow ...
ubot =squeeze(nc{'u'}(tindex0,3,yindex,:)); % init time

% ... hrms ...  
hrms =squeeze(nc{'hrm'}(tindex0,yindex,:)); % init time

close(nc)

zeta(D<Dcrit)=NaN;
u(Du2d<Dcrit)=NaN;
sup(D<=max(0.1,Dcrit))=NaN;
sup(hr<0)=NaN;
%
%======================================================================
%     FLUME DATA (LIP: ROELVINK & RENIER 1995)
%======================================================================
%
xd=linspace(0,200,64);

% --- LIP-1B bathy at 18h on xd grid ---
hd1B = -[ ...
-4.0935  -4.0931  -4.0926  -4.0921  -4.0917  -4.0912  -4.0908  -3.9686 ...
-3.8283  -3.6879  -3.5487  -3.4097  -3.2665  -3.1167  -2.9668  -2.8105 ...
-2.6449  -2.5004  -2.4013  -2.3178  -2.2750  -2.2422  -2.2418  -2.2274 ...
-2.1987  -2.1699  -2.1234  -2.0759  -2.0055  -1.9297  -1.8476  -1.7596 ...
-1.7340  -1.7175  -1.6960  -1.6517  -1.6024  -1.5488  -1.4808  -1.4050 ...
-1.2944  -1.1205  -0.8597  -0.8404  -0.9656  -0.9895  -1.0150  -1.0028 ...
-0.9712  -0.8307  -0.5613  -0.4577  -0.4202  -0.4011  -0.3702  -0.3169 ...
-0.2464  -0.1274  -0.0041   0.1719   0.3584   0.6120   0.7816   0.8273];

% --- LIP-1C bathy at 13h on xd grid ---
hd1C = -[ ...
-4.1117  -4.1068  -4.1019  -4.0970  -4.0921  -4.0872  -4.0772  -3.9639 ...
-3.8228  -3.6655  -3.5212  -3.3955  -3.2532  -3.0899  -2.9443  -2.8001 ...
-2.6459  -2.4907  -2.3907  -2.3071  -2.2576  -2.2382  -2.2436  -2.1997 ...
-2.1558  -2.1240  -2.0958  -2.0408  -1.9742  -1.9114  -1.8500  -1.7977 ...
-1.7617  -1.7219  -1.6751  -1.6315  -1.5914  -1.5515  -1.5093  -1.4539 ...
-1.3661  -1.2169  -1.0205  -0.8062  -0.7222  -0.9540  -1.0534  -1.0534 ...
-1.0044  -0.8524  -0.5634  -0.5000  -0.4464  -0.3899  -0.3659  -0.3586 ...
-0.2218  -0.2098  -0.1132   0.1382   0.3922   0.6086   0.7579   0.8107];

% ----------

x_hrms_d1B=[20 65 100 115 130 138 145 152 160 170]; % --- LIP-1B 8hr
hrms_d1B  =[.85 .8 0.7 0.65 0.6 0.5 0.38 0.35 0.35 0.25];

x_ubot_d1B= [65 102 130 138 145 152 160 170];
  ubot_d1B=-[13  15  18  30  32  18  18  13];

x_Cbot_d1B= [65 102 130 138 145 152 160 170];
  Cbot_d1B= [0.3 0.2 0.9  3 1.9 0.8 0.9 0.9];

% -----------

x_hrms_d1C=[20  40  65 100 115 130 132 138 145 152 160 170]; % --- LIP-1C 7hr
hrms_d1C  =[.4 .41 .43 .44 .43 .43 .46 .43 .35 .33 .32 .22];

x_ubot_d1C= [65 102 115 125 130 134 152 160];
  ubot_d1C=-[ 1  1   1   2   2   3  13  11];

x_Cbot_d1C= [65  102 115 125 130  134 152 160];
  Cbot_d1C= [0.1 0.1 0.3 0.2 0.35 0.5 0.3 0.8];

if Expname=='1B',
 hd=hd1B;
 x_hrms_d=x_hrms_d1B;
 hrms_d=hrms_d1B;
 x_ubot_d=x_ubot_d1B;
 ubot_d=ubot_d1B;
 x_Cbot_d=x_Cbot_d1B;
 Cbot_d=Cbot_d1B;
else
 hm=interp1(xr,hr,xd);
 hd=hm+hd1C-hd1B;
 %hd=hd1C;
 x_hrms_d=x_hrms_d1C;
 hrms_d=hrms_d1C;
 x_ubot_d=x_ubot_d1C;
 ubot_d=ubot_d1C;
 x_Cbot_d=x_Cbot_d1C;
 Cbot_d=Cbot_d1C;
end
%============================================================
% --- plot ---
%=============================================================
%
figure('Position',[50 50 600 800])
xmin=80; xmax=190;
%
% section u ...
%
h1=subplot(4,1,1);
h1.Position = h1.Position + [0 -0.12 0 0.15];
%zmin=-Inf; zmax=Inf;
zmin=-2.5; zmax=0.5;
cmin=-0.7; cmax=0.7; nbcol=14;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xu2d,zru,u,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar('h');
leg(1)=plot(xr,-hr,  'k',  'LineWidth',2);
leg(2)=plot(xr,-hnew,'k--','LineWidth',3);
leg(3)=plot(xd,-hd  ,'r--','LineWidth',3);
plot(xr,zeta,'color','g','LineWidth',3);
legend(leg(1:3),'Initial','Final Model','Final data', ...
                'location','southeast');
ylabel('Depth [m]','Fontsize',15)
grid on
axis([xmin xmax zmin zmax])
caxis([cmin cmax])
thour=floor(time*24);
if Expname=='1B',
 title(['SANDBAR EROSION   LIP-1B - U at Time ',num2str(thour),' hour'])
else
 title(['SANDBAR ACCRETION LIP-1C - U at Time ',num2str(thour),' hour'])
end
set(gca,'Fontsize',15)
hold off
%
% hrms ...
%
h2=subplot(4,1,2);
h2.Position = h2.Position + [0 -0.05 0 -0.05];
%zmin=-Inf; zmax=Inf;
zmin=0; zmax=1;
plot(xr,hrms,'b',x_hrms_d,hrms_d,'b*','LineWidth',2); hold on
legend('Model','Flume','Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
ylabel('Hrms [m]','Fontsize',15)
set(gca,'Fontsize',15)
hold off
%
% Undertow ...
%
h3=subplot(4,1,3);
h3.Position = h3.Position + [0 -0.025 0 -0.05];
%zmin=-Inf; zmax=Inf;
zmin=-40; zmax=0;
plot(xu,100*ubot,'b',x_ubot_d,ubot_d,'b*','LineWidth',2); hold on;
legend('Model','Flume','Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
ylabel('U [cm/s]','Fontsize',15)
set(gca,'Fontsize',15)
hold off
%
% Sand Concentration ...
%
h4=subplot(4,1,4);
h4.Position = h4.Position + [0 0.0 0 -0.05];
%zmin=-Inf; zmax=Inf;
zmin=0; zmax=4;
plot(xr,Cbot,'b',x_Cbot_d,Cbot_d,'b*','LineWidth',2); hold on;
legend('Model','Flume','Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
xlabel('X [m]','Fontsize',15)
ylabel('C [g/l]','Fontsize',15)
set(gca,'Fontsize',15)
hold off

if makepdf
 if Expname=='1B',
  export_fig -transparent sandbar_LIP_1B.pdf
 else 
  export_fig -transparent sandbar_LIP_1C.pdf
 end
end

return
%

