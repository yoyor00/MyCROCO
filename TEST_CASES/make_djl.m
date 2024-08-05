%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a CROCO initial file from SOLITON_DJL configuration
%
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
%  Copyright (c) 2024 by Jean-Baptiste ROUSTAN & Laurent ROBLOU 
%  e-mail:laurent.roblou@cnrs.fr
%
%  Version of 25-Jul-2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

addpath('./djles')

%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
%  Title 
%
CROCO_title  = 'SOLITON_DJL Run';
%
grdname  = 'croco_grd.nc';
ininame  = 'croco_ini.nc';
%
% Common parameters
%
crocotools_param_djl
djles_common
%
%
INTEGRATION_TEST=0;     % Set to 1 to activate debug tests  
%
CORRECT_UBAR=1;         % Set to 1 to correct u from vertical-averaged current
%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%

if ~exist('croco_tools_path','var') | (croco_tools_path == "")
    error('ERROR: do set your croco_tools_path in crocotools_param_djl.m file')
end

addpath([croco_tools_path,'/Preprocessing_tools'])
addpath([croco_tools_path,'/UTILITIES/mexcdf'])
addpath([croco_tools_path,'/UTILITIES/mexcdf/mexnc'])
addpath([croco_tools_path,'/UTILITIES/mexcdf/netcdf_toolbox/netcdf'])
addpath([croco_tools_path,'/UTILITIES/mexcdf/netcdf_toolbox/netcdf/ncsource'])
addpath([croco_tools_path,'/UTILITIES/mexcdf/netcdf_toolbox/netcdf/nctype'])
addpath([croco_tools_path,'/UTILITIES/mexcdf/netcdf_toolbox/netcdf/ncutility'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
disp([' Title: ',CROCO_title])
disp(' ')
%
%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;
% Find the solution of the DJL equation
djles_refine_solution

% Increase the resolution, and iterate to convergence
epsilon=1e-6;
NX=200; NZ=50;
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute (and plot the diagnostics)
%  pattern shif IS NOT applied at this stage
djles_diagnostics
djles_plot          


%%%% SAVE INPUT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp([' Making the grid: ',grdname])
disp(' ')
disp([' Resolution: 1/',num2str(L/NX),' m'])    

%
% Create the grid file
%
disp(' ')
disp(' Create the grid file...')
[Np,Lp]=size(w);
Mp=1;
L=Lp-1; M=Mp-1; N=Np-1;

disp([' LLm = ',num2str(L-1)])
%disp([' MMm = ',num2str(M-1)])
create_grid_djl(L,M,grdname,CROCO_title)

%
% Fill the grid file
%
disp(' ')
disp(' Fill the grid file...')
%
nc=netcdf(grdname,'write');
nc{'x_rho'}(:)=x;
nc{'y_rho'}(:)=x.*0;    
nc{'angle'}(:)=0.;       
nc{'spherical'}(:)='F';
nc{'mask_rho'}(:)=1+0.*x;
nc{'h'}(:)=x.*0+H;              
close(nc);

%
% reguraly-spaced vertical grid
%
% see https://www.myroms.org/wiki/Vertical_S-coordinate
%
vtransform=2;
theta_s=0;
theta_b=0;
hc=1.e16; 

if INTEGRATION_TEST == 1
    % test sigma coordinates stretching parameters :
    % compute ZC equiv. levels for rho grid,
    %    with Np requested rho levels
    % and compare to ZC (print tolerance)
    myzr=zlevs_1d(H, 0, theta_s, theta_b, hc, Np, 'r', vtransform)';
    tol=min(squeeze(ZC(:,1))-squeeze(myzr(1,:))');
    % -> stretching parameters are VALIDATED
end

% update ranks for vertical dimension
N=Np;
Np=Np+1;

%
% Create the initial condition file
%
disp(' ')
disp(' Create the ini file...')
create_inifile_djl(ininame,grdname,CROCO_title,...
                theta_s,theta_b,hc,N,...
                0.,'clobber',vtransform); % all physical variables preset to zero
%
% Fill the initial condition file
%
disp(' ')
disp(' Fill the ini file...')

nc=netcdf(ininame,'write');

% density and temperature variables, 
%   already on inner grid, equiv. to rho grid
%   to be shifted
cff=circshift(density*ref0, [0 -f_shift]);  
nc{'rho'}(1,:,:,:)=reshape(cff, [N Mp Lp]);
cff1=(cff-ref0)/TCOEFF + T0;     
nc{'temp'}(1,:,:,:)=reshape(cff1, [N Mp Lp]);

% horizontal velocity, 
%   ubar corrected (TBC)
%   already on inner grid
%   to be shifted
cff=circshift(u, [0 -f_shift]);
cff1=rho2u_2d(cff);
if CORRECT_UBAR == 1
   disp(' Reset u removing depth-averaged current...')
   myzw=zlevs_1d(H, 0, theta_s, theta_b, hc, N, 'w', vtransform)'; 
   cff=myzw(2:end)-myzw(1:end-1); % w level heights
   mydz=repmat(cff, Lp,1)';       % conform along X direction 
   dzu=rho2u_2d(mydz);            % and move to u points  
   hu=squeeze(sum(dzu.*cff1));     
   D_u=squeeze(sum(dzu));          
   ubar=squeeze(hu./D_u);
   %mystats=sprintf('min=%e ; avg=%e ; max=%e',min(ubar),mean(ubar),max(ubar));
   %disp(['correct u from ubar: ', mystats]);
   cff1=cff1-ubar; 
end
nc{'u'}(1,:,:,1:L)=reshape(cff1, [N Mp L]);

% vertical velocity, 
%   from inner grid to exterior grid
%   bottom/surface already set to zero
%   to be shifted
cff=circshift(w, [0 -f_shift]);   
for k=2:N-1
    nc{'w'}(1,k,:,:)=0.5*(cff(k,:)+cff(k+1,:));  
end
close(nc)

disp('Completed ')
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


