%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a CROCO SOLITON_DJL configuration.
%  Create a grid file and an intial file for the vortex experiment
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
%  Ref:    Stastna M. and K.G. Lamb, (2002).
%          Large fully nonlinear internal solitary waves: 
%          The effect of background current. Physics of fluids.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SOLITON_DJL matlab package
%
% initialisation of the SOLITON_DJL croco config
% You have to set few parameter from the matlab file to adap to your domain size/depth. 
% H: the depth of you domain
% L:The lenght of your domain
% NX : number of x grid point
% NZ : number of vertical level (Z-grid, the interpolation to sigma coordinate is performed 
%      at the initilaisation of the croco model)
% 
% Note that an increase of resolution is asked later in the code, there you should provide the good NX, NZ 
%      for your croco grid 
% 
% The DJL algorithm ask an APE for the wave as initilisation. 
%
%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = 6*10e3;                            % APE for wave (m^4/s^2)
L  = 25000;                             % domain width (m)
H  = 600;                               % domain depth (m)
NX = 100;                               % grid
NZ = 50;                                % grid
%
f_shift = floor(NX*5/4);                  % shift to the left along X-direction
% 
% Equation of Seawater parameters, should be identical to your croco.in
% file
% lin_EOS_cff:  R0 [kg/m3], T0 [Celsius], S0 [PSU], TCOEF [1/Celsius], SCOEF [1/PSU] 
%              30.         11.5            0.        0.28                0.
R0=30.;
T0=11.5;
S0=0.;
TCOEFF=0.28;
SCOEFF=0.;
% 
%%%  
%
% The unitless density profile (normalized by a reference density rho0)
ref0=1028.0;
a_d=1.5/ref0; 
z0_d=100.0; 
d_d=30.0;
%
rho =@(z) 1-a_d*tanh((z+z0_d)/d_d);
%% go 
rhoz =@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;
%rho0 =1028;

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; 
Ubgz=@(z) 0*z; 
Ubgzz=@(z) 0*z;
%
%
croco_tools_path=''
