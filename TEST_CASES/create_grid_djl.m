function  create_grid(L,M,grdname,title)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf gridfile
%       L: total number of psi points in x direction  
%       M: total number of psi points in y direction  
%       grdname: name of the grid file
%       title: title in the netcdf file  
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
%  Copyright (c) 2024 by Laurent ROBLOU 
%  e-mail:laurent.roblou@cnrs.fr
%
%  Version of 25-Jul-2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lp=L+1;
Mp=M+1;

nw = netcdf(grdname, 'clobber');
%result = redef(nw);

%
%  Create dimensions
%

nw('xi_rho') = Lp;
nw('eta_rho') = Mp;
nw('one') = 1;

%
%  Create variables and attributes
%

nw{'xl'} = ncdouble('one');
nw{'xl'}.long_name = ncchar('domain length in the XI-direction');
nw{'xl'}.long_name = 'domain length in the XI-direction';
nw{'xl'}.units = ncchar('meter');
nw{'xl'}.units = 'meter';

nw{'el'} = ncdouble('one');
nw{'el'}.long_name = ncchar('domain length in the ETA-direction');
nw{'el'}.long_name = 'domain length in the ETA-direction';
nw{'el'}.units = ncchar('meter');
nw{'el'}.units = 'meter';
 
nw{'spherical'} = ncchar('one');
nw{'spherical'}.long_name = ncchar('Grid type logical switch');
nw{'spherical'}.long_name = 'Grid type logical switch';
nw{'spherical'}.option_T = ncchar('spherical');
nw{'spherical'}.option_T = 'spherical';

nw{'h'} = ncdouble('eta_rho', 'xi_rho');
nw{'h'}.long_name = ncchar('Final bathymetry at RHO-points');
nw{'h'}.long_name = 'Final bathymetry at RHO-points';
nw{'h'}.units = ncchar('meter');
nw{'h'}.units = 'meter';

nw{'f'} = ncdouble('eta_rho', 'xi_rho');
nw{'f'}.long_name = ncchar('Coriolis parameter at RHO-points');
nw{'f'}.long_name = 'Coriolis parameter at RHO-points';
nw{'f'}.units = ncchar('second-1');
nw{'f'}.units = 'second-1';

nw{'pm'} = ncdouble('eta_rho', 'xi_rho');
nw{'pm'}.long_name = ncchar('curvilinear coordinate metric in XI');
nw{'pm'}.long_name = 'curvilinear coordinate metric in XI';
nw{'pm'}.units = ncchar('meter-1');
nw{'pm'}.units = 'meter-1';

nw{'pn'} = ncdouble('eta_rho', 'xi_rho');
nw{'pn'}.long_name = ncchar('curvilinear coordinate metric in ETA');
nw{'pn'}.long_name = 'curvilinear coordinate metric in ETA';
nw{'pn'}.units = ncchar('meter-1');
nw{'pn'}.units = 'meter-1';

nw{'x_rho'} = ncdouble('eta_rho', 'xi_rho');
nw{'x_rho'}.long_name = ncchar('x location of RHO-points');
nw{'x_rho'}.long_name = 'x location of RHO-points';
nw{'x_rho'}.units = ncchar('meter');
nw{'x_rho'}.units = 'meter';

nw{'y_rho'} = ncdouble('eta_rho', 'xi_rho');
nw{'y_rho'}.long_name = ncchar('y location of RHO-points');
nw{'y_rho'}.long_name = 'y location of RHO-points';
nw{'y_rho'}.units = ncchar('meter');
nw{'y_rho'}.units = 'meter';

nw{'mask_rho'} = ncdouble('eta_rho', 'xi_rho');
nw{'mask_rho'}.long_name = ncchar('mask on RHO-points');
nw{'mask_rho'}.long_name = 'mask on RHO-points';
nw{'mask_rho'}.option_0 = ncchar('land');
nw{'mask_rho'}.option_0 = 'land';
nw{'mask_rho'}.option_1 = ncchar('water');
nw{'mask_rho'}.option_1 = 'water';

%result = endef(nw);

%
% Create global attributes
%

nw.title = ncchar(title);
nw.title = title;
nw.date = ncchar(date);
nw.date = date;
nw.type = ncchar('CROCO grid file');
nw.type = 'CROCO grid file';

close(nw);
