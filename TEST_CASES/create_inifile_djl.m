function create_inifile(inifile,gridfile,title,...
                         theta_s,theta_b,hc,N,...
                         time,clobber, vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function nc=create_inifile(inifile,gridfile,theta_s,...
%                  theta_b,hc,N,ttime,stime,utime,... 
%                  cycle,clobber)
%
%   This function create the header of a Netcdf initial 
%   file.
%
%   Input: 
% 
%   inifile      Netcdf initial file name (character string).
%   gridfile     Netcdf grid file name (character string).
%   title        Netcdf initial file title (character string).
%   N            Number of w vertical levels.(Integer)  
%   time         Initial time.(Real) 
%   clobber      Switch to allow or not writing over an existing
%                file.(character string) 
%
%   Output
%
%   nc       Output netcdf object.
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
disp(' ')
disp([' Creating the file : ',inifile])
%
%  Read the grid file
%
nc=netcdf(gridfile,'r');
h=nc{'h'}(:);  
% mask=nc{'mask_rho'}(:);
close(nc);
hmin=min(min(h));
[Mp,Lp]=size(h);
L=Lp-1;
M=Mp-1;    
Np=N+1;
%
%  Create the initial file
%
type = 'INITIAL file' ; 
history = 'CROCO' ;
disp(['!rm ',inifile])
eval(['!rm ',inifile])
ncid=netcdf.create(inifile,'netcdf4');   

%
%  Create dimensions
%
xi_udimid = netcdf.defDim(ncid,'xi_u',L); 	
xi_vdimid = netcdf.defDim(ncid,'xi_v',Lp); 	
xi_rhodimid = netcdf.defDim(ncid,'xi_rho',Lp); 	
eta_udimid = netcdf.defDim(ncid,'eta_u',Mp); 	
eta_vdimid = netcdf.defDim(ncid,'eta_v',M); 	
eta_rhodimid = netcdf.defDim(ncid,'eta_rho',Mp); 
s_rhodimid = netcdf.defDim(ncid,'s_rho',N); 	
s_wdimid = netcdf.defDim(ncid,'s_w',Np); 	
tracerdimid = netcdf.defDim(ncid,'tracer',2); 	
timedimid = netcdf.defDim(ncid,'time',0); 
onedimid = netcdf.defDim(ncid,'one',1); 	
%
%  Create variables
%
netcdf.defVar(ncid,'spherical','char',onedimid); 
netcdf.defVar(ncid,'Vtransform','int',onedimid); 
netcdf.defVar(ncid,'Vstretching','int',onedimid); 
netcdf.defVar(ncid,'tstart','double',onedimid); 
netcdf.defVar(ncid,'tend','double',onedimid); 	
netcdf.defVar(ncid,'theta_s','double',onedimid); 
netcdf.defVar(ncid,'theta_b','double',onedimid); 
netcdf.defVar(ncid,'Tcline','double',onedimid); 
netcdf.defVar(ncid,'hc','double',onedimid); 	
netcdf.defVar(ncid,'s_rho','double',s_rhodimid); 
netcdf.defVar(ncid,'Cs_rho','double',s_rhodimid); 
netcdf.defVar(ncid,'ocean_time','double',timedimid); 
netcdf.defVar(ncid,'scrum_time','double',timedimid); 
varid = netcdf.defVar(ncid,'u','double',[xi_udimid,eta_udimid,s_rhodimid,timedimid]); 
varid = netcdf.defVar(ncid,'v','double',[xi_vdimid,eta_vdimid,s_rhodimid,timedimid]); 
varid = netcdf.defVar(ncid,'w','double',[xi_rhodimid,eta_rhodimid,s_wdimid,timedimid]); 
varid = netcdf.defVar(ncid,'ubar','double',[xi_udimid,eta_udimid,timedimid]); 	
varid = netcdf.defVar(ncid,'vbar','double',[xi_vdimid,eta_vdimid,timedimid]); 	
varid = netcdf.defVar(ncid,'zeta','double',[xi_rhodimid,eta_rhodimid,timedimid]); 	
varid = netcdf.defVar(ncid,'temp','double',[xi_rhodimid,eta_rhodimid,s_rhodimid,timedimid]); 
varid = netcdf.defVar(ncid,'salt','double',[xi_rhodimid,eta_rhodimid,s_rhodimid,timedimid]); 
varid = netcdf.defVar(ncid,'rho','double',[xi_rhodimid,eta_rhodimid,s_rhodimid,timedimid]); 
%
% Leave define mode
%
netcdf.close(ncid) 
%
%  Create attributes
%
ncwriteatt(inifile,'spherical','long_name','grid type logical switch');
ncwriteatt(inifile,'spherical','flag_values','T, F');
ncwriteatt(inifile,'spherical','flag_meanings','spherical Cartesian');
% 
ncwriteatt(inifile,'Vtransform','long_name','vertical terrain-following transformation equation');
% 
ncwriteatt(inifile,'Vstretching','long_name','vertical terrain-following stretching function');
%
ncwriteatt(inifile,'tstart','long_name','start processing day');
ncwriteatt(inifile,'tstart','units','day');
%
ncwriteatt(inifile,'tend','long_name','end processing day');
ncwriteatt(inifile,'tend','units','day');
%
ncwriteatt(inifile,'theta_s','long_name','S-coordinate surface control parameter');
ncwriteatt(inifile,'theta_s','units','nondimensional');
%
ncwriteatt(inifile,'theta_b','long_name','S-coordinate bottom control parameter');
ncwriteatt(inifile,'theta_b','units','nondimensional');
%
ncwriteatt(inifile,'Tcline','long_name','S-coordinate surface/bottom layer width');
ncwriteatt(inifile,'Tcline','units','meter');
%
ncwriteatt(inifile,'hc','long_name','S-coordinate parameter, critical depth');
ncwriteatt(inifile,'hc','units','meter');
%
ncwriteatt(inifile,'s_rho','long_name','S-coordinate at RHO-points');
ncwriteatt(inifile,'s_rho','units','nondimensional');
ncwriteatt(inifile,'s_rho','valid_min',-1);
ncwriteatt(inifile,'s_rho','valid_max',0);
ncwriteatt(inifile,'s_rho','positive','up');
if (vtransform ==1)
    ncwriteatt(inifile,'s_rho','standard_name','ocena_s_coordinate_g1');
elseif (vtransform ==2)
    ncwriteatt(inifile,'s_rho','standard_name','ocena_s_coordinate_g2');
end
ncwriteatt(inifile,'s_rho','formula_terms','s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
%
ncwriteatt(inifile,'Cs_rho','long_name','S-coordinate stretching curves at RHO-points');
ncwriteatt(inifile,'Cs_rho','units','nondimensional');
ncwriteatt(inifile,'Cs_rho','valid_min',-1);
ncwriteatt(inifile,'Cs_rho','valid_max',0);
%
ncwriteatt(inifile,'ocean_time','long_name','time since initialization');
ncwriteatt(inifile,'ocean_time','units','second');
%
ncwriteatt(inifile,'scrum_time','long_name','time since initialization');
ncwriteatt(inifile,'scrum_time','units','second');
%
ncwriteatt(inifile,'u','long_name','u-momentum component');
ncwriteatt(inifile,'u','units','meter second-1');
%
ncwriteatt(inifile,'v','long_name','v-momentum component');
ncwriteatt(inifile,'v','units','meter second-1');
%
ncwriteatt(inifile,'w','long_name','w-momentum component');
ncwriteatt(inifile,'w','units','meter second-1');
%
ncwriteatt(inifile,'ubar','long_name','vertically integrated u-momentum component');
ncwriteatt(inifile,'ubar','units','meter second-1');
%
ncwriteatt(inifile,'vbar','long_name','vertically integrated v-momentum component');
ncwriteatt(inifile,'vbar','units','meter second-1');
%
ncwriteatt(inifile,'zeta','long_name','free-surface');
ncwriteatt(inifile,'zeta','units','meter');
%
ncwriteatt(inifile,'temp','long_name','potential temperature');
ncwriteatt(inifile,'temp','units','Celsius');
%
ncwriteatt(inifile,'salt','long_name','salinity');
ncwriteatt(inifile,'salt','units','PSU');
%
ncwriteatt(inifile,'rho','long_name','density anomaly');
ncwriteatt(inifile,'rho','units','kilogram meter-3');
%
% Create global attributes
%
ncwriteatt(inifile,'/','title', title); 	
ncwriteatt(inifile,'/','date', date);
ncwriteatt(inifile,'/','clim_file', inifile);
ncwriteatt(inifile,'/','grd_file',gridfile);
ncwriteatt(inifile,'/','nc_type',type); 
ncwriteatt(inifile,'/','history',history); 
%
% Leave define mode
%
%%result = endef(nc);
%
% Compute S coordinates
%
[s_rho,Cs_rho,s_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%
% Write variables
%
ncwrite(inifile,'spherical','F');
ncwrite(inifile,'Vtransform',vtransform);   
ncwrite(inifile,'Vstretching',1);           
ncwrite(inifile,'tstart',time);              
ncwrite(inifile,'tend',time);                
ncwrite(inifile,'theta_s',theta_s);         
ncwrite(inifile,'theta_b',theta_b);         
ncwrite(inifile,'Tcline',hc);               
ncwrite(inifile,'hc',hc);                   
ncwrite(inifile,'s_rho',s_rho);               
ncwrite(inifile,'Cs_rho',Cs_rho);               
ncwrite(inifile,'ocean_time',time*24*3600);  
ncwrite(inifile,'scrum_time',time*24*3600);  
ncwrite(inifile,'u',zeros(L,Mp,N,1));           
ncwrite(inifile,'v',zeros(Lp,M,N,1));         
ncwrite(inifile,'w',zeros(Lp,Mp,Np,1));         
ncwrite(inifile,'zeta',zeros(Lp,Mp,1));         
ncwrite(inifile,'ubar',zeros(L,Mp,1));          
ncwrite(inifile,'vbar',zeros(Lp,M,1));         
ncwrite(inifile,'temp',zeros(Lp,Mp,N,1));       
ncwrite(inifile,'salt',zeros(Lp,Mp,N,1));       
ncwrite(inifile,'rho',zeros(Lp,Mp,N,1));       
return


