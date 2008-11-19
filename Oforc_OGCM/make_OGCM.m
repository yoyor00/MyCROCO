%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create and fill ROMS clim and bry files with OGCM data.
%
% The on-line reference to SODA is at
% http://iridl.ldeo.columbia.edu./SOURCES/.CARTON-GIESE/.SODA/
% The on-line reference to ECCO is at
% http://ecco.jpl.nasa.gov/cgi-bin/nph-dods/datasets/
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Contributions of P. Marchesiello (IRD), J. Lefevre (IRD),
%                   and F. Colberg (UCT)
%
%  Updated    6-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start % to be used in batch mode %
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param
%
if strcmp(OGCM,'SODA')
%
%  SODA DODS URL
%
% SODA 1.2:   Jan 1958 to Dec 2001 (now obsolete)
%  url = 'http://iridl.ldeo.columbia.edu./SOURCES/.CARTON-GIESE/.SODA/.v1p2'; 
% SODA 1.4.2: Jan 1958 to Dec 2001
%  url = 'http://iridl.ldeo.columbia.edu./SOURCES/.CARTON-GIESE/.SODA/.v1p4p2'; 
% SODA 1.4.3: Jan 2000 to Dec 2004 (QuickSCAT)
%  url =  'http://iridl.ldeo.columbia.edu./SOURCES/.CARTON-GIESE/.SODA/.v1p4p3'
% SODA 2.0.2-3.
%  url=    'http://iridl.ldeo.columbia.edu/SOURCES/.CARTON-GIESE/.SODA/.v2p0p2-3'
% SODA 2.0.2-4.
url= 'http://iridl.ldeo.columbia.edu/SOURCES/.CARTON-GIESE/.SODA/.v2p0p2-4'

disp({'!! Change from 1.4.3 version';
      '   New variables : taux, tauy, u, v';
      '   in m and m/s and no more in cm and cm/s'})

elseif strcmp(OGCM,'ECCO')
%
%  ECCO DODS URL
%
% Kalman filter 
%
%  url = 'http://ecco.jpl.nasa.gov/cgi-bin/nph-dods/datasets/kf049f/kf049f_'; 
  url = 'http://ecco.jpl.nasa.gov/cgi-bin/nph-dods/datasets/kf066b/kf066b_'; 
%
else
  error(['Unknown OGCM: ',OGCM])
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user input  parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if level==0
  nc_suffix='.nc';
else
  nc_suffix=['.nc.',num2str(level)];
  grdname=[grdname,'.',num2str(level)];
end
%
% Get the model grid
%
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
close(nc)
%
% Extract data over the internet
%
if Download_data==1
%
% Get the model limits
%
  lonmin=min(min(lon));
  lonmax=max(max(lon));
  latmin=min(min(lat));
  latmax=max(max(lat));
%
% Download data with DODS (the download matlab routine depends on the OGCM)
% 
  disp('Download data...')
  eval(['download_',OGCM,'(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,',...
                         'OGCM_dir,OGCM_prefix,url,Yorig)'])
%
end
%
% Get the OGCM grid 
% 
nc=netcdf([OGCM_dir,OGCM_prefix,'Y',num2str(Ymin),'M',num2str(Mmin),'.cdf']);
lonT=nc{'lonT'}(:);
latT=nc{'latT'}(:);
lonU=nc{'lonU'}(:);
latU=nc{'latU'}(:);
lonV=nc{'lonV'}(:);
latV=nc{'latV'}(:);
Z=-nc{'depth'}(:);
NZ=length(Z);
NZ=NZ-rmdepth;
Z=Z(1:NZ);
close(nc)
%
% Initial file 
% (the strategy is to start at the begining of a month)
% it is possible to do some temporal interpolation... 
% but I am too lazy. lets start the first day of
% month Mmin of year Ymin... with the first data available.
%
if makeini==1
  ininame=[ini_prefix,'Y',num2str(Ymin),'M',num2str(Mmin),nc_suffix];
%
% Process the time in Yorig time (i.e days since Yorig-01-01)
%
  tini=datenum(Ymin,Mmin,1)-datenum(Yorig,1,1);
  disp(['Create an initial file for ',datestr(tini+datenum(Yorig,1,1));])
  create_inifile(ininame,grdname,ROMS_title,...
                 theta_s,theta_b,hc,N,...
                 tini,'clobber');
  nc_ini=netcdf(ininame,'write');
  interp_OGCM(OGCM_dir,OGCM_prefix,Ymin,Mmin,Roa,interp_method,...
              lonU,latU,lonV,latV,lonT,latT,Z,1,...
              nc_ini,[],lon,lat,angle,h,1)
  close(nc_ini)
end
%
% Clim and Bry files 
%
if makeclim==1 | makebry==1
%
% Loop on the years and the months
%
  for Y=Ymin:Ymax
    if Y==Ymin 
      mo_min=Mmin;
    else
      mo_min=1;
    end
    if Y==Ymax
      mo_max=Mmax;
    else
    mo_max=12;
    end
    for M=mo_min:mo_max
      disp(' ')
      disp(['Processing  year ',num2str(Y),...
          ' - month ',num2str(M)])
      disp(' ')
%
      Mm=M-1;Ym=Y;
      if Mm==0
        Mm=12;
        Ym=Y-1;
      end
      Mp=M+1;Yp=Y;
      if Mp==13
        Mp=1;
        Yp=Y+1;
      end
%
% Add 2 times step in the ROMS files: 1 at the beginning and 1 at the end 
%
        nc=netcdf([OGCM_dir,OGCM_prefix,'Y',num2str(Y),'M',num2str(M),'.cdf']);
	OGCM_time=nc{'time'}(:);
	ntimes=length(OGCM_time);
	if ntimes==1
	  dt=30; % monthly files (SODA..)
	else
          dt=max(gradient(OGCM_time)); 
	end
	roms_time=0*(1:ntimes+2);
	roms_time(2:end-1)=OGCM_time;
	roms_time(1)=roms_time(2)-dt;
	roms_time(end)=roms_time(end-1)+dt;
	close(nc)
%
% Create and open the ROMS files
%
      if makebry==1
        bryname=[bry_prefix,'Y',num2str(Y),...
               'M',num2str(M),nc_suffix];
        create_bryfile(bryname,grdname,ROMS_title,[1 1 1 1],...
                       theta_s,theta_b,hc,N,...
                       roms_time,0,'clobber');
        nc_bry=netcdf(bryname,'write');
      else
        nc_bry=[];
      end
      if makeclim==1
        clmname=[clm_prefix,'Y',num2str(Y),...
                 'M',num2str(M),nc_suffix];
        create_climfile(clmname,grdname,ROMS_title,...
                        theta_s,theta_b,hc,N,...
                        roms_time,0,'clobber');
        nc_clm=netcdf(clmname,'write');
      else
        nc_clm=[];
      end
%
% Check if there are OGCM files for the previous Month
%
      fname=[OGCM_dir,OGCM_prefix,'Y',num2str(Ym),'M',num2str(Mm),'.cdf'];
      if exist(fname)==0
        disp(['   No data for the previous month: using current month'])
        Mm=M;
        Ym=Y;
	tndx_OGCM=1;
      else
        nc=netcdf(fname);
        tndx_OGCM=length(nc('time'));
        close(nc)
      end
%
% Perform the interpolations for the previous month
%
      disp(' Previous month :')
      interp_OGCM(OGCM_dir,OGCM_prefix,Ym,Mm,Roa,interp_method,...
                  lonU,latU,lonV,latV,lonT,latT,Z,tndx_OGCM,...
	  	  nc_clm,nc_bry,lon,lat,angle,h,1)
%
% Perform the interpolations for the current month
%
      for tndx_OGCM=1:ntimes
        disp([' Time step : ',num2str(tndx_OGCM),' of ',num2str(ntimes),' :'])
        interp_OGCM(OGCM_dir,OGCM_prefix,Y,M,Roa,interp_method,...
                    lonU,latU,lonV,latV,lonT,latT,Z,tndx_OGCM,...
		    nc_clm,nc_bry,lon,lat,angle,h,tndx_OGCM+1)
      end
%
% Read the OGCM file for the next month
%
      fname=[OGCM_dir,OGCM_prefix,'Y',num2str(Yp),'M',num2str(Mp),'.cdf'];
      if exist(fname)==0
         disp(['   No data for the next month: using current month'])
        Mp=M;
        Yp=Y;
	tndx_OGCM=ntimes;
      else
        nc=netcdf(fname);
	tndx_OGCM=1;
        close(nc)
      end
%
% Perform the interpolations for the next month
%
      disp(' Next month :')
      interp_OGCM(OGCM_dir,OGCM_prefix,Yp,Mp,Roa,interp_method,...
                  lonU,latU,lonV,latV,lonT,latT,Z,tndx_OGCM,...
		  nc_clm,nc_bry,lon,lat,angle,h,ntimes+2)
%
% Close the ROMS files
%
      if ~isempty(nc_clm)
        close(nc_clm);
      end
      if ~isempty(nc_bry)
        close(nc_bry);
      end
%
    end
  end
end
%
% Spin-up: (reproduce the first year 'SPIN_Long' times)
% just copy the files for the first year and change the time
%
if SPIN_Long>0
%
% Initial file
%
  if makeini==1
%
% Copy the file
%
    ininame=[ini_prefix,'Y',num2str(Ymin),'M',num2str(Mmin),nc_suffix];
    ininame2=[ini_prefix,'Y',num2str(Ymin-SPIN_Long),'M',num2str(Mmin),nc_suffix];
    disp(['Create ',ininame2]) 
    eval(['!cp ',ininame,' ',ininame2])
%
% Change the time
%
    nc=netcdf(ininame2,'write');
    time=nc{'scrum_time'}(:);
    time=time/(24*3600)+datenum(Yorig,1,1);
    [y,m,d,h,mi,s]=datevec(time);
    time=datenum(y-SPIN_Long,m,d,h,mi,s)-datenum(Yorig,1,1);
%    disp(datestr(time+datenum(Yorig,1,1)))
    nc{'scrum_time'}(:)=time*(24*3600);
    close(nc)
  end
%
  M=Mmin-1;
  Y=Ymin-SPIN_Long;
  for month=1:12*SPIN_Long
    M=M+1;
    if M==13
      M=1; 
      Y=Y+1;
    end
%
% Climatology files
%
    if makeclim==1
%
% Copy the file
%
      clmname=[clm_prefix,'Y',num2str(Ymin),'M',num2str(M),nc_suffix];
      clmname2=[clm_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
      disp(['Create ',clmname2]) 
      eval(['!cp ',clmname,' ',clmname2]) 
%
% Change the time
%
      nc=netcdf(clmname2,'write');
      time=nc{'tclm_time'}(:)+datenum(Yorig,1,1);
      [y,m,d,h,mi,s]=datevec(time);
      dy=Ymin-Y;
      y=y-dy;
      time=datenum(y,m,d,h,mi,s)-datenum(Yorig,1,1);
%      disp(datestr(time+datenum(Yorig,1,1)))
      nc{'tclm_time'}(:)=time;
      nc{'temp_time'}(:)=time;
      nc{'sclm_time'}(:)=time;
      nc{'salt_time'}(:)=time;
      nc{'uclm_time'}(:)=time;
      nc{'vclm_time'}(:)=time;
      nc{'v2d_time'}(:)=time;
      nc{'v3d_time'}(:)=time;
      nc{'ssh_time'}(:)=time;
      nc{'zeta_time'}(:)=time;
      close(nc)
    end
%
% Boundary files
%
    if makebry==1
%
% Copy the file
%
      bryname=[bry_prefix,'Y',num2str(Ymin),'M',num2str(M),nc_suffix];
      bryname2=[bry_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
      disp(['Create ',bryname2]) 
      eval(['!cp ',bryname,' ',bryname2]) 
%
% Change the time
%
      nc=netcdf(bryname2,'write');
      time=nc{'bry_time'}(:)+datenum(Yorig,1,1);
      [y,m,d,h,mi,s]=datevec(time);
      dy=Ymin-Y;
      y=y-dy;
      time=datenum(y,m,d,h,mi,s)-datenum(Yorig,1,1);
%      disp(datestr(time+datenum(Yorig,1,1)))
      nc{'bry_time'}(:)=time;
      close(nc)
    end
  end
end
%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
if makeplot==1
  disp(' ')
  disp(' Make a few plots...')
  if makeini==1
    ininame=[ini_prefix,'Y',num2str(Ymin),'M',num2str(Mmin),nc_suffix];
    figure
    test_clim(ininame,grdname,'temp',1,coastfileplot)
    figure
    test_clim(ininame,grdname,'salt',1,coastfileplot)
  end
  if makeclim==1
    clmname=[clm_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
    figure
    test_clim(clmname,grdname,'temp',1,coastfileplot)
    figure
    test_clim(clmname,grdname,'salt',1,coastfileplot)
  end
  if makebry==1
    bryname=[bry_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
    figure
    test_bry(bryname,grdname,'temp',1,obc)
    figure
    test_bry(bryname,grdname,'salt',1,obc)
    figure
    test_bry(bryname,grdname,'u',1,obc)
    figure
    test_bry(bryname,grdname,'v',1,obc)
  end
end


