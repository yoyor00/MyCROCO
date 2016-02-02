function interp_ECMWF(NCEP_dir,Y,M,Roa,interp_method,...
                     lon1,lat1,mask1,tin,...
		     nc_frc,nc_blk,lon,lat,angle,tout)
%
% Read the local ECMWF files and perform the space interpolations
%
%  Illig, 2010, from interp_NCEP
%  Updated    January 2016 (E. Cadier and S. Illig)
%---------------------------------------------------------------------------------
%
% 1: Air temperature: Convert from Kelvin to Celsius
%
vname='T2M';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
tair=squeeze(nc{vname}(tin,:,:));
close(nc);
tair=get_missing_val(lon1,lat1,mask1.*tair,nan,Roa,nan);
tair=tair-273.15;
tair=interp2(lon1,lat1,tair,lon,lat,interp_method);
%
%
% 2: Relative humidity: Convert from % to fraction
%
% Get Specific Humidity [Kg/Kg]
%
vname='Q';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
shum=squeeze(nc{vname}(tin,:,:));
close(nc);
shum=get_missing_val(lon1,lat1,mask1.*shum,nan,Roa,nan);
shum=interp2(lon1,lat1,shum,lon,lat,interp_method);

%
% computes specific humidity at saturation (Tetens  formula)
% (see air_sea tools, fonction qsat)
%
rhum=shum./qsat(tair);
%
% 3: Precipitation rate: Convert from [kg/m^2/s] to cm/day
%
vname='TP';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
prate=squeeze(nc{vname}(tin,:,:));
close(nc);
prate=get_missing_val(lon1,lat1,mask1.*prate,nan,Roa,nan);
prate=prate*0.1*(24.*60.*60.0);
prate=interp2(lon1,lat1,prate,lon,lat,interp_method);
prate(prate<1.e-4)=0;
%
% 4: Shortwave flux: [W/m^2]
%      ROMS convention: downward = positive
%
%  Solar shortwave
%
vname='SSR';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
dswrf=squeeze(nc{vname}(tin,:,:));
close(nc);
radsw=get_missing_val(lon1,lat1,mask1.*dswrf,nan,Roa,nan);
radsw=interp2(lon1,lat1,radsw,lon,lat,interp_method);
radsw(radsw<1.e-10)=0;
%
%
% 5: Longwave flux:  [W/m^2]
%      ROMS convention: positive upward.
%
%  5.1 Get the net longwave flux [W/m^2]
%
vname='STR';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
dlwrf=squeeze(nc{vname}(tin,:,:));
close(nc);
radlw=get_missing_val(lon1,lat1,mask1.*dlwrf,nan,Roa,nan);
radlw=interp2(lon1,lat1,radlw,lon,lat,interp_method);
%
%  5.2 get the downwrad longwave flux [W/m^2]
%
vname='STRD';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
dlwrf_in=squeeze(nc{vname}(tin,:,:));
close(nc);
radlw_in=get_missing_val(lon1,lat1,mask1.*dlwrf_in,nan,Roa,nan);
radlw_in=interp2(lon1,lat1,radlw_in,lon,lat,interp_method);
%
%
% 6: Wind  [m/s]
%
vname='U10M';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
uwnd=squeeze(nc{vname}(tin,:,:));
close(nc)
uwnd=get_missing_val(lon1,lat1,mask1.*uwnd,nan,Roa,nan);
uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
%
%
vname='V10M';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
vwnd=squeeze(nc{vname}(tin,:,:));
close(nc)
vwnd=get_missing_val(lon1,lat1,mask1.*vwnd,nan,Roa,nan);
vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);
%
wspd=sqrt(uwnd.^2+vwnd.^2);
%
%
% 7: Wind & Wind stress [m/s]
%
vname='EWSS';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
tx=squeeze(nc{vname}(tin,:,:));
close(nc)
tx=get_missing_val(lon1,lat1,mask1.*tx,nan,0,nan);
tx=interp2(lon1,lat1,tx,lon,lat,interp_method);
%
vname='NSSS';
nc=netcdf([NCEP_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc']);
ty=squeeze(nc{vname}(tin,:,:));
close(nc)
ty=get_missing_val(lon1,lat1,mask1.*ty,nan,0,nan);
ty=interp2(lon1,lat1,ty,lon,lat,interp_method);
%
% Compute the stress
%
%[Cd,uu]=cdnlp(wspd,10.);
%rhoa=air_dens(tair,rhum*100);
%tx2=Cd.*rhoa.*uwnd.*wspd;
%ty2=Cd.*rhoa.*vwnd.*wspd;

%
% Rotations on the ROMS grid
%
cosa=cos(angle);
sina=sin(angle);
%
sustr=rho2u_2d(tx.*cosa+ty.*sina);
svstr=rho2v_2d(ty.*cosa-tx.*sina);
%
% uwnd et vwnd sont aux points 'rho'
%
u10=rho2u_2d(uwnd.*cosa+vwnd.*sina);
v10=rho2v_2d(vwnd.*cosa-uwnd.*sina);
%
% Fill the ROMS files
%
if ~isempty(nc_frc)
  nc_frc{'sustr'}(tout,:,:)=sustr;
  nc_frc{'svstr'}(tout,:,:)=svstr;
end
if ~isempty(nc_blk)
  nc_blk{'tair'}(tout,:,:)=tair;
  nc_blk{'rhum'}(tout,:,:)=rhum;
  nc_blk{'prate'}(tout,:,:)=prate;
  nc_blk{'wspd'}(tout,:,:)=wspd;
  nc_blk{'radlw'}(tout,:,:)=radlw;
  nc_blk{'radlw_in'}(tout,:,:)=radlw_in;
  nc_blk{'radsw'}(tout,:,:)=radsw;
  nc_blk{'uwnd'}(tout,:,:)=u10;
  nc_blk{'vwnd'}(tout,:,:)=v10;
  %nc_blk{'sustr'}(tout,:,:)=sustr;
  %nc_blk{'svstr'}(tout,:,:)=svstr;
end


