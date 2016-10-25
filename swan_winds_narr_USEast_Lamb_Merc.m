% Adjust year
clear all ; close all ; clc ; 

uwindfile='uwnd.10m.1994.nc'; 
vwindfile='vwnd.10m.1994.nc';

ncu = netcdf.open('uwnd.10m.1994.nc','NC_NOWRITE');
ncv = netcdf.open('vwnd.10m.1994.nc','NC_NOWRITE');

file_name='wind_USeast_1994_Taran.dat';%file_name='wind_USeast_2012_format.dat';

% Create a 1/4 degree lat-lon grid
lon_rho=-102:.25:-53;
lat_rho=11:.25:49;

varid = netcdf.inqVarID(ncu,'time'); 
NARR_time = netcdf.getVar(ncu,varid);
NARR_time = NARR_time./24+julian(1800,1,1,0);
NARR_time = datenum(gregorian(NARR_time));

narr_wind_start = 1;
narr_wind_end   = 2;  % length(NARR_time);
ntimes=narr_wind_end-narr_wind_start+1;

%varid = netcdf.inqVarID(ncu,'lon'); 
%nlon = netcdf.getVar(ncu,varid);

%varid = netcdf.inqVarID(ncu,'lat'); 
%nlat = netcdf.getVar(ncu,varid);

%numx=size(nlon,1);
%numy=size(nlon,2);
nlon = ncread(uwindfile,'lon');   
nlat = ncread(vwindfile,'lat');
numx=size(nlon,1);          
numy=size(nlon,2);
 
varname='uwnd'; 
%varid = netcdf.inqVarID(ncu,varname);
%uwnd = netcdf.getVar(ncu,varid,[0 0 narr_wind_start-1],[numx numy ntimes]);
%attname = netcdf.inqAttName(ncu,varid,11)
% uwnd_factor= netcdf.getAtt(ncu,varid,'scale_factor');
%uwnd_factor=1; 
uwind_lamb=ncread(uwindfile,varname,[1 1 narr_wind_start],[numx numy ntimes]);

varname='vwnd'; 
%varid = netcdf.inqVarID(ncv,varname);
vwind_lamb=ncread(vwindfile,varname,[1 1 narr_wind_start],[numx numy ntimes]);
%vwnd = netcdf.getVar(ncv,varid,[0 0 narr_wind_start-1],[numx numy ntimes]);
%vwnd_factor= netcdf.getAtt(ncv,varid,'scale_factor');
%vwnd_factor=1; 

lon_rho=repmat(lon_rho',1,length(lat_rho));
lat_rho=repmat(lat_rho,size(lon_rho,1),1);

[Lp,Mp]=size(lon_rho);
L=Lp-1;
M=Mp-1;
%u_merc=zeros(Lp,Mp,ntimes);
%v_merzeros(Lp,Mp,ntimes);

netcdf.close(ncu);
netcdf.close(ncv);

%   Rotate winds to earth lon lat based on http://ruc.noaa.gov/RUC.faq.html%
%   ROTCON_P          R  WIND ROTATION CONSTANT, = 1 FOR POLAR STEREO
%                         AND SIN(LAT_TAN_P) FOR LAMBERT CONFORMAL
%   LON_XX_P          R  MERIDIAN ALIGNED WITH CARTESIAN X-AXIS(DEG)
%   LAT_TAN_P         R  LATITUDE AT LAMBERT CONFORMAL PROJECTION
%                         IS TRUE (DEG)
 clat=50.0   ; 
 clon=-107.0 ; 
 lat_tan_p  =  clat;                     % 50.0 for NARR;
 lon_xx_p   =  clon;                     % -107.0 for NARR;
 rotcon_p   =  sin(lat_tan_p*pi/180);
 deg2rad=2*pi/360;
 angle2 = rotcon_p*(lon_rho-lon_xx_p)*deg2rad;
 sinx2 = sin(angle2);
 cosx2 = cos(angle2);
 
count = 0;
tic;
for i = 1:ntimes
    i
    count=count+1;
    disp(['Winds for USER grid at ',datestr(NARR_time(narr_wind_start-1+i))])
%    un = double(uwnd(:,:,i)).*uwnd_factor;
    un=double(uwind_lamb(:,:,i)); 
    
    A=un>999;
    un(A)=0;
    F  = TriScatteredInterp(double(nlon(:)),double(nlat(:)),double(un(:)));
    Uint_lamb(:,:,count) = F(double(lon_rho),double(lat_rho));
    
    clear F A
    
    vn = double(vwind_lamb(:,:,i));
    B=vn>999;
    vn(B)=0;
    F  = TriScatteredInterp(double(nlon(:)),double(nlat(:)),double(vn(:)));
    Vint_lamb(:,:,count) = F(double(lon_rho),double(lat_rho));
    
    clear F B
    
    % Convert the interpolated Lamber Conformal velocities to Mercator Projection for
    % User grids
      %

      u_merc(:,:,count)= cosx2.*Uint_lamb(:,:,count)+......
                          sinx2.*Vint_lamb(:,:,count);
      v_merc(:,:,count)=-sinx2.*Uint_lamb(:,:,count)+.....
                         cosx2.*Vint_lamb(:,:,count);
 
    % Store the Mercator projected velocity variables at each time
    

end; 
toc;
 
% Interpolate wind information from NARR onto the USER gri
fid = fopen(file_name,'w');
for i=1:count 
    uswan=squeeze(u_merc(:,:,i)'); 
    vswan=squeeze(v_merc(:,:,i)');
    fprintf(fid,'%10.2f\n',uswan');
    fprintf(fid,'%10.2f\n',vswan');
end 
fclose(fid);  

