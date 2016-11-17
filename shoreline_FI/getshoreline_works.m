%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
% Guidance from Ilgar Safak
% USAGE:  This code is used to get lon,lat (shoreline coordinates)
% from Fire Island's Strip Grid 
%
% DESCRIPTION:
%
% INPUT:
% 1. Strip grid netcdf file
% 2. deg2utm.m -> convert degree to UTM 
% 3. get coordinates corresponding to a contour value
%
% OUTPUT:
%   'shoreline_FI_strip.txt' containing -> lon_rho,lat_rho,angle of FI strip grid
%   
%
% DISCLAIMER:
%  This software is provided "as is" without warranty of any kind.  
%
% REFERENCE:
%  
%  
clear all ; close all ; clc ; 

% Read strip grid netcdf file 
grd_file='/home/gadar/Documents/Development/Q_deep_FireIsland/strip_grid_FI/FI_strip_meanshore_extended.nc';

mask_rho=ncread(grd_file,'mask_rho');
% Get the bathymetry  
h_bathy=ncread(grd_file,'h').*(mask_rho); 
lon_rho=ncread(grd_file,'lon_rho');
lat_rho=ncread(grd_file,'lat_rho');
 % Get the dimensions of grid
[imax,jmax]=size(lon_rho);

% chosen depth for coastline extraction
depth=-0.46; 
figure(1)
pcolorjw(lon_rho,lat_rho,h_bathy)
caxis([0 30])
colorbar
hold on  

% ILGAR comments
% For each cross-grid line (ETA axis?) -- 
%and this should roughly be the gridlines along ~North-South orientation
%for Fire Island where we have 1280 cross-grid lines on strip grid) --
%I was finding the first point where the depth was greater than -0.46 
%and interpolating to the exact lat&lon values where the elevation is z=-0.46 (1280 lat&lon combo in total).
depth=-0.46; 
sten=2; % 2 point stencil for interpolation
for i=1:imax
    index=find(h_bathy(i,:)<depth,1);
     
    lon_rho_new(i)=interp1(h_bathy(i,index-sten:index+sten),......
                        lon_rho(i,index-sten:index+sten),depth);
    
    lat_rho_new(i)=interp1(h_bathy(i,index-sten:index+sten),......
                        lat_rho(i,index-sten:index+sten),depth);
end 
 
plot(lon_rho_new(:),lat_rho_new(:),'go')
hold on 

% Convert to DEG tO UTM 
[x,y,utmzone] = deg2utm(lat_rho_new,lon_rho_new);

% Get the angles
for i=1:imax-1
  ang(i)=(180/pi*atan((y(i+1)-y(i))/(x(i+1)-x(i))));  
end

% Fill the last angle
ang(imax)=ang(imax-1);

C=[lon_rho_new',lat_rho_new',ang'];
dlmwrite('shoreline_FI_strip.txt',C,'delimiter',',','precision',9)
% DIRECT EXTRACTION  
%[c,h]=contour(lon_rho,lat_rho,h_bathy,[depth,depth] ,'b*-');
%s=getcontourlines(c);
% Get the first vector of x and y 
%xc_d=s(1).x;  
%yc_d=s(1).y;
%plot(xc_d,yc_d,'r*')
