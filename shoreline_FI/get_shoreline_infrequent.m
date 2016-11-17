%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
%
% USAGE:  This code is used to get lon,lat (shoreline coordinates)
% from Fire Island's Strip Grid 
%
% DESCRIPTION:
%
% INPUT:
% 1. Strip grid netcdf file
% 2. getcontourlines.m -> get coordinates corresponding to a contour value
% 3. The code requires an angle calculation (azi) that is done using
% the mapping toolbox "azimuth" function which calls several other functions
%
% OUTPUT:
%   'shoreline_FI_strip.txt' containing -> lon_rho,lat_rho of FI strip grid
%   'lon_shore', 'lat_shore' -> shoreline data 
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
netcdf_load(grd_file);

mask_rho=ncread(grd_file,'mask_rho');
% Get the bathymetry  
h_bathy=ncread(grd_file,'h').*(mask_rho); 
lon_rho=ncread(grd_file,'lon_rho');
lat_rho=ncread(grd_file,'lat_rho');
 % Get the dimensions of grid
[imax,jmax]=size(lon_rho);

% chosen depth for coastline extraction
depth=0; 
figure(1)
pcolorjw(lon_rho,lat_rho,h_bathy)
caxis([0 30])
colorbar
hold on 
[c,h]=contour(lon_rho,lat_rho,h_bathy,[depth,depth] ,'b*-');
% This code gets you the coordinates corresponding to "depth"
% in a structure 
s=getcontourlines(c);
% Get the first vector of x and y 
xc_d=s(1).x;  
yc_d=s(1).y;

% Get the coordinates in the grid that correspond to the
% above found contour depth points 
index_matrix1 = 1:imax*jmax; 
index_matrix1 = reshape(index_matrix1,size(lon_rho));
lin_ind = griddata(lon_rho,lat_rho,index_matrix1,.....
                   xc_d,yc_d,'nearest') ; 
[ipt,jpt]=ind2sub(size(lon_rho),lin_ind)   ;

% Get the points in lon lat in a new array and a file  
for i=1:length(xc_d)
  lon_xc_d(i)=lon_rho(ipt(i),jpt(i));
  lat_yc_d(i)=lat_rho(ipt(i),jpt(i));
end 
% These files are containing several points that 
% are repeated and also some points that are parallel
% to longtidunal axis 
% Visualize 
figure(2)
% Only take last few points because you know that coastline
% is around
plot(lon_rho(:,(end-150:end)),lat_rho(:,(end-150:end)),'go')
hold on 
plot(lon_xc_d,lat_yc_d,'k+')

% Now remove the points that are exactly repeated 
A=unique([lon_xc_d' lat_yc_d'], 'rows', 'stable');
%dlmwrite('filename_unique.txt',A,'delimiter',',','precision',9) ; 

% Now remove the points that are paralleel 
% on a longitudnal axis so
% Find angle between two line segments using "azimuth" functions
% from the 
lon_pt=A(:,1); lat_pt=A(:,2); 
for i =1:length(A)-1
  azi(i)=azimuth(lat_pt(i),lon_pt(i),lat_pt(i+1),lon_pt(i+1));  
end
% because we didn't chose the last point in above azi calculation
% artificial manipulation of data 
azi(length(A))=20.0 ;  

%B=[lon_pt,lat_pt,azi'];
%dlmwrite('file_angle.txt',B,'delimiter',',','precision',9)

count=1;
% Now chose a threshold angle to omit the parallel points
thresh_ang=100.0; 
for i=1:length(A)
  if(azi(i)<thresh_ang)
   i
   lon_shore(count)=lon_pt(i);
   lat_shore(count)=lat_pt(i);
   azinew(count)=azi(i);
   count=count+1;
  end 
end
display('Total number of shoreline points')
count 
% Add the last point for the shoreline data because it was left out it
% in angle calculation. For now manually check if it alignes with the rest
% % of the points 
% lon_shore(count)=lon_pt(end); 
% lat_shore(count)=lat_pt(end);
hold on
%lon_shore(count+1)=lon_pt(1436);
%lat_shore(count+1)=lat_pt(1436); 
 plot(lon_shore,lat_shore,'rx')
legend('Last few FI pts.','Extracted pts.','Final Refined pts.')

%save('lonlat_FI_strip.mat' 'lon_shore' 'lat_shore')
 