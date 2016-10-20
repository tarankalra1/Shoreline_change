%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
%
% USAGE:
% Reads lon,lat,shoreline angle from a file 
% Reads grid and bathy data from US East Coast solution 
% Requires get_wavedata.m that gets wave data interpolated 
%
% DESCRIPTION:
%
% INPUT:  (all must have same dimensions)
%
% OUTPUT:
%   Q = Along shore sediment flux  [m^3 s]
%
% DISCLAIMER:
%  This software is provided "as is" without warranty of any kind.  
%
% REFERENCE:
%  Code to find the 40 m depth contour from US East Coast grid, smoothen bathymetry if
%  needed, then interpolate wave properties for the new bathymetric points 
%  and feed them into Ashton and Murray Calculation (JGR, 2006, doi:10.1029/2005JF000422) 
%=========================================================================
clear all ; close all ; clc ;

% Get the grid parameters from a solution of US East Coast
grd_file='/media/gadar/DATADRIVE1/shoreline_change/coawst_forecast_feb2014/coawst_us_20140214_01.nc';

mask_rho=ncread(grd_file,'mask_rho');
h_bathy=ncread(grd_file,'h').*(mask_rho) ; % Get the bathymetry  
lon_rho=ncread(grd_file,'lon_rho');
lat_rho=ncread(grd_file,'lat_rho');
[imax,jmax]=size(lon_rho); % Get the dimensions of grid

depth=40;  % chosen depth
% draw the chosen depth contour line 
figure(1)
pcolorjw(lon_rho,lat_rho,h_bathy.*mask_rho)
caxis([0 300])
colorbar
title('US East Coast with 40m bathy')
saveFigure(gcf,'US_East_40m.jpg')
% hold on 
% [c,h]=contour(lon_rho,lat_rho,h,[depth,depth],'r*-');
hold on 
xfire=-73.165069; yfire=40.651961; % FI location 
plot(xfire,yfire,'go')
%ylim([40.15 41.45])
%xlim([-74.5 -69.0])

%%
%Step 1 Fetch the x,y points where the contour is of certain depth
% There were some junk values in matrix c of (x,y) so limit it
% to be within a certain longitude of -80.0 and -65.0 
% This is performed outside in another code "get_shoreline_lonlat.m" 
shore_FI=importdata('shoreline_FI_strip.txt');
lon_shore_FI=shore_FI(:,1); lat_shore_FI=shore_FI(:,2); 

shoreline_angle=importdata('shoreline_FI/FI_shoreline_angle.1'); 
tot_ang=length(shoreline_angle); 

[c,h]=contour(lon_rho,lat_rho,h_bathy,[depth,depth] ,'b*-');
hold on 
s=getcontourlines(c);
xcontour_depth=s(1).x;  %Get the first vector of x and y 
ycontour_depth=s(1).y;

figure(2)
plot(xcontour_depth,ycontour_depth,'g+')

hold on 
% hold on
% plot(lon_rho,lat_rho,'ko') 
% hold on 
% plot(lon_rho',lat_rho','ko')
%xlim([-73.53 -71.0])
%ylim([39.2 40.0])

%% Step 2
%%Now get the Fire Island coastline data from Ilgar's Strip grid and
%%get the nearest bathymetric location from all these points
XI_FI_grid=[lon_shore_FI,lat_shore_FI]; 
X=[xcontour_depth',ycontour_depth'];
% Returns indices of nearest bathymetric location from X array 
K=dsearchn(X,XI_FI_grid); 
xnear_contour_depth=xcontour_depth(K); 
ynear_contour_depth=ycontour_depth(K);
hold on 
plot(lon_shore_FI,lat_shore_FI,'ro')
title('Extracted 40m bathy with FI in red')
saveFigure(gcf,'40m_US_East_FI.jpg')

figure(3)
plot(lon_shore_FI,lat_shore_FI,'ro')
hold on 
plot(xnear_contour_depth,ynear_contour_depth,'bx')
hold on 
plot(xcontour_depth,ycontour_depth,'go')
legend('shoreline-data','bathy-40m','nearbathy-40m')
xlim([-73.4 -72.7])
ylim([40.2 40.8])
title('Shoreline with nearest bathy')
saveFigure(gcf,'Shoreline_near_bathy_40m.jpg')
% A total of 1280 bathymetry points
% For each point of Ilgar strip grid, find closest 40 m contour point 
%x_strip,y_grid,'nearest')

%% Step 3
% Smoothen the bathymetry in the location that are at certain chosen
% depth above to Fire Island   
% % find the index of the bathymetric locations 
% index_matrix1 = 1:imax*jmax; 
% index_matrix1 = reshape(index_matrix1,size(lon_rho));
% % for i=1:5%length(xcontour_depth)
% lin_ind = griddata(lon_rho,lat_rho,index_matrix1,.....
%                    xnear_contour_depth,ynear_contour_depth,'nearest') ; %
% % Store the indices for all chosen bathymetric locations 
% [ipt,jpt]=ind2sub(size(lon_rho),lin_ind)   ;
%   
%% 
%Step 4
%Call an interpolation function to fit the values of wwave properties
%on the smoothened 40 meters contour line 
url='http://geoport.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd'; 
nc=ncgeodataset(url)  ; 
% Get time series --> jd contains the entire time series 
% and now find index of that time series 
% using nj_time to get time series of entire time
jd=nj_time(nc,('time')) ;
% corresponds to 01-Feb to 04-May 
% %itime=date_index(jd,[2014 02 01 0 0 0],[2014 05 04 0 0 0]) ;
itime=date_index(jd,[2014 02 01 0 0 0]);  
grd_size=imax*jmax; 
% reshape in column vectors
lon_rho_col=reshape(lon_rho,[grd_size 1]); 
lat_rho_col=reshape(lat_rho,[grd_size 1]);
ntotal=length(itime) 

for t=1:ntotal
  [F_Hwave,F_Lwave,F_DEwave,F_DNwave]=get_wavedata(itime(t),nc,grd_size,.....
                                                   lon_rho_col,lat_rho_col);
  for i=1:length(shoreline_angle)
% Number of points here is equal to the number of Fire Island stations
% Each station of FI needs a contour depth
    H=F_Hwave(xnear_contour_depth(i),ynear_contour_depth(i))  ; 
    L=F_Lwave(xnear_contour_depth(i),ynear_contour_depth(i))  ;
    DE=F_DEwave(xnear_contour_depth(i),ynear_contour_depth(i)); 
    DN=F_DNwave(xnear_contour_depth(i),ynear_contour_depth(i)); 
% Calculate the director vector from East And North components     
    D=atan2(DE,DN) ;
    if(D<0.0) 
% Add 2*pi as atan2 lies between -pi to pi 
      D=D+2.0*pi   ; 
    end 
% Convert radians to degrees for Ashton calculation     
    D=D*180.0/pi  ;
% Wave length to wave period 
    T=sqrt(L/1.56)  ;   
% Get Ashton calculation 
    [Q(i,t),angle_rel(i,t)]=Q_calc_Ashton(H,T,D,shoreline_angle(i));

  end 
end 
