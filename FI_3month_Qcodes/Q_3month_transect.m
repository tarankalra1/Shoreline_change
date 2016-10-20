echo off  
% Code to extract wave conditions from 01-Feb to 04-May COAWST forecast data 
% and feed them into Ashton and Murray (JGR, 2006, doi:10.1029/2005JF000422)
% Written by Tarandeep S. Kalra Oct 2016
  
clear all ;
close all ; 
clc ;
% Get the grid parameters from a single COAWST forecast file on local machine
% find the index of the depth that you use 
grd_file='/media/gadar/DATADRIVE1/shoreline_change/coawst_forecast_feb2014/coawst_us_20140214_01.nc';

% to get the nearest 300 meters point 
xfire=-73.133268; yfire=40.653496; % FI location 

% find the nearest offshore point in FI that is also 300 m deep 
mask_rho=ncread(grd_file,'mask_rho');
h=ncread(grd_file,'h').*(mask_rho) ; % Get the bathymetry  
lon_rho=ncread(grd_file,'lon_rho');
lat_rho=ncread(grd_file,'lat_rho');
[imax,jmax]=size(lon_rho); % Get the dimensions of grid

colors = 'ygbmkcr';
markers='+o.<>x^'; 
% Find the nearest point to Fire Island on grid 
dist=0.1; 
for j=1:jmax
  for i=1:imax
       % if (h(i,j)>=depth);
    dist_actual=sqrt((lon_rho(i,j)-xfire).^2+(lat_rho(i,j)-yfire).^2);
    if(dist_actual<=dist);
      dist=dist_actual ;
      istore=i ;
      jstore=j ;
    end
        %else
      % do nothing if h is not greater than depth(dp) m
   end      
      %end 
end

for dp=1:5
  jeff(dp)=jstore-(dp-1)*6;
  hh(dp)=h(istore,jeff(dp));
end 
% skip in space along latitude 30 kms first and then decrease the
% incremental distance to stay within 300 meters 
for dp=6:6
  jeff(dp)=jeff(5)-(dp-5)*2;
  hh(dp)=h(istore,jeff(dp));
end
for dp=7:7
  jeff(dp)=jeff(6)-(dp-6)*1;
  hh(dp)=h(istore,jeff(dp));
end  
% 
% % % Non straightened shoreline of Fire Island from Ilgar
shoreline_angle=importdata('shoreline_FI_Ilgar/FI_shoreline_angle.1'); 
tot_ang=length(shoreline_angle); 
alongshore=0:.04:51.16; % This is in kms., so a step for 40 meters
 

url='http://geoport.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd'; 
nc2=ncgeodataset(url)  ;  % corresponds to 01-Feb to 04-May 
jd=nj_time(nc2,'time') ;  % using nj_time to get time series of entire time
% Get time series --> jd contains the entire time series and now find index
% of that time series 
itime=date_index(jd,[2014 02 01 0 0 0],[2014 05 04 0 0 0]) ;
ntotal=length(itime);
date_time(:)=jd(itime(:));  % to store date_time 
%
for dp=1:length(jeff)
  for t=1:ntotal
    D=nc2{'Dwave'}(itime(t),jeff(dp),istore); 
    H=nc2{'Hwave'}(itime(t),jeff(dp),istore); 
    L=nc2{'Lwave'}(itime(t),jeff(dp),istore);
    T=sqrt(L/1.56); 
    for ialong=1:length(alongshore);      
      [Q(ialong,t,dp),angle_rel(ialong,t,dp)]=......
          Q_calc_Ashton(H,T,D,shoreline_angle(ialong));
    end
  end  
end 
save('Q_Ashton_vary_h.mat','Q','date_time','hh')