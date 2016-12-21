clear all ; close all ; clc ;
%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
%
% USAGE:
% Save time integrated fluxes from .mat files 
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
%=========================================================================
date_frmt='dd-mm-yyyy';
beginTime='01-01-1997'; finishTime='01-12-1997'; 
%% 
% LOAD ILGAR's DATA
load('/media/gadar/DATADRIVE2/CSI/USEast_1994_1999/flux_FI_1994_1999_COAWST.mat');
for i=1:length(x)
  Q_timeint_Ilgar_1997(i)=0.0; 
  for t=1:length(t97)
    Q_timeint_Ilgar_1997(i)=Q_timeint_Ilgar_1997(i)+flux97(i,t); 
  end 
end 
%% 
% Load Ashton calculations data 
shore_FI=importdata('/home/gadar/Documents/Development/Q_deep_FireIsland/shoreline_FI_strip.txt');
lon_shore_FI=shore_FI(:,1); lat_shore_FI=shore_FI(:,2);

[dist,phaseangle]=sw_dist(lat_shore_FI,lon_shore_FI,'km'); 

x_Ashton(1)=0.0 ; 
for i=2:length(lon_shore_FI)
    x_Ashton(i)=x_Ashton(i-1)+dist(i-1);  
end

load('/media/gadar/DATADRIVE2/CSI/Ashton/Q_Ashton_1997.mat');
factor=1; % Change this factor for matching units  
for i=1:length(lon_shore_FI)
  Q_timeint_Ashton_1997(i)=0.0 ;
  for t=1:length(time_Ashton_1997)
    Q_timeint_Ashton_1997(i)=Q_Ashton_1997(i,t)*factor+Q_timeint_Ashton_1997(i); 
  end   
  Q_timeint_Ashton_1997(i)=Q_timeint_Ashton_1997(i); 
end 
save('/media/gadar/DATADRIVE2/CSI/Ashton/timeint_fluxes/Q_97_timint',......
    'x','Q_timeint_Ilgar_1997','x_Ashton','Q_timeint_Ashton_1997')
 
 