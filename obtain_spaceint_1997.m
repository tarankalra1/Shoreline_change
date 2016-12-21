clear all ; close all ; clc ;
%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
%
% USAGE:
% Save space integrated fluxes from .mat files 
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
i_13=326 ; 
% LOAD ILGAR's DATA
load('/media/gadar/DATADRIVE2/CSI/USEast_1994_1999/flux_FI_1994_1999_COAWST.mat');
for t=1:length(t97)
  Q_spaceint_Ilgar_1997(t)=0.0; 
%  for i=1:length(x)
  for i=i_13:i_13
    Q_spaceint_Ilgar_1997(t)=Q_spaceint_Ilgar_1997(t)+flux97(i,t); 
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
for t=1:length(time_Ashton_1997)
  Q_spaceint_Ashton_1997(t)=0.0 ;
  for i=i_13:i_13 % 1:length(lon_shore_FI)
    Q_spaceint_Ashton_1997(t)=Q_Ashton_1997(i,t)*factor+Q_spaceint_Ashton_1997(t); 
  end   
  Q_spaceint_Ashton_1997(t)=Q_spaceint_Ashton_1997(t); 
end 
save('/media/gadar/DATADRIVE2/CSI/Ashton/spaceint_fluxes/Q_97_spaceint',......
    't97','Q_spaceint_Ilgar_1997','time_Ashton_1997','Q_spaceint_Ashton_1997')
 
 
