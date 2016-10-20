
function [F_H,F_L,F_DE,F_DN ] = get_wavedata(tidx,nc,grd_size,.....
                                             lon_rho_col,lat_rho_col)
%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
%
% USAGE:  This function is called at each time of time series 
%
% DESCRIPTION:
%
% INPUT: 
%   tidx= time index 
%   nc = handle for getting ncgeodataset wave variables 
%   lon_rho_col= column vector of lon
%   lat_rho_col= column vector of lat
% OUTPUT:
%   F_H  =  Interpolated function with wave height at chosen depth
%   F_L  =  Interpolated function with wave length at chosen depth 
%   F_DE =  Interpolated function with wave dir x at chosen depth
%   F_DN =  Interpolated function with wave dir y at chosen depth
% 
% DISCLAIMER:
%  This software is provided "as is" without warranty of any kind.  
%
% REFERENCE:
%  Code to interpolate wave properties for the new bathymetric points 
%  and feed them into Ashton and Murray Calculation (JGR, 2006, doi:10.1029/2005JF000422)

H=double(squeeze(nc{'Hwave'}(tidx,:,:))'); echo off ;
L=double(squeeze(nc{'Lwave'}(tidx,:,:))'); echo off ;
D=double(squeeze(nc{'Dwave'}(tidx,:,:))'); echo  off ; 

%Make NaN=0
H(isnan(H))=0; 
L(isnan(L))=0;
D(isnan(D))=0; 

%Before interpolation of the wave direction, split into East and North
%to avoid issues if (360 deg, 0deg), interpolated to 180 so split into
%two components and interpolate separately and then combine two vectors
%again
%Convert degrees to radians  
DEwave=1.0*sin(D*pi/180.0);
DNwave=1.0*cos(D*pi/180.0);

%Reshape all the data in column vector 
H=reshape(H,[grd_size 1]);
L=reshape(L,[grd_size 1]);
DEwave=reshape(DEwave,[grd_size 1]); 
DNwave=reshape(DNwave,[grd_size 1]); 
   
%Get the interpolant function 
F_H=scatteredInterpolant(lon_rho_col,lat_rho_col,H);
F_L=scatteredInterpolant(lon_rho_col,lat_rho_col,L);
F_DE=scatteredInterpolant(lon_rho_col,lat_rho_col,DEwave);
F_DN=scatteredInterpolant(lon_rho_col,lat_rho_col,DNwave);

 

