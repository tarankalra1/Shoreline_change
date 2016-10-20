% Ilgar Safak, USGS, June 2016
% Modified by Tarandeep S. Kalra, Sept 2016

% This routine computes the alongshore sediment flux, Q (in m^3/s) 
% based on Ashton and Murray (JGR, 2006, doi:10.1029/2005JF000422), Eq.5
% for given:
% Inputs 
% (H in meters)  - significant wave height (H, in m)
% (D in degrees) - wave direction (D, in deg, where waves 
%                  come from, clockwise from N which is zero)
% (theta in degrees) - shoreline orientation (theta) in degrees
%    with respect to East, measured counterclockwise. 
%    For Fire Island, NY, this would increase from ~5 degrees (near Fire Island Inlet at west) to ~25 degrees (near Moriches Inlet at east)

% 
%  -
%  - wave period (T, in s)
%   all from deep water (could be buoy data, or COAWST or other model output)

%  A location is deep water for a wave with period T if h / Lo > 0.5 
%  h:  depth (m)
%  Lo: deep water wavelength (m) = 1.56*T^2    ;    1.56 comes from g/(2*pi) where g=9.806
%  An 8 s wave is in deep water at h > 50 m
%  A 14 s wave is in deep water at h > 306 m

% Notes:
% 1) The relative angle is being calculated for Fire Island
% 2) Negative means ~westward, positive means ~eastward transport for Fire Island

function [Q,angle_rel]=Q_calc_Ashton(H,T,D,theta)

if D<=(90-theta)
    angle_rel=0;
elseif D>=(270-theta)
    angle_rel=0;
elseif D>(90-theta) && D<=(180-theta)
    angle_rel=D-(180-theta);
elseif D>(180-theta) && D<(270-theta)
    angle_rel=D-(180-theta);
end        

% based on Ashton and Murray (JGR, 2006, doi:10.1029/2005JF000422), Eq.5
Q=.34*T^.2*H^2.4*sin(angle_rel*pi/180)*(cos(angle_rel*pi/180)^1.2);


