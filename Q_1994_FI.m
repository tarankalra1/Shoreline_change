%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
%
% USAGE:
% Reads lon,lat,shoreline angle from a file that is separately generated
% Reads grid and bathy data from US East Coast solution 
% Requires get_wavedata.m that gets wave data interpolated 
% Requires read_swan_data.m - function that reads yearly SWAN model output
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
grd_file='/media/gadar/DATADRIVE2/CSI/USEast_1994_1999/USEast_grd31.nc';

mask_rho=ncread(grd_file,'mask_rho');
h_bathy=ncread(grd_file,'h').*(mask_rho) ; % Get the bathymetry  
lon_rho=ncread(grd_file,'lon_rho');
lat_rho=ncread(grd_file,'lat_rho');
[imax,jmax]=size(lon_rho); % Get the dimensions of grid

depth=40;  % chosen depth
% % draw the chosen depth contour line 
% figure(1)
% pcolorjw(lon_rho,lat_rho,h_bathy.*mask_rho)
% caxis([0 300])
% colorbar
% title('US East Coast with 40m bathy')
% saveFigure(gcf,'US_East_40m.jpg')
% % hold on 
% % [c,h]=contour(lon_rho,lat_rho,h,[depth,depth],'r*-');
% hold on 
% xfire=-73.165069; yfire=40.651961; % FI location 
% plot(xfire,yfire,'go')
% %ylim([40.15 41.45])
% %xlim([-74.5 -69.0])
% 
%%
% Step 1 Fetch the x,y points where the contour is of certain depth
% There were some junk values in matrix c of (x,y) so limit it
% to be within a certain longitude of -80.0 and -65.0 
% This is performed outside in another code "get_shoreline_lonlat.m" 
shore_FI=importdata('shoreline_FI/shoreline_FI_strip.txt');
lon_shore_FI=shore_FI(:,1); lat_shore_FI=shore_FI(:,2);
shoreline_angle=shore_FI(:,3); 
% 
% %shoreline_angle=importdata('shoreline_FI/FI_shoreline_angle.1'); 
% %tot_ang=length(shoreline_angle); 
% 
%[c,h]=contour(lon_rho,lat_rho,h_bathy,[depth,depth] ,'b*-');
[c,h]=contour(lon_rho,lat_rho,h_bathy,[depth,depth]);
% hold on 
s=getcontourlines(c);
xcontour_depth=s(1).x;  %Get the first vector of x and y 
ycontour_depth=s(1).y;
% 
%figure(2)
%plot(xcontour_depth,ycontour_depth,'g+')
%hold on   
%plot(lon_shore_FI,lat_shore_FI,'ro')
%title('Extracted 40m bathy with FI in red')
%saveFigure(gcf,'40m_US_East_FI.jpg')
% 
%% Step 2
% Now get the Fire Island coastline data from Ilgar's Strip grid and
% get the nearest bathymetric location from all these points
XI_FI_grid=[lon_shore_FI,lat_shore_FI]; 
X=[xcontour_depth',ycontour_depth'];
% Returns indices of nearest bathymetric location from X array 
K=dsearchn(X,XI_FI_grid); 
xnear_contour_depth=xcontour_depth(K); 
ynear_contour_depth=ycontour_depth(K);
% 
% figure(3)
% plot(lon_shore_FI,lat_shore_FI,'ro')
% hold on 
% plot(xnear_contour_depth,ynear_contour_depth,'bx')
% hold on 
% plot(xcontour_depth,ycontour_depth,'go')
% legend('shoreline-data','bathy-40m','nearbathy-40m')
% %xlim([-73.4 -72.7])
% %ylim([40.2 40.8])
% title('Shoreline with nearest bathy')
% saveFigure(gcf,'Output/Shoreline_near_bathy_40m.jpg')
% A total of 1280 bathymetry points
% For each point of Ilgar strip grid, find closest 40 m contour point 
%x_strip,y_grid,'nearest')

%% Step 3
% Smoothen the bathymetry in the location that are at certain chosen
% depth above to Fire Island    
%   
%% 
%Step 4
%Read the SWAN model output from .mat files
%Call an interpolation function to fit the values of wwave properties
% on the smoothened 40 meters contour line 

grd_size=imax*jmax; 
% reshape in column vectors
lon_rho_col=reshape(lon_rho,[grd_size 1]);
lat_rho_col=reshape(lat_rho,[grd_size 1]);
%  
year=1994; 
path_str='/media/gadar/DATADRIVE2/CSI/swan_1994_3hr';
% The months correspond to last time when data from SWAN output was stored
% For 1994, 04-01 months stored together, then individual months until 12
months={'12'};%,'05','06','07','08','09','10','11','12'};
ntotal=0; 

for i=1:length(months)

% Use these filenames before the math begins 
  year_month_str=strcat(num2str(year),'_',months(i));
  filename_wavetxt=char(strcat('Output/wave_char_',year_month_str,'.txt')); 
  filename_qtxt=char(strcat('Output/Q_',year_month_str,'.txt')); 
  filename_wavemat=char(strcat('Output/wave_char_',year_month_str','.mat')); 
  filename_qmat=char(strcat('Output/Q_',year_month_str,'.mat')); 
  
  Output_str=strcat('Output_',year_month_str)   ;
  wdir_str=strcat('wdir_',year_month_str,'.mat');
  hsig_str=strcat('hsig_',year_month_str,'.mat');
  wlen_str=strcat('wlen_',year_month_str,'.mat');
  
  wdir_files=strcat(path_str,'/',Output_str,'/',wdir_str); 
  hsig_files=strcat(path_str,'/',Output_str,'/',hsig_str);
  wlen_files=strcat(path_str,'/',Output_str,'/',wlen_str);
  
% Load all the data from each file 
  F_wdir=load(wdir_files{1});
  F_hsig=load(hsig_files{1}); 
  F_wlen=load(wlen_files{1}); 

% get the name of the wave characteristics 
  vname_wdir=fieldnames(F_wdir);
  vname_hsig=fieldnames(F_hsig);
  vname_wlen=fieldnames(F_wlen); 
%  
  numsteps=length(vname_wdir);
% Because last time instance at the end of each .mat is repeated
  numsteps=numsteps-1;
% numsteps for the 12th month should include all time snapshots
  if strcmp(months(i), '12')
    numsteps=numsteps+1;
  end
  
% Get the Julian time steps from .mat files
  for t=1:numsteps
% Based on the file data at all time steps, get the date 
    date=vname_wdir{t};  
    yr(t)=str2num(date(5:8)); 
    mo(t)=str2num(date(9:10));
    dy(t)=str2num(date(11:12));
    hr(t)=str2num(date(14:15));
    mm(t)=str2num(date(16:17));
    ss(t)=str2num(date(18:19));
    jul_times(t+ntotal)=datenum(yr(t),mo(t),dy(t),hr(t),mm(t),ss(t));
 
 % Fetch the wave charactersitics at each time step
    Wdir=(getfield(F_wdir,char(vname_wdir(t))))';
    Hsig=(getfield(F_hsig,char(vname_hsig(t))))';
    Wlen=(getfield(F_wlen,char(vname_wlen(t))))';   
    
    [F_H,F_L,F_DE,F_DN]=interp_wavedata_swan(lon_rho_col,.....
                                      lat_rho_col,grd_size,Hsig,Wlen,Wdir);
                                  
    for i=1:length(shoreline_angle)
% Number of points here is equal to the number of Fire Island stations
% Each station of FI needs a contour depth
      H(i,t+ntotal)=F_H(xnear_contour_depth(i),ynear_contour_depth(i))  ; 
      L(i,t+ntotal)=F_L(xnear_contour_depth(i),ynear_contour_depth(i))  ;
      DE=F_DE(xnear_contour_depth(i),ynear_contour_depth(i)); 
      DN=F_DN(xnear_contour_depth(i),ynear_contour_depth(i)); 
% Calculate the direction vector from East And North components     
      D=atan2(DE ,DN ) ;
      if(D<0.0) 
% Add 2*pi as atan2 lies between -pi to pi 
         D=D+2.0*pi   ; 
      end 
%% Convert radians to degrees for Ashton calculation     
      D1(i,t+ntotal)=D*180.0/pi  ;

%% Wave length to wave period 
      T=sqrt(L/1.56)  ;   
           
% Get Ashton calculation 
      [Q(i,t+ntotal),angle_rel(i,t+ntotal)]=Q_calc_Ashton(H(i,t+ntotal),....
            T(i,t+ntotal),D1(i,t+ntotal),shoreline_angle(i));

% Write the data in text format
      save_time(i,t+ntotal)=jul_times(t+ntotal); 
      
      C1=[H(i,t+ntotal),T(i,t+ntotal),D1(i,t+ntotal),save_time(i,t+ntotal)];                                       
      dlmwrite(filename_wavetxt,C1,'delimiter',',','precision',9,'-append') 
      C2=[Q(i,t+ntotal),angle_rel(i,t+ntotal),save_time(i,t+ntotal)];
      dlmwrite(filename_qtxt,C2,'delimiter',',','precision',9,'-append');

    end   
    clear F_H;F_L;F_DE;F_DN;Hsig;Wlen;Wdir;
  end
% Total steps  
  ntotal=ntotal+numsteps; 
end 
% Write the data in .mat
save(filename_wavemat,'H','T','D1','jul_times')  
save(filename_qmat,'Q','angle_rel','jul_times');
 
 clear Q; angle_rel; jul_times; H; T; D1; save_time;
