clear all ; close all ; clc ;
%=========================================================================
% AUTHOR:  Tarandeep S. Kalra Oct 18-2014  (tkalra@usgs.gov)
%
% USAGE:
% this code will check the wave characteristics at different stations
% to know if they give different values 
% for 1994 .mat files 
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

load('/media/gadar/DATADRIVE2/CSI/Ashton/wave_char_Ashton_1997.mat');
 [Hunique, ind97]=unique(Hwave_Ashton_1997(:,1));
  % At different stations
  
  figure(1)
for i=1:length(ind97)
 %subplot(4,3,i)
 plot(time_Ashton_1997,Hwave_Ashton_1997(ind97(i),:));
 xlabel('time')
 ylabel('Hwave')
 title('1997 Hwave at 11 different stations')
 datetick('x',6)
 set(gca,'XLim',[datenum(beginTime,date_frmt),datenum(finishTime,date_frmt)])
 hold on 
end 
legend('1','2','4','5','6','7','8','9','10','11')

saveFigure(gcf,'Hwave_1997.jpg')
hold off 
figure(2)
for i=1:length(ind97)
% subplot(4,3,i)
 plot(time_Ashton_1997,Dwave_Ashton_1997(ind97(i),:));
 xlabel('time')
 ylabel('Dwave')
 title('1997 Dwave at 11 different stations')
 datetick('x',6)
 set(gca,'XLim',[datenum(beginTime,date_frmt),datenum(finishTime,date_frmt)])
 hold on 
end
legend('1','2','4','5','6','7','8','9','10','11')
saveFigure(gcf,'Dwave_1997.jpg')
% 
% 
% figure(3)
% for i=1:length(ind95)
%  subplot(4,3,i)
%  plot(time_Ashton_1995,Twave_Ashton_1995(ind95(i),:));
%  xlabel('time')
%  ylabel('Twave')
%  title('1995 Twave at different stations')
%  datetick('x',6)
%  set(gca,'XLim',[datenum(beginTime,date_frmt),datenum(finishTime,date_frmt)])
%  hold on 
% end 
% saveFigure(gcf,'Twave_1995.jpg')
% 
% %hold on   
% %end 
