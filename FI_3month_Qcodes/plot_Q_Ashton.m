clear all ; close all ; clc ;
load('COAWST_Forecast_FI_Q_Ashton_300mdeep.mat'); 
flux=Q ; 
  
% % % Non straightened shoreline of Fire Island from Ilgar
shoreline_angle=importdata('shoreline_FI_Ilgar/FI_shoreline_angle.1'); 
tot_ang=length(shoreline_angle); 
alongshore=0:.04:51.16; % This is in kms., so a step for 40 meters

t=datestr(date_time,6); % convert index to date string  

colormap(jet)
pcolorjw(alongshore(:)',date_time(:),flux) 
set(gca,'fontWeight','bold','Xtick',[0 5 10 15 20 25 30 35 40 45 50]);
xlabel('alonshore distance from W to E (km)')
datetick('y',6)
date_frmt='dd-mm-yyyy';
beginTime='01-02-2014'; finishTime='04-05-2014'; 
title('Alongshore sediment flux --negative is ~westward from Forecast Ashton form.') 
set(gca,'YLim',[datenum(beginTime,date_frmt),datenum(finishTime,date_frmt)]) 
caxis([-10 10])
colorbar
saveFigure(gcf,'Ashton formulation.jpg')

for ialong=1:length(alongshore)
    for it=1:length(t)
        if isnan(flux(ialong,it))
            flux(ialong,it)=0.0;
        end
    end
end

%flux_avg(:)=0.0 
for ialong=1:length(alongshore)
    flux_tavg(ialong)=0.0; 
    for it=1:length(t)
        %if(flux(ialong,it)~=NaN)
         flux_tavg(ialong)=flux(ialong,it)*3600+flux_tavg(ialong);
        %end 
    end
    flux_tavg(ialong)=flux_tavg(ialong)/length(t);
end 
figure(2)
plot(alongshore(:)',flux_tavg(:),'r*') 
set(gca,'fontWeight','bold','Xtick',[0 5 10 15 20 25 30 35 40 45]);
xlim([0 45])
xlabel('alonshore distance from W to E (km)')
ylabel('Time weighted sediment flux for 3 months Feb-May,2014')
title('Alongshore sediment flux --negative is ~westward, Ashton results') 
saveFigure(gcf,'sed_flux_Ashton.jpg')

% %% Ilgar data 
% load('/home/gadar/Documents/Development/Q_deep_FireIsland/LST_FI_2014/LST_FI_2014.mat');
% Q_Ilgar=flux;
% 
% for ialong=1:length(x)
%     flux_tavg_Ilgar(ialong)=0.0; 
%     for it=1:length(t)
%         %if(flux(ialong,it)~=NaN)
%          flux_tavg_Ilgar(ialong)=Q_Ilgar(ialong,it)+flux_tavg_Ilgar(ialong);
%         %end 
%     end
%     flux_tavg_Ilgar(ialong)=flux_tavg_Ilgar(ialong)/length(t);
% end 
% hold on 
% plot(x(:),flux_tavg_Ilgar(:)')
%     
%     
% % datestr(date_time(785)) --  08-Mar-2014 16:00:00