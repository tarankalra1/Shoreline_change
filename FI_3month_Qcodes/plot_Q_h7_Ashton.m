clear all ; close all ; clc ;
load('Q_Ashton_vary_h.mat'); 
flux=Q ; 

[imax tmax hmax]=size(flux); 
hh=[7.3, 40.5, 53.76, 71.13, 124.46, 178.64, 292.57]; 

% % % Non straightened shoreline of Fire Island from Ilgar
shoreline_angle=importdata('shoreline_FI_Ilgar/FI_shoreline_angle.1'); 
tot_ang=length(shoreline_angle); 
alongshore=0:.04:51.16; % This is in kms., so a step for 40 meters

t=datestr(date_time,6); % convert index to date string  
%
for dp=1:hmax
  for ialong=1:imax
    for it=1:tmax
      if isnan(flux(ialong,it,dp))
        flux(ialong,it,dp)=0.0; 
      end
    end 
  end
end 

colors = 'ygbmkcr'; markers='+o.<>x^'; 
for dp=1:hmax
  for ialong=1:imax
    flux_tavg(ialong,dp)=0.0; 
    for it=1:tmax
      flux_tavg(ialong,dp)=flux(ialong,it,dp)*3600+flux_tavg(ialong,dp);
    end 
    flux_tavg(ialong,dp)=flux_tavg(ialong,dp)/length(t);
  end 
    h_line=plot(alongshore(1:imax),flux_tavg(1:imax,dp));
    set(h_line,'Color', colors(dp),'Marker',markers(dp),'MarkerSize',4) ;
    hold on
end
xlabel('Alongshore distance from Fire Island Inlet to Moriches Inlet (km)');
ylabel(['Alongshore sediment flux (Q,m^3/s)(positive: ~Eastward transport ;negative: "Westward transport)'])
set(gca,'fontWeight','bold','Xtick',[0 5 10 15 20 25 30 35 40 45]);
xlim([0 45])
title('Alongshore Sediment Flux, Ashton Equation, using USEast waves, varying depths, Feb-May 2014')
h_legend=legend(num2str(hh(1)),num2str(hh(2)),num2str(hh(3)),num2str(hh(4)).....
               ,num2str(hh(5)),num2str(hh(6)),num2str(hh(7)));
set(h_legend,'FontSize',20)
set(gca,'fontweight','bold')
saveFigure(gcf,'Varying_transect_integrated_3month_Q.jpg')


% figure(2)
% plot(alongshore(:)',flux_tavg(:),'r*') 
% set(gca,'fontWeight','bold','Xtick',[0 5 10 15 20 25 30 35 40 45]);
% xlim([0 45])
% xlabel('alonshore distance from W to E (km)')
% ylabel('Time weighted sediment flux for 3 months Feb-May,2014')
% title('Alongshore sediment flux --negative is ~westward, Ashton results') 
% saveFigure(gcf,'sed_flux_Ashton.jpg')

% colormap(jet)
% pcolorjw(alongshore(:)',date_time(:),flux) 
% set(gca,'fontWeight','bold','Xtick',[0 5 10 15 20 25 30 35 40 45 50]);
% xlabel('alonshore distance from W to E (km)')
% datetick('y',6)
% date_frmt='dd-mm-yyyy';
% beginTime='01-02-2014'; finishTime='04-05-2014'; 
% title('Alongshore sediment flux --negative is ~westward from Forecast Ashton form.') 
% set(gca,'YLim',[datenum(beginTime,date_frmt),datenum(finishTime,date_frmt)]) 
% caxis([-10 10])
% colorbar
% saveFigure(gcf,'Ashton formulation.jpg')
