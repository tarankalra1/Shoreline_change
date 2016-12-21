% This code is to group all the data from various months 
clear all ; close all ; clc ;

months={'04_01','07_04','10_07','12_10'};
year=1997; 
path_str='/home/gadar/Documents/Development/Q_deep_FireIsland/Output/Q_output_1997'

ntotal=0;

for m=1:length(months)
    
  year_month_str=strcat(num2str(year),'_',months(m)); 
  file_get=char(strcat('Q_',year_month_str,'.mat')); 
 
  filename_Q=strcat(path_str,'/',file_get)  
  
  load(filename_Q)
  [jmax,tmax]=size(Q);
  
  for t=1:tmax
    for j=1:jmax
      Q_Ashton_1997(j,t+ntotal)=Q(j,t);
      anglerel_Ashton_1997(j,t+ntotal)=angle_rel(j,t);
      time_Ashton_1997(t+ntotal)=jul_times(t); 
    end
  end
  ntotal=ntotal+tmax ;
end 
save('/media/gadar/DATADRIVE2/CSI/Ashton/Q_Ashton_1997.mat','Q_Ashton_1997',......
                'anglerel_Ashton_1997','time_Ashton_1997')
            
            
 
