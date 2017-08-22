function [TotalVariation,TotalVariationwind,MeanUpdate,StdDevUpdate,MeanHorizon,StdDevHorizon] = JoeCapacity(a1,a2,s)
%JOECAPACITY: We vary the capacity of the battery and see its effects on
%the total variation, mean/stdev of time horizons and update times, and on
%max charge/discharge required

tic

M=2:2:700; %varying the store capacities in kWh

for i=1:length(M);
    [TotalVariation(i),TotalVariationwind(i),MeanUpdate(i),StdDevUpdate(i),MeanHorizon(i),StdDevHorizon(i),Charge(i),Discharge(i)]=JoeOptimiser(M(i),a1,a2,s);
end

figure
plot([0,M],[max(TotalVariationwind),TotalVariation]);
title('Effect of Capacity on Total Variaton');
xlabel('Capacity (kWh)')
ylabel('Total Variation (kWh)')

figure
plot(M,MeanUpdate);
title('Effect of Capacity on Mean Update Time');
xlabel('Capacity (kWh)')
ylabel('Mean Update Time (kWh)')

figure
plot(M,StdDevUpdate);
title('Effect of Capacity on Standard Deviation of Update Times');
xlabel('Capacity (kWh)')
ylabel('Standard Deviation (kWh)')

figure
plot(M,MeanHorizon);
title('Effect of Capacity on Mean Time Horizon');
xlabel('Capacity (kWh)')
ylabel('Mean Time Horizon (kWh)')

figure
plot(M,StdDevHorizon);
title('Effect of Capacity on Standard Deviation of Time Horizons');
xlabel('Capacity (kWh)')
ylabel('Standard Deviation (kWh)')

figure
plot(M,Charge);
title('Effect of Capacity on Max Charging Power Required');
xlabel('Capacity (kWh)')
ylabel('Power (kW)')

figure
plot(M,Discharge);
title('Effect of Capacity on Max Discharging Power Required');
xlabel('Capacity (kWh)')
ylabel('Power (kW)')

toc
end