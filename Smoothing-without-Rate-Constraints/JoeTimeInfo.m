function [MeanUpdate,StdDevUpdate,MeanHorizon,StdDevHorizon] = JoeTimeInfo(tauvec, Tvec,time)
%JOETIMEINFO: Calculates information about the mean and standard deviation 
%of the distribution of time horizons and the updates times.

%% Update Times

updates= [];
for i=1:length(tauvec)-1
    updates(i) = tauvec(i+1)-tauvec(i);
end
updates = updates*time;

MeanUpdate = mean(updates); %in hours
StdDevUpdate = std(updates);
    
% figure
% histogram(updates,'BinWidth',0.5);
% title('Histogram of the difference between updates times at each step in the algorithm');
% xlabel('Difference between updates times (h)')
% ylabel('Frequency')


%% Time Horizons

%tauvec is always one index bigger than Tvec, since tauvec includes the
%final time value. We remove this.
tauvec = tauvec(1:length(tauvec)-1);

horizons = (Tvec-tauvec)*time; %in hours

% figure
% histogram(horizons,'BinWidth',0.5);
% title('Histogram of the time ahead required of energy predictions to complete each step in the algorithm');
% xlabel('Time ahead of energy predictions required (h)')
% ylabel('Frequency')

MeanHorizon = mean(horizons);
StdDevHorizon = std(horizons);

end

