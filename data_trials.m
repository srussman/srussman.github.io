%% Find stim points
clear all
close all
load('subdural1_LU_15_filtered_30-300')

[pks,locs] = findpeaks(data_orig(100,:),'MinPeakProminence',55,'MinPeakDistance',3000);

%% Remove any additional incorrect peaks
pks([28])=[];
locs([28])=[];

%% Plot all captured peaks
figure()
hold on
plot(t,data_orig(100,:))
plot(t(locs),pks,'o')

%% Save file
save('subdural1_LU_15_trials','pks','locs')

