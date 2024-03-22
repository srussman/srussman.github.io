% This code (1) performs spectral analysis on SSEPs, (2) generates heatmaps
% of the recorded SSEPs, and (3) calculates SSEP velocity


% Author: srussman@ucsd.edu

%% Spectral analysis
% a. spectrogram of responses
% b. example response
load('post_below_LU_35_filtered_30-300') % load stimulation segment
load('post_below_LU_35_filtered_trials') % load stimulation trial timepoints

% Spectrogram
data = data_orig;
nfft = [1:5:600];
clear powers
%for ch=1:10 %length(goodchs)
%
ch=100; % channel of choice
for i = 1:length(locs)
[S,F,T,P] = spectrogram(data(ch,locs(i)-50*20:locs(i)+200*20),500,400,nfft,20000, 'yaxis');
for k=1:length(P(:,1))
   logged_P = 10*log10(abs(P));
   db = logged_P(k,1:8);
   base_db(k) = mean(logged_P(k,1:8));
   P(k,:) = logged_P(k,:)-base_db(k); % frequency vs. time
end
powers(i,:,:)=P;
end
av_powers_ch = squeeze(mean(powers,1));

figure
hold on
% panel a
subplot(1,2,1)
h = imagesc((T-0.05).*1000,F,av_powers_ch,[0 6]); % offset time of stim at 10 ms
c = colorbar;
c.Color = 'k';
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Spectrogram of response from ch. 100')
set(gca,'FontSize',16)
c.Label.String = 'Power (dB/Hz)';
c.Label.FontSize = 16
              
load('post_below_LU_35_data_30-300')

% panel b
subplot(1,2,2)
hold on
plot((1:2201)./20-10,averages(100,:),'b','Linewidth',2)
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title('Response from ch. 100')
plot([0 0], [-10 10],'k')
set(gca,'FontSize',16)

set(gcf,'Position', [160 585 1094 393])


%% Heatmap plots
% Panels b-d: left side
load('subdural2_goodchs') % load low impedance channels
addpath '/Users/samantha/Downloads'
load('subdural2_LU_50_averages_30-300') % load averaged data
addpath '/Users/samantha/Downloads/electrodeimpedancemappingcode_SC';
load('padCoords');
load('diameters_neuromon_Sam.mat');
elecCoords = diameters_neuromon_Sam(:,1:2);
diams = diameters_neuromon_Sam(:,3);
elecX = elecCoords(:,1);
elecY = elecCoords(:,2);

% Plotting spatial responses
goodchs2=goodchs;
indexMap = [];
for i=1:length(goodchs)
   if elecY(goodchs(i))<1500 
       indexMap = [indexMap;i];
   end
end
goodchs2(indexMap) = [];
averages2=averages;
averages2(indexMap,:) =[];
close all

% Run through multiple timesteps after stimulation and plot heatmaps
for t=[10:0.5:15]
figure()
hold on
for timestep = t*20
    s=scatter(-elecX(goodchs2)./1000,-elecY(goodchs2)./1000,600,averages2(:,t+10*20))
    axis([-2.500 2.500 -12.500 -1.600])
    s.Marker = 's'
    set(gca,'FontSize',24)
    color_Range = [0  2];
s.MarkerEdgeColor = 'flat';
s.MarkerFaceColor = 'flat';
    s.LineWidth = 0.75;
    caxis(color_Range)
    c = colorbar;
    c.Color = 'k';
    %title(['Response at ' num2str( timestep./20-10 ) 'ms'], 'FontSize', 14 )
    set(gcf,'position',[680,117,width,height])
    %c.Label.String = 'Voltage \muV';
    %c.Label.FontSize = 24
    xlabel('x-position (mm)','FontSize', 14)
    ylabel('y-position (mm)','FontSize', 14)
end
set(gca,'Fontsize',14)
set(gcf,'Position',[680   218   366   567])
set(gca,'XColor', 'none','YColor','none')
end





%% Midline mapping from 03-15-2021 case
% panel b: plot scatter overlay of left and right responses
clear all
close all
figure()
load('post_below_LU_35_data_30-300')

addpath '/Users/samantha/Downloads/electrodeimpedancemappingcode_SC';
load('padCoords');
load('diameters_neuromon_Sam.mat');
elecCoords = diameters_neuromon_Sam(:,1:2);
diams = diameters_neuromon_Sam(:,3);
elecX = elecCoords(:,1);
elecY = elecCoords(:,2);
goodchs2=goodchs;
indexMap = [];
for i=1:length(goodchs)
   if elecY(goodchs(i))<1500 
       indexMap = [indexMap;i];
   end
end
goodchs2(indexMap) = [];
averages2=averages;
averages2(indexMap,:) =[];

offset_x = 20*40;
offset_y = 12;
averages30 = averages2(:,12*20:42*20-1);
subplot(2,2,[1 3])
hold on
for i=1:length(goodchs2)
    ch = goodchs2(i)
            multiple_x = diameters_neuromon_Sam(ch,1)/350;
            multiple_y = diameters_neuromon_Sam(ch,2)/400;
            plot((1:20*30)+multiple_x*offset_x,averages30(i,:)+multiple_y*offset_y,'Color',[0 0 0]+0.05*2,'Linewidth',2)
end
set(gcf,'Position',[309 1 488 984])
xlabel('x-position (separation 350 um)')
ylabel('y-position (separation 400 um)')
xlim([-4500 5100])
ylim([40 420])
title('Overlay of left and right upper responses')

clear all
load('post_below_RU_30_data_30-300')
for i=1:length(goodchs)
    if find(averages(i,:)>20)>0
        averages(i,:)=zeros(1,2201);
    end
end

%% Velocity calculation
% Both SSEPs and MEPs around 50 m/s based on literature
% https://www.sciencedirect.com/science/article/pii/S0022510X07004777

load('subdural1_LU_30_data_30-300');
load('diameters_neuromon_Sam.mat');
elecCoords = diameters_neuromon_Sam(:,1:2);
diams = diameters_neuromon_Sam(:,3);
elecX = elecCoords(:,1);
elecY = elecCoords(:,2);

for i=1:length(goodchs)
    if find(averages(i,:)>12)>0
        averages(i,:)=zeros(1,2201);
    end
end
goodchs2=goodchs;
indexMap = [];
for i=1:length(goodchs)
   if elecY(goodchs(i))<1500 
       indexMap = [indexMap;i];
   end
end
goodchs2(indexMap) = [];
averages2=averages;
averages2(indexMap,:) =[];

% Compare time delays is rightmost column from 1st electrode
idx1=find(goodchs2==20); % 1st electrode

idxs=find(elecX(goodchs2)==elecX(goodchs2(idx1)));

for i=1:length(idxs)
D(i) = finddelay(averages2(idxs(i),17*20:26*20),averages2(idx1,17*20:26*20));
delay(i)=D(i)./20; % in ms
end

figure()
plot(elecY(goodchs2(idxs)),delay,'k*')
xlabel('Electrode y-position (\mum)')
ylabel('Delay (ms)')
set(gca,'FontSize',14)
