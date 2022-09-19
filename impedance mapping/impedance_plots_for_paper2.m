%% Impedance mapping script for spinal cord grids
% Original author: Youngbin Tchoe
% Modified by: Samantha Russman
clear all
close all

% Change the following two paths:
addpath '/Users/samantha/Downloads';
addpath '/Users/samantha/Downloads/Impedance files'% path to impedance file
addpath '/Users/samantha/Downloads/electrodeimpedancemappingcode_SC'; % path to folder

%% load mapping filesv
load('impMagOrder');
load('eaglePadOrder');
load('padCoords');
load('diameters_neuromon_Sam.mat');
elecCoords = diameters_neuromon_Sam(:,1:2);
diams = diameters_neuromon_Sam(:,3);
load('ribbonCoords');

%% subplot e
filepath = 'neuromon_box2_2_imp_shorted.csv'; %neuromon_box2_1_imp.csv: 3/15,
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter = 0;
 
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];

 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter = counter+1;
         lowImpedanceMag(counter) = impedanceMag(i);
     end
 end
c = log10(lowImpedanceMag);
dat_subj4 = c;
avg = mean(lowImpedanceMag)
stdev = std(lowImpedanceMag)
counter

% Histograms of low impedances

f = figure('position', [100,100,1300,650]);

movegui('center');

edges = [4:0.05:5.3979];
h1 = histogram(c, edges);
h1.FaceColor = [0.3010 0.7450 0.9330];
h1.EdgeColor = 'none'
h1.FaceAlpha=0.5
h1.EdgeColor=[0.3010 0.7450 0.9330]
h1.LineWidth = 2
xticks([2 3 4 5 6 7])
xticklabels({'100','1K','10K','100K', '1M', '10M'})
xticks([4 4.6990 5 5.3979])
xticklabels({'10K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
set(gca,'Fontweight','Bold')

xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')
subj4_pre1 = c;
%% subplot f

% in saline
filepath = 'neuromon_box1_3_imp_ground_open.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter = counter+1;
         lowImpedanceMag(counter) = impedanceMag(i);
     end
 end
c = log10(lowImpedanceMag);
avg = mean(lowImpedanceMag)
stdev = std(lowImpedanceMag)
counter
movegui('center');

filepath = 'neuromon_box2_1_imp.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter2 = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter2 = counter2+1;
         lowImpedanceMag(counter2) = impedanceMag(i);
     end
 end
c2 = log10(lowImpedanceMag);
avg2 = mean(lowImpedanceMag)
stdev2 = std(lowImpedanceMag)
counter2
pre_1 = c2;
av_1 = avg2;
std_1 = stdev2;

filepath = 'neuromon_box2_2_impB.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter3 = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = 'r';
     else
         counter3 = counter3+1;
         lowImpedanceMag(counter3) = impedanceMag(i);
     end
 end
c3 = log10(lowImpedanceMag);
avg3 = mean(lowImpedanceMag)
stdev3 = std(lowImpedanceMag)
counter3

% Plotting
figure()
edges = [4:0.05:5.3979];
%h1 = histogram(c, edges);
%h1.FaceColor = 'r';
%h1.EdgeColor = 'none'
hold on
h2 = histogram(c2, edges);
h2.FaceColor = [0.3010 0.7450 0.9330];
h2.EdgeColor = 'none'
hold on
h3 = histogram(c3, edges);
h3.FaceColor = [0.4660 0.6740 0.1880];
h3.EdgeColor = 'none'
%h1.FaceAlpha=0.5
h2.FaceAlpha=0.5
h3.FaceAlpha=0.5
%h1.EdgeColor='r'
h2.EdgeColor = [0.3010 0.7450 0.9330]
h3.EdgeColor = [0.4660 0.6740 0.1880]
%h1.LineWidth = 2
h2.LineWidth = 2
h3.LineWidth = 2
xticks([4 4.6990 5 5.3979])
xticklabels({'10K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')
set(gca,'Fontweight','Bold')

%% post-sterilization
filepath = 'epidural_pre_resection_imp_saline.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter = counter+1;
         lowImpedanceMag(counter) = impedanceMag(i);
     end
 end
c = log10(lowImpedanceMag);
avg = mean(lowImpedanceMag)
stdev = std(lowImpedanceMag)
counter
movegui('center');

filepath = 'Neuromon3_pre-imp.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter2 = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter2 = counter2+1;
         lowImpedanceMag(counter2) = impedanceMag(i);
     end
 end
c2 = log10(lowImpedanceMag);
avg2 = mean(lowImpedanceMag)
stdev2 = std(lowImpedanceMag)
counter2

filepath = 'Neuromon_3_pre-implant.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter3 = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = 'r';
     else
         counter3 = counter3+1;
         lowImpedanceMag(counter3) = impedanceMag(i);
     end
 end
c3 = log10(lowImpedanceMag);
avg3 = mean(lowImpedanceMag)
stdev3 = std(lowImpedanceMag)
counter3
pre_2 = c2;
av_2 = avg2;
std_2 = stdev2;

% Plotting
figure()
edges = [4:0.05:5.3979];
h1 = histogram(c, edges);
h1.FaceColor = 'r';
h1.EdgeColor = 'none'
hold on
h2 = histogram(c2, edges);
h2.FaceColor = [0.3010 0.7450 0.9330];
h2.EdgeColor = 'none'
hold on
h3 = histogram(c3, edges);
h3.FaceColor = [0.4660 0.6740 0.1880];
h3.EdgeColor = 'none'
h1.FaceAlpha=0.5
h2.FaceAlpha=0.5
h3.FaceAlpha=0.5
h1.EdgeColor='r'
h2.EdgeColor = [0.3010 0.7450 0.9330]
h3.EdgeColor = [0.4660 0.6740 0.1880]
h1.LineWidth = 2
h2.LineWidth = 2
h3.LineWidth = 2
xticks([4 4.6990 5 5.3979])
xticklabels({'10K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')
set(gca,'Fontweight','Bold')

%% in tissue
filepath = 'epidural_pre_resection_imp_tissue.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter = counter+1;
         lowImpedanceMag(counter) = impedanceMag(i);
     end
 end
c = log10(lowImpedanceMag);
avg = mean(lowImpedanceMag)
stdev = std(lowImpedanceMag)
counter
movegui('center');

filepath = 'Neuromon3_Epidural.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter2 = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter2 = counter2+1;
         lowImpedanceMag(counter2) = impedanceMag(i);
     end
 end
c2 = log10(lowImpedanceMag);
avg2 = mean(lowImpedanceMag)
stdev2 = std(lowImpedanceMag)
counter2

filepath = 'pre-epidural2.csv';
[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
 colorVect = zeros(1027,3);
 elecX = elecCoords(:,1);
 elecY = elecCoords(:,2);
 counter3 = 0;
 indexMap = [];
for i=1:1024
   if elecY(i)<1500 
       indexMap = [indexMap;i];
   end
end
impedanceMag(indexMap) = [];
elecX(indexMap) = [];
elecY(indexMap) = [];
 lowImpedanceMag = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)>120e3 | impedanceMag(i)<10e3)% || impedanceMag(i)<10e3)
        %assume open or short
        colorVect(i,:) = 'r';
     else
         counter3 = counter3+1;
         lowImpedanceMag(counter3) = impedanceMag(i);
     end
 end
c3 = log10(lowImpedanceMag);
avg3 = mean(lowImpedanceMag)
stdev3 = std(lowImpedanceMag)
counter3
pre_3 = c2;
av_3 = avg2;
std_3 = stdev2;

% Plotting
figure()
edges = [4:0.05:5.3979];
h1 = histogram(c, edges);
h1.FaceColor = 'r';
h1.EdgeColor = 'none'
hold on
h2 = histogram(c2, edges);
h2.FaceColor = [0.3010 0.7450 0.9330];
h2.EdgeColor = 'none'
hold on
h3 = histogram(c3, edges);
h3.FaceColor = [0.4660 0.6740 0.1880];
h3.EdgeColor = 'none'
h1.FaceAlpha=0.5
h2.FaceAlpha=0.5
h3.FaceAlpha=0.5
h1.EdgeColor='r'
h2.EdgeColor = [0.3010 0.7450 0.9330]
h3.EdgeColor = [0.4660 0.6740 0.1880]
h1.LineWidth = 2
h2.LineWidth = 2
h3.LineWidth = 2
xticks([4 4.6990 5 5.3979])
xticklabels({'10K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')
set(gca,'Fontweight','Bold')

%% New impedance plot

figure()
subplot(1,3,1)
edges = [4:0.05:5.3979];
h1 = histogram(pre_1, edges);
title('As fabricated, saline')
h1.EdgeColor = 'none'
xticks([4 4.3010 4.6990 5 5.3979])
xticklabels({'10K','20K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
%xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')
set(gca,'Fontweight','Bold')
subplot(1,3,2)
h2 = histogram(pre_2, edges);
title('Post-sterilization, saline')
h1.FaceColor = [0.3010 0.7450 0.9330];
h2.EdgeColor = 'none'
xticks([4 4.3010 4.6990 5 5.3979])
xticklabels({'10K','20K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
xlabel('Impedance Magnitude @ 1 kHz')
%ylabel('Counts')
set(gca,'Fontweight','Bold')
subplot(1,3,3)
h3 = histogram(pre_3, edges);
title('Post-sterilization, tissue')
h3.FaceColor = [0.4660 0.6740 0.1880];
h3.EdgeColor = 'none'
h1.FaceAlpha=0.5
h2.FaceAlpha=0.5
h3.FaceAlpha=0.5
h2.EdgeColor='r'
h2.FaceColor = 'r';
h1.EdgeColor = [0.3010 0.7450 0.9330]
h3.EdgeColor = [0.4660 0.6740 0.1880]
h1.LineWidth = 2
h2.LineWidth = 2
h3.LineWidth = 2
xticks([4 4.3010 4.6990 5 5.3979])
xticklabels({'10K','20K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
%xlabel('Impedance Magnitude @ 1 kHz')
%ylabel('Counts')
set(gca,'Fontweight','Bold')
set(gcf,'Position',[306         394        1185         411])

%% New impedance plot 2/24

figure()
subplot(1,3,1)
edges = [4:0.05:5.3979];
h1 = histogram(pre_1, edges);
title('As fabricated, saline')
h1.EdgeColor = 'none'
xticks([4 4.3010 4.6990 5 5.3979])
xticklabels({'10K','20K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
%xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')
set(gca,'Fontweight','Bold')
subplot(1,3,2)
h2 = histogram(pre_2, edges);
title('Post-sterilization, saline')
h1.FaceColor = [0.3010 0.7450 0.9330];
h2.EdgeColor = 'none'
xticks([4 4.3010 4.6990 5 5.3979])
xticklabels({'10K','20K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
xlabel('Impedance Magnitude @ 1 kHz')
%ylabel('Counts')
set(gca,'Fontweight','Bold')
subplot(1,3,3)
h3 = histogram(pre_3, edges);
title('Post-sterilization, tissue')
h3.FaceColor = [0.4660 0.6740 0.1880];
h3.EdgeColor = 'none'
h1.FaceAlpha=0.5
h2.FaceAlpha=0.5
h3.FaceAlpha=0.5
h2.EdgeColor='r'
h2.FaceColor = 'r';
h1.EdgeColor = [0.3010 0.7450 0.9330]
h3.EdgeColor = [0.4660 0.6740 0.1880]
h1.LineWidth = 2
h2.LineWidth = 2
h3.LineWidth = 2
xticks([4 4.3010 4.6990 5 5.3979])
xticklabels({'10K','20K','50K','100K', '250K'})
set(gcf,'Position',[680   218   366   567])
box off
set(gca,'Fontsize',20)
%xlabel('Impedance Magnitude @ 1 kHz')
%ylabel('Counts')
set(gca,'Fontweight','Bold')
set(gcf,'Position',[306         394        1185         411])

%% Welch's test to compare pre_1 and pre_2
dat1=10.^(pre_1);
dat2=10.^(pre_2);
[h,p,ci,stats] = ttest2(dat1,dat2,'Vartype','unequal') % Welch's test


