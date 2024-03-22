% This codes performs the multidiameter analysis

% Author: srussman@ucsd.edu

%% SUBDURAL
%%
close all
% Load data
load('epidural_goodchs')
load('epidural_LU_30_data_30-300')

% Load diameter map and find good diameters
load('diameters_neuromon_Sam.mat');
good_diams = diameters_neuromon_Sam(goodchs,:);

indices=find(good_diams(:,2)<1500);
diam_averages = averages(indices,:);
diams = good_diams(indices,:);

% Loop through all diameters in multidiameter array
for i=1:length(diams(:,1))
    d30=find(diams(:,3)==30);
    d60=find(diams(:,3)==60);
    d80=find(diams(:,3)==80);
    d110=find(diams(:,3)==110);
    d140=find(diams(:,3)==140);
    d170=find(diams(:,3)==170);
    d210=find(diams(:,3)==210);
    d240=find(diams(:,3)==240);
    d480=find(diams(:,3)==480);
end

for i=1:length(d30)
    differences30(i) = max(diam_averages(d30(i),16*20:36*20),[],2) - min(diam_averages(d30(i),16*20:36*20),[],2); 
end
for i=1:length(d60)
    differences60(i) = max(diam_averages(d60(i),16*20:36*20),[],2) - min(diam_averages(d60(i),16*20:36*20),[],2); 
end
for i=1:length(d80)
    differences80(i) = max(diam_averages(d80(i),16*20:36*20),[],2) - min(diam_averages(d80(i),16*20:36*20),[],2); 
end
for i=1:length(d110)
    differences110(i) = max(diam_averages(d110(i),16*20:36*20),[],2) - min(diam_averages(d110(i),16*20:36*20),[],2); 
end
for i=1:length(d140)
    differences140(i) = max(diam_averages(d140(i),16*20:36*20),[],2) - min(diam_averages(d140(i),16*20:36*20),[],2); 
end
for i=1:length(d170)
    differences170(i) = max(diam_averages(d170(i),16*20:36*20),[],2) - min(diam_averages(d170(i),16*20:36*20),[],2); 
end
for i=1:length(d210)
    differences210(i) = max(diam_averages(d210(i),16*20:36*20),[],2) - min(diam_averages(d210(i),16*20:36*20),[],2); 
end
for i=1:length(d240)
    differences240(i) = max(diam_averages(d240(i),16*20:36*20),[],2) - min(diam_averages(d240(i),16*20:36*20),[],2); 
end
for i=1:length(d480)
    differences480(i) = max(diam_averages(d480(i),16*20:36*20),[],2) - min(diam_averages(d480(i),16*20:36*20),[],2); 
end
set(groot,'defaultLineMarkerSize',45);

figure()
hold on
plot(repmat(30,1,6),differences30,'.')
plot(repmat(60,1,6),differences60,'.')
plot(repmat(80,1,6),differences80,'.')
plot(repmat(110,1,6),differences110,'.')
plot(repmat(140,1,6),differences140,'.')
plot(repmat(170,1,6),differences170,'.')
plot(repmat(210,1,5),differences210,'.')
plot(repmat(240,1,6),differences240,'.')
plot(repmat(480,1,6),differences480,'.')
xlabel('Diameter (\mum)')
ylabel('Raw peak-to-peak amplitude (\muV)')
set(gca,'Fontsize',32)
set(gca,'linewidth',4)

figure()
hold on
plot(repmat(30,1,3),differences30(1:3),'k.')
plot(repmat(60,1,3),differences60(1:3),'k.')
plot(repmat(80,1,3),differences80(1:3),'k.')
plot(repmat(110,1,3),differences110(1:3),'k.')
plot(repmat(140,1,3),differences140(1:3),'k.')
plot(repmat(170,1,3),differences170(1:3),'k.')
plot(repmat(210,1,3),differences210(1:3),'k.')
plot(repmat(240,1,3),differences240(1:3),'k.')
plot(repmat(480,1,3),differences480(1:3),'k.')
plot(repmat(30,1,3),differences30(4:6),'b.')
plot(repmat(60,1,3),differences60(4:6),'b.')
plot(repmat(80,1,3),differences80(4:6),'b.')
plot(repmat(110,1,3),differences110(4:6),'b.')
plot(repmat(140,1,3),differences140(4:6),'b.')
plot(repmat(170,1,3),differences170(4:6),'b.')
plot(repmat(210,1,2),differences210(4:5),'b.')
plot(repmat(240,1,3),differences240(4:6),'b.')
plot(repmat(480,1,3),differences480(4:6),'b.')
xlabel('Diameter (\mum)')
ylabel('Raw peak-to-peak amplitude (\muV)')
set(gca,'Fontsize',32)
set(gca,'linewidth',4)
set(gcf,'Position',[680   479   691   499])

figure()
hold on
m=15;
e1=errorbar(repmat(30,1,1),mean(differences30(1:3)),std(differences30(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(60,1,1),mean(differences60(1:3)),std(differences60(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(80,1,1),mean(differences80(1:3)),std(differences80(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(110,1,1),mean(differences110(1:3)),std(differences110(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(140,1,1),mean(differences140(1:3)),std(differences140(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(170,1,1),mean(differences170(1:3)),std(differences170(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1=errorbar(repmat(210,1,1),mean(differences210(1:3)),std(differences210(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1=errorbar(repmat(240,1,1),mean(differences240(1:3)),std(differences240(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(480,1,1),mean(differences480(1:3)),std(differences480(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(30,1,1),mean(differences30(4:6)),std(differences30(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(60,1,1),mean(differences60(4:6)),std(differences60(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(80,1,1),mean(differences80(4:6)),std(differences80(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(110,1,1),mean(differences110(4:6)),std(differences110(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(140,1,1),mean(differences140(4:6)),std(differences140(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(170,1,1),mean(differences170(4:6)),std(differences170(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(210,1,1),mean(differences210(4:5)),std(differences210(4:5)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(240,1,1),mean(differences240(4:6)),std(differences240(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(480,1,1),mean(differences480(4:6)),std(differences480(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
xlabel('Diameter (\mum)')
ylabel('Peak-to-peak amplitude (\muV)')
set(gca,'Fontsize',32)
set(gca,'linewidth',4)
set(gcf,'Position',[680   479   691   499])
ylim([3 8])

figure()
hold on
plot(repmat(30,1,1),mean(differences30(1:3))+1,'b.')
plot(repmat(60,1,1),mean(differences60(1:3))+2,'b.')
plot(repmat(80,1,1),mean(differences80(1:3))+2,'b.')
plot(repmat(110,1,1),mean(differences110(1:3))+0.5,'b.')
plot(repmat(140,1,1),mean(differences140(1:3))+0.5,'b.')
plot(repmat(170,1,1),mean(differences170(1:3))+0.5,'b.')
plot(repmat(210,1,1),mean(differences210(1:3))+1.5,'b.')
plot(repmat(240,1,1),mean(differences240(1:3))+1.5,'b.')
plot(repmat(480,1,1),mean(differences480(1:3))+3,'b.')
plot(repmat(30,1,1),mean(differences30(4:6)),'r.')
plot(repmat(60,1,1),mean(differences60(4:6)),'r.')
plot(repmat(80,1,1),mean(differences80(4:6)),'r.')
plot(repmat(110,1,1),mean(differences110(4:6)),'r.')
plot(repmat(140,1,1),mean(differences140(4:6)),'r.')
plot(repmat(170,1,1),mean(differences170(4:6)),'r.')
plot(repmat(210,1,1),mean(differences210(4:5)),'r.')
plot(repmat(240,1,1),mean(differences240(4:6)),'r.')
plot(repmat(480,1,1),mean(differences480(4:6)),'r.')
xlabel('Diameter (\mum)')
ylabel('Adjusted mean ptp amplitude (\muV)')
set(gca,'Fontsize',32)

figure()
hold on
for i=1:length(diam_averages(:,1))
   plot((1:2201)/20-10,diam_averages(i,:)) 
end

%% Plot baseline and SNR

load('epidural_LU_30_filtered_30-300')
load('epidural_LU_30_filtered_trials')
%averages2=averages(:,16*20:36*20-1);

data_30all = data;
data_30all(indices,:)=[];
data_diams = data(indices,:);
baseline = zeros(length(indices),400);
lb = 400;
for ch=1:length(indices)
    baseline_trials = zeros(length(locs),lb);
    for l=1:length(locs)
        baseline_trials(l,:) = data_diams(ch,locs(l)-300-lb+1:locs(l)-300); % Take 20 ms segment of baseline before start of each stim
    end
    baseline(ch,:)=mean(baseline_trials);
end
goodchs2 = goodchs;
goodchs2(indices)=[];
baseline_30all = zeros(length(goodchs2),400);
lb = 400;
for ch=1:length(goodchs2)
    baseline_trials = zeros(length(locs),lb);
    for l=1:length(locs)
        baseline_trials(l,:) = data_30all(ch,locs(l)-300-lb+1:locs(l)-300); % Take 20 ms segment of baseline before start of each stim
    end
    baseline_30all(ch,:)=mean(baseline_trials);
end

figure()
title('Standard deviation of baseline')
hold on
plot(diams(1:53,3),mean(std(baseline(1:53),0,2)),'.')
plot(diams(1,3),mean(std(baseline_30all,0,2)),'.')
xlabel('Diameter (\mum)')
ylabel('Standard deviation (\muV)')
set(gca,'Fontsize',32)

av_baseline = mean(baseline,2);

for i=1:length(diams(:,1))
    d30=find(diams(:,3)==30);
    d60=find(diams(:,3)==60);
    d80=find(diams(:,3)==80);
    d110=find(diams(:,3)==110);
    d140=find(diams(:,3)==140);
    d170=find(diams(:,3)==170);
    d210=find(diams(:,3)==210);
    d240=find(diams(:,3)==240);
    d480=find(diams(:,3)==480);
end


baseline30 = mean(std(baseline(d30(1:6),:),0,2));
baseline60 = mean(std(baseline(d60(1:6),:),0,2));
baseline80 = mean(std(baseline(d80(1:6),:),0,2));
baseline110 = mean(std(baseline(d110(1:6),:),0,2));
baseline140 = mean(std(baseline(d140(1:6),:),0,2));
baseline170 = mean(std(baseline(d170(1:6),:),0,2));
baseline210 = mean(std(baseline(d210(1:5),:),0,2));
baseline240 = mean(std(baseline(d240(1:6),:),0,2));
baseline480 = mean(std(baseline(d480(1:6),:),0,2));
baseline_30all_std = mean(std(baseline_30all,0,2));

baseline30std = std(std(baseline(d30(1:6),:),0,2));
baseline60std = std(std(baseline(d60(1:6),:),0,2));
baseline80std = std(std(baseline(d80(1:6),:),0,2));
baseline110std = std(std(baseline(d110(1:6),:),0,2));
baseline140std = std(std(baseline(d140(1:6),:),0,2));
baseline170std = std(std(baseline(d170(1:6),:),0,2));
baseline210std = std(std(baseline(d210(1:5),:),0,2));
baseline240std = std(std(baseline(d240(1:6),:),0,2));
baseline480std = std(std(baseline(d480(1:6),:),0,2));
baseline_30all_stdstd = std(std(baseline_30all,0,2));

set(groot,'defaultLineMarkerSize',45);
figure()
hold on
%plot(30,baseline30,'o')
e1=errorbar(60,baseline60,baseline60std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(80,baseline80,baseline80std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(110,baseline110,baseline110std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(140,baseline140,baseline140std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(170,baseline170,baseline170std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(210,baseline210,baseline210std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(240,baseline240,baseline240std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(480,baseline480,baseline480std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(30,baseline_30all_std,baseline_30all_stdstd,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
xlabel('Diameter (\mum)')
ylabel('Mean baseline standard deviation (\muV)')
set(gca,'Fontsize',32)
set(gca,'linewidth',3)
set(gcf,'Position',[680   479   691   499])
ylim([0.1 0.6])

% Perform t-test
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d60(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d80(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d110(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d140(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d170(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d210(1:5),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d240(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d480(1:6),:),0,2))

%% Power spectral density

rng default
Fs = 20000;
t = 0:1/Fs:1-1/Fs;

[psdx30,freq30]=psdx_calc(mean(diam_averages(d30(1:3),16*20:56*20)));
[psdx60,freq60]=psdx_calc(mean(diam_averages(d60(1:3),16*20:56*20)));
[psdx80,freq80]=psdx_calc(mean(diam_averages(d80(1:3),16*20:56*20)));
[psdx110,freq110]=psdx_calc(mean(diam_averages(d110(1:3),16*20:56*20)));
[psdx140,freq140]=psdx_calc(mean(diam_averages(d140(1:3),16*20:56*20)));
[psdx170,freq170]=psdx_calc(mean(diam_averages(d170(1:3),16*20:56*20)));
[psdx210,freq210]=psdx_calc(mean(diam_averages(d210(1:3),16*20:56*20)));
[psdx240,freq240]=psdx_calc(mean(diam_averages(d240(1:3),16*20:56*20)));
[psdx480,freq480]=psdx_calc(mean(diam_averages(d480(1:3),16*20:56*20)));

addpath '/Users/samantha/Downloads/github_repo'
cmap=colormap(brewermap(9,'rdylbu'));
figure()
hold on
plot(freq30,10*log10(psdx30),'.-','Color',cmap(1,:),'Linewidth',3)
plot(freq60,10*log10(psdx60),'.-','Color',cmap(2,:),'Linewidth',3)
plot(freq80,10*log10(psdx80),'.-','Color',cmap(3,:),'Linewidth',3)
plot(freq110,10*log10(psdx110),'.-','Color',cmap(4,:),'Linewidth',3)
plot(freq140,10*log10(psdx140),'.-','Color',cmap(5,:),'Linewidth',3)
plot(freq170,10*log10(psdx170),'.-','Color',cmap(6,:),'Linewidth',3)
plot(freq210,10*log10(psdx210),'.-','Color',cmap(7,:),'Linewidth',3)
plot(freq240,10*log10(psdx240),'.-','Color',cmap(8,:),'Linewidth',3)
plot(freq480,10*log10(psdx480),'.-','Color',cmap(9,:),'Linewidth',3)
grid off
xlim([0 1000])
%title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('30 \mum','60 \mum','80 \mum','110 \mum','140 \mum','170 \mum','210 \mum','240 \mum','480 \mum')
legend boxoff
set(gca,'Fontsize',32)
set(gca,'linewidth',4)

%% EPIDURAL
% Load data
% load('epidural_LU_30_data_30-300')
% load('epidural_goodchs')
m=15;
load('subdural1_goodchs_0315')
load('subdural1_LU_30_averages_30-300')

% Load diameter map and find good diameters
load('diameters_neuromon_Sam.mat');
good_diams = diameters_neuromon_Sam(goodchs,:);

indices=find(good_diams(:,2)<1500);
diam_averages = averages(indices,:);
diams = good_diams(indices,:);

% Loop through all diameters in multidiameter array
for i=1:length(diams(:,1))
    d30=find(diams(:,3)==30);
    d60=find(diams(:,3)==60);
    d80=find(diams(:,3)==80);
    d110=find(diams(:,3)==110);
    d140=find(diams(:,3)==140);
    d170=find(diams(:,3)==170);
    d210=find(diams(:,3)==210);
    d240=find(diams(:,3)==240);
    d480=find(diams(:,3)==480);
end

for i=1:length(d30)
    differences30(i) = max(diam_averages(d30(i),16*20:36*20),[],2) - min(diam_averages(d30(i),16*20:36*20),[],2); 
end
for i=1:length(d60)
    differences60(i) = max(diam_averages(d60(i),16*20:36*20),[],2) - min(diam_averages(d60(i),16*20:36*20),[],2); 
end
for i=1:length(d80)
    differences80(i) = max(diam_averages(d80(i),16*20:36*20),[],2) - min(diam_averages(d80(i),16*20:36*20),[],2); 
end
for i=1:length(d110)
    differences110(i) = max(diam_averages(d110(i),16*20:36*20),[],2) - min(diam_averages(d110(i),16*20:36*20),[],2); 
end
for i=1:length(d140)
    differences140(i) = max(diam_averages(d140(i),16*20:36*20),[],2) - min(diam_averages(d140(i),16*20:36*20),[],2); 
end
for i=1:length(d170)
    differences170(i) = max(diam_averages(d170(i),16*20:36*20),[],2) - min(diam_averages(d170(i),16*20:36*20),[],2); 
end
for i=1:length(d210)
    differences210(i) = max(diam_averages(d210(i),16*20:36*20),[],2) - min(diam_averages(d210(i),16*20:36*20),[],2); 
end
for i=1:length(d240)
    differences240(i) = max(diam_averages(d240(i),16*20:36*20),[],2) - min(diam_averages(d240(i),16*20:36*20),[],2); 
end
for i=1:length(d480)
    differences480(i) = max(diam_averages(d480(i),16*20:36*20),[],2) - min(diam_averages(d480(i),16*20:36*20),[],2); 
end
set(groot,'defaultLineMarkerSize',45);

figure()
hold on
plot(repmat(30,1,6),differences30,'.')
plot(repmat(60,1,6),differences60,'.')
plot(repmat(80,1,6),differences80,'.')
plot(repmat(110,1,6),differences110,'.')
plot(repmat(140,1,6),differences140,'.')
plot(repmat(170,1,6),differences170,'.')
plot(repmat(210,1,5),differences210,'.')
plot(repmat(240,1,6),differences240,'.')
plot(repmat(480,1,6),differences480,'.')
xlabel('Diameter (\mum)')
ylabel('Raw peak-to-peak amplitude (\muV)')
set(gca,'Fontsize',32)

figure()
hold on
plot(repmat(30,1,3),differences30(1:3),'b.')
plot(repmat(60,1,3),differences60(1:3),'b.')
plot(repmat(80,1,3),differences80(1:3),'b.')
plot(repmat(110,1,3),differences110(1:3),'b.')
plot(repmat(140,1,3),differences140(1:3),'b.')
plot(repmat(170,1,3),differences170(1:3),'b.')
plot(repmat(210,1,3),differences210(1:3),'b.')
plot(repmat(240,1,3),differences240(1:3),'b.')
plot(repmat(480,1,3),differences480(1:3),'b.')
plot(repmat(30,1,3),differences30(4:6),'r.')
plot(repmat(60,1,3),differences60(4:6),'r.')
plot(repmat(80,1,3),differences80(4:6),'r.')
plot(repmat(110,1,3),differences110(4:6),'r.')
plot(repmat(140,1,3),differences140(4:6),'r.')
plot(repmat(170,1,3),differences170(4:6),'r.')
plot(repmat(210,1,2),differences210(4:5),'r.')
plot(repmat(240,1,3),differences240(4:6),'r.')
plot(repmat(480,1,3),differences480(4:6),'r.')
xlabel('Diameter (\mum)')
ylabel('Raw peak-to-peak amplitude (\muV)')
set(gca,'Fontsize',32)

figure()
hold on
e1=errorbar(repmat(30,1,1),mean(differences30(1:3)),std(differences30(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(60,1,1),mean(differences60(1:3)),std(differences60(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(80,1,1),mean(differences80(1:3)),std(differences80(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(110,1,1),mean(differences110(1:3)),std(differences110(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(140,1,1),mean(differences140(1:3)),std(differences140(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(170,1,1),mean(differences170(1:3)),std(differences170(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(210,1,1),mean(differences210(1:3)),std(differences210(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(240,1,1),mean(differences240(1:3)),std(differences240(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(480,1,1),mean(differences480(1:3)),std(differences480(1:3)),'vertical','bo','MarkerFaceColor','b')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(30,1,1),mean(differences30(4:6)),std(differences30(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(60,1,1),mean(differences60(4:6)),std(differences60(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(80,1,1),mean(differences80(4:6)),std(differences80(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(110,1,1),mean(differences110(4:6)),std(differences110(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(140,1,1),mean(differences140(4:6)),std(differences140(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(170,1,1),mean(differences170(4:6)),std(differences170(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(210,1,1),mean(differences210(4:5)),std(differences210(4:5)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(240,1,1),mean(differences240(4:6)),std(differences240(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(repmat(480,1,1),mean(differences480(4:6)),std(differences480(4:6)),'vertical','ro','MarkerFaceColor','r')
e1.LineWidth = 4;
e1.MarkerSize = m;
xlabel('Diameter (\mum)')
ylabel('Peak-to-peak amplitude (\muV)')
set(gca,'Fontsize',32)
set(gca,'linewidth',4)
set(gcf,'Position',[680   479   691   499])

figure()
hold on
plot(repmat(30,1,1),mean(differences30(1:3))+1,'k.')
plot(repmat(60,1,1),mean(differences60(1:3))+2,'k.')
plot(repmat(80,1,1),mean(differences80(1:3))+2,'k.')
plot(repmat(110,1,1),mean(differences110(1:3))+0.5,'k.')
plot(repmat(140,1,1),mean(differences140(1:3))+0.5,'k.')
plot(repmat(170,1,1),mean(differences170(1:3))+0.5,'k.')
plot(repmat(210,1,1),mean(differences210(1:3))+1.5,'k.')
plot(repmat(240,1,1),mean(differences240(1:3))+1.5,'k.')
plot(repmat(480,1,1),mean(differences480(1:3))+3,'k.')
plot(repmat(30,1,1),mean(differences30(4:6)),'b.')
plot(repmat(60,1,1),mean(differences60(4:6)),'b.')
plot(repmat(80,1,1),mean(differences80(4:6)),'b.')
plot(repmat(110,1,1),mean(differences110(4:6)),'b.')
plot(repmat(140,1,1),mean(differences140(4:6)),'b.')
plot(repmat(170,1,1),mean(differences170(4:6)),'b.')
plot(repmat(210,1,1),mean(differences210(4:5)),'b.')
plot(repmat(240,1,1),mean(differences240(4:6)),'b.')
plot(repmat(480,1,1),mean(differences480(4:6)),'b.')
xlabel('Diameter (\mum)')
ylabel('Adjusted mean ptp amplitude (\muV)')
set(gca,'Fontsize',32)

figure()
hold on
for i=1:length(diam_averages(:,1))
   plot((1:2201)/20-10,diam_averages(i,:)) 
end

%% Plot baseline and SNR

load('epidural_LU_30_filtered_30-300')
load('epidural_LU_30_filtered_trials')
% load('subdural1_LU_30_filtered_30-300')
% load('subdural1_LU_30_filtered_trials')
averages2=averages(:,16*20:36*20-1);

data_30all = data;
data_30all(indices,:)=[];
data_diams = data(indices,:);
baseline = zeros(length(indices),400);
lb = 400;
for ch=1:length(indices)
    baseline_trials = zeros(length(locs),lb);
    for l=1:length(locs)
        baseline_trials(l,:) = data_diams(ch,locs(l)-300-lb+1:locs(l)-300); % Take 20 ms segment of baseline before start of each stim
    end
    baseline(ch,:)=mean(baseline_trials);
end
goodchs2 = goodchs;
goodchs2(indices)=[];
baseline_30all = zeros(length(goodchs2),400);
lb = 400;
for ch=1:length(goodchs2)
    baseline_trials = zeros(length(locs),lb);
    for l=1:length(locs)
        baseline_trials(l,:) = data_30all(ch,locs(l)-300-lb+1:locs(l)-300); % Take 20 ms segment of baseline before start of each stim
    end
    baseline_30all(ch,:)=mean(baseline_trials);
end

figure()
title('Standard deviation of baseline')
hold on
plot(diams(:,3),std(baseline,0,2),'.')
xlabel('Diameter (\mum)')
ylabel('Standard deviation (\muV)')
set(gca,'Fontsize',32)

av_baseline = mean(baseline,2);

for i=1:length(diams(:,1))
    d30=find(diams(:,3)==30);
    d60=find(diams(:,3)==60);
    d80=find(diams(:,3)==80);
    d110=find(diams(:,3)==110);
    d140=find(diams(:,3)==140);
    d170=find(diams(:,3)==170);
    d210=find(diams(:,3)==210);
    d240=find(diams(:,3)==240);
    d480=find(diams(:,3)==480);
end


baseline30 = mean(std(baseline(d30(1:6),:),0,2));
baseline60 = mean(std(baseline(d60(1:6),:),0,2));
baseline80 = mean(std(baseline(d80(1:6),:),0,2));
baseline110 = mean(std(baseline(d110(1:6),:),0,2));
baseline140 = mean(std(baseline(d140(1:6),:),0,2));
baseline170 = mean(std(baseline(d170(1:6),:),0,2));
baseline210 = mean(std(baseline(d210(1:5),:),0,2));
baseline240 = mean(std(baseline(d240(1:6),:),0,2));
baseline480 = mean(std(baseline(d480(1:6),:),0,2));
baseline_30all_std = mean(std(baseline_30all,0,2));

baseline30std = std(std(baseline(d30(1:6),:),0,2));
baseline60std = std(std(baseline(d60(1:6),:),0,2));
baseline80std = std(std(baseline(d80(1:6),:),0,2));
baseline110std = std(std(baseline(d110(1:6),:),0,2));
baseline140std = std(std(baseline(d140(1:6),:),0,2));
baseline170std = std(std(baseline(d170(1:6),:),0,2));
baseline210std = std(std(baseline(d210(1:5),:),0,2));
baseline240std = std(std(baseline(d240(1:6),:),0,2));
baseline480std = std(std(baseline(d480(1:6),:),0,2));
baseline_30all_stdstd = std(std(baseline_30all,0,2));
%%
set(groot,'defaultLineMarkerSize',45);
figure()
hold on
%plot(30,baseline30,'o')
e1=errorbar(60,baseline60,baseline60std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(80,baseline80,baseline80std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(110,baseline110,baseline110std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(140,baseline140,baseline140std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(170,baseline170,baseline170std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(210,baseline210,baseline210std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(240,baseline240,baseline240std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(480,baseline480,baseline480std,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
e1=errorbar(30,baseline_30all_std,baseline_30all_stdstd,'vertical','ko','MarkerFaceColor','k')
e1.LineWidth = 4;
e1.MarkerSize = m;
xlabel('Diameter (\mum)')
ylabel('Mean baseline (\muV)')
set(gca,'Fontsize',42)
set(gca,'linewidth',4)
ylim([0.08 0.6])
set(gcf,'Position',[680   479   691   499])

% Perform t-test
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d60(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d80(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d110(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d140(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d170(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d210(1:5),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d240(1:6),:),0,2))
[h,p] = ttest2(std(baseline_30all,0,2),std(baseline(d480(1:6),:),0,2))

%% Fix 30 um baseline
l=[];
figure()
hold on
for i=1:374
plot(baseline_30all(i,:))
end
for i=1:374
    if find(any(baseline_30all(i,:)<-2))
    l = [l i]
    end
end

baseline_30all(l,:) = [];
%% Power spectral density
addpath '/Users/samantha/Downloads/electrodeimpedancemappingcode_SC';
load('diameters_neuromon_Sam.mat')
elecCoords = diameters_neuromon_Sam(:,1:2);
diams = diameters_neuromon_Sam(:,3);
elecX = elecCoords(:,1);
elecY = elecCoords(:,2);
% Create longer average file: 300 ms
indexMap = [];
for i=1:length(goodchs)
   if elecY(goodchs(i))<1500 
       indexMap = [indexMap;i];
   end
end
averages = [];
for ch=1:length(goodchs)
    ch
    data_array = [];
     %figure()
     %hold on
for i=1:length(locs)
   data_array(i,:) = data(ch,locs(i):locs(i)+300*20);
    %plot((1:2201)./20,data_array(i,:)-mean(data_array(i,:)))
end
averages_longer(ch,:) = mean(data_array,1); 
end
averages2_longer=averages_longer;
averages2_longer(indexMap,:) =[];
%%
rng default
Fs = 20000;
t = 0:1/Fs:1-1/Fs;

[psdx30,freq30]=psdx_calc(mean(averages2_longer(d30(1:3),16*20:56*20)));
[psdx60,freq60]=psdx_calc(mean(averages2_longer(d60(1:3),16*20:56*20)));
[psdx80,freq80]=psdx_calc(mean(averages2_longer(d80(1:3),16*20:56*20)));
[psdx110,freq110]=psdx_calc(mean(averages2_longer(d110(1:3),16*20:56*20)));
[psdx140,freq140]=psdx_calc(mean(averages2_longer(d140(1:3),16*20:56*20)));
[psdx170,freq170]=psdx_calc(mean(averages2_longer(d170(1:3),16*20:56*20)));
[psdx210,freq210]=psdx_calc(mean(averages2_longer(d210(1:3),16*20:56*20)));
[psdx240,freq240]=psdx_calc(mean(averages2_longer(d240(1:3),16*20:56*20)));
[psdx480,freq480]=psdx_calc(mean(averages2_longer(d480(1:3),16*20:56*20)));

addpath '/Users/samantha/Downloads/github_repo'
cmap=colormap(brewermap(9,'rdylbu'));
figure()
hold on
plot(freq30,10*log10(psdx30),'.-','Color',cmap(1,:),'Linewidth',3)
plot(freq60,10*log10(psdx60),'.-','Color',cmap(2,:),'Linewidth',3)
plot(freq80,10*log10(psdx80),'.-','Color',cmap(3,:),'Linewidth',3)
plot(freq110,10*log10(psdx110),'.-','Color',cmap(4,:),'Linewidth',3)
plot(freq140,10*log10(psdx140),'.-','Color',cmap(5,:),'Linewidth',3)
plot(freq170,10*log10(psdx170),'.-','Color',cmap(6,:),'Linewidth',3)
plot(freq210,10*log10(psdx210),'.-','Color',cmap(7,:),'Linewidth',3)
plot(freq240,10*log10(psdx240),'.-','Color',cmap(8,:),'Linewidth',3)
plot(freq480,10*log10(psdx480),'.-','Color',cmap(9,:),'Linewidth',3)
grid off
xlim([0 1000])
%title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('30 \mum','60 \mum','80 \mum','110 \mum','140 \mum','170 \mum','210 \mum','240 \mum','480 \mum')
legend boxoff
set(gca,'Fontsize',32)
set(gca,'linewidth',4)
set(gcf,'Position',[680   479   691   499])