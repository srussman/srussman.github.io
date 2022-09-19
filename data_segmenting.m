%% Load all files and align to see
%clear all
%close all
dat = [];
t = [];
ch=10;

% Load files of interest
load('Subdural1_preresection_210315_115251_goodchannels')
dat = [dat amplifier_data(ch,:)];
t = [t t_amplifier];
load('Subdural1_preresection_210315_115402_goodchannels')
dat = [dat amplifier_data(ch,:)];
t = [t t_amplifier];
load('Subdural1_preresection_210315_115513_goodchannels')
dat = [dat amplifier_data(ch,:)];
t = [t t_amplifier];
load('Subdural1_preresection_210315_115624_goodchannels')
dat = [dat amplifier_data(ch,:)];
t = [t t_amplifier];
load('Subdural1_preresection_210315_115734_goodchannels')
dat = [dat amplifier_data(ch,:)];
t = [t t_amplifier];
load('Subdural1_preresection_210315_115844_goodchannels')
dat = [dat amplifier_data(ch,:)];
t = [t t_amplifier];

% Plot to see stimulus artifacts
figure()
hold on
plot(t,dat)
%plot(t,dat2)
%plot(t,dat3)

%% Segment data: Left upper 35 mA
load('Subdural1_preresection_210319_115140_goodchannels')
load('Subdural1_preresection_210319_115140_filtered_30-300')

figure()
hold on
plot(t_amplifier,amplifier_data(10,:))
plot(t_amplifier(3*20000:38*20000),filtered_data(10,3*20000:38*20000))

data = filtered_data(:,3*20000:38*20000);
data_orig = amplifier_data(:,3*20000:38*20000);
t = t_amplifier(3*20000:38*20000);

save('subdural1_RU_30_filtered_30-300','data','t','data_orig')

%% Left upper 18 mA
dat = [];
t = [];
f = [];
load('Subdural2_post-resection_210319_131558_goodchannels')
load('Subdural2_post-resection_210319_131558_filtered_30-300')
dat = [dat amplifier_data];
t = [t t_amplifier];
f = [f filtered_data];
load('Subdural2_post-resection_210319_131704_goodchannels')
load('Subdural2_post-resection_210319_131704_filtered_30-300')
dat = [dat amplifier_data];
t = [t t_amplifier];
f = [f filtered_data];

data = f(:,42*20000:74*20000);
data_orig = dat(:,42*20000:74*20000);
t = t(42*20000:74*20000);

save('subdural2_RU_25_filtered_30-300','data','t','data_orig','-v7.3')

%% Left upper 9 mA - time 00:34
clear all
close all
load('Subdural1_preresection_210315_115251_goodchannels')
load('Subdural1_preresection_210315_115251_filtered_30-300')

figure()
hold on
plot(t_amplifier,amplifier_data(10,:))
plot(t_amplifier(2.5*20000:22*20000),filtered_data(10,2.5*20000:22*20000))

data = filtered_data(:,2.5*20000:22*20000);
data_orig = amplifier_data(:,2.5*20000:22*20000);
t = t_amplifier(2.5*20000:22*20000);

save('subdural1_LU_75_filtered_70-170','data','t','data_orig')

%% Left upper 4 mA
figure()
hold on
plot(t_amplifier,amplifier_data(10,:))
plot(t_amplifier(23*20000:39*20000),filtered_data(10,23*20000:39*20000))

data = filtered_data(:,23*20000:39*20000);
data_orig = amplifier_data(:,23*20000:39*20000);
t = t_amplifier(23*20000:39*20000);

save('subdural1_LU_4_filtered_70-170','data','t','data_orig')

%% Right upper 30 mA
clear all
close all
load('Subdural1_preresection_210315_115402_goodchannels')
load('Subdural1_preresection_210315_115402_filtered_30-300')

figure()
hold on
plot(t_amplifier,amplifier_data(10,:))
plot(t_amplifier(7*20000:22*20000),filtered_data(10,7*20000:22*20000))

data = filtered_data(:,7*20000:22*20000);
data_orig = amplifier_data(:,7*20000:22*20000);
t = t_amplifier(7*20000:22*20000);

save('subdural1_RU_30_filtered_30-300','data','t','data_orig')

%% Right upper 15 mA

figure()
hold on
plot(t_amplifier,amplifier_data(40,:))
ylim([-400 400])
plot(t_amplifier(22*20000:38*20000),filtered_data(40,22*20000:38*20000))

data = filtered_data(:,22*20000:38*20000);
data_orig = amplifier_data(:,22*20000:38*20000);
t = t_amplifier(22*20000:38*20000);

save('subdural1_RU_15_filtered_30-300','data','t','data_orig')

%% Right upper 7 mA
clear all
close all
load('Subdural1_preresection_210315_115402_goodchannels')
load('Subdural1_preresection_210315_115402_filtered_30-300')

figure()
hold on
plot(t_amplifier,amplifier_data(10,:))
ylim([-400 400])
plot(t_amplifier(40*20000:55*20000),filtered_data(10,40*20000:55*20000))

data = filtered_data(:,40*20000:55*20000);
data_orig = amplifier_data(:,40*20000:55*20000);
t = t_amplifier(40*20000:55*20000);

save('subdural1_RU_75_filtered_30-300','data','t','data_orig')

%% Right upper 4 mA
dataa = filtered_data(:,55*20000:end);
data_origa = amplifier_data(:,55*20000:end);
ta = t_amplifier(55*20000:end);

load('Subdural1_preresection_210315_115513_goodchannels')
load('Subdural1_preresection_210315_115513_filtered_30-300')

figure()
hold on
plot(t_amplifier,amplifier_data(10,:))
ylim([-400 400])
plot(t_amplifier(1:9*20000),filtered_data(10,1:9*20000))

datab = filtered_data(:,1:9*20000);
data_origb = amplifier_data(:,1:9*20000);
tb = t_amplifier(1:9*20000);

data = [dataa datab];
data_orig = [data_origa data_origb];
t = [ta tb];

save('subdural1_RU_4_filtered_30-300','data','t','data_orig')
