%% Clinical data analysis script
%% L UE C4-FZ
clear all
% Load data file
addpath '/Users/samantha/Desktop/Clinical_data'

% Import csv file set - this is the only thing to change
trace_data = csvread('L UEC4-FZ_trace_data.csv');
timestamp_labels = csvread('L UEC4-FZ_timestamp_labels.csv');
%messages = csvread('L UEC4-FZ_messages.csv',1,1,[1 1 83 1]);
trace_data_scalar = csvread('L UEC4-FZ_trace_data_scalar.csv');
timestamp_list = csvread('L UEC4-FZ_timestamp_list.csv');

% Filtering
Fs = 640/0.05;
time = (1:length(trace_data(1,:)))./Fs*1000; % in ms

% Filter data and plot
for ch=1:length(trace_data(:,1))
data_seg = trace_data(ch,:);
RawFilt = data_seg;
     for nf = [60,120,180,240,300,360,420,480]
        Wo = nf/(Fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        RawFilt = filtfilt(b,a,RawFilt')';
     end
     filtered_data1(ch,:)= RawFilt;
end

%Run bandpass filter
[b,a] = butter(4,[30 750]./(Fs/2));
filtered_data=filtfilt(b,a,filtered_data1')';

% Apply smoothing filter
upper_limit = Fs*0.0156; % SF5 = 200 Hz
[b,a] = butter(4,upper_limit./(Fs/2),'high');
filtered_data_smooth=filtfilt(b,a,filtered_data')';

%
figure()
hold on
for i=160 %:40:length(filtered_data(:,1))
plot(filtered_data(i,:))
plot(filtered_data_smooth(i,:))
end

% Converting to uV with trace data scalar
% Voltage units: 2*10^4*10^-10 in uV
% Trace Data Scalar: 3.6439453833736478E-10
filtered_data_uV = zeros(length(filtered_data_smooth(:,1)),length(filtered_data_smooth(1,:)));
for row = 1:length(trace_data_scalar)
   filtered_data_uV(row,:) = filtered_data_smooth(row,:).*trace_data_scalar(row)*10^6;
end

% Important time points - compare to time stamp labels
% R UE baseline
stamp = timestamp_labels(13);
[val,idx]=min(abs(timestamp_list-stamp));

% R
stamp2 = timestamp_labels(61);
[val,idx2]=min(abs(timestamp_list-stamp2));

% Difficulty with right SSEP
stamp3 = timestamp_labels(62);
[val,idx3]=min(abs(timestamp_list-stamp3));

% Important plots
figure()
hold on
plot(time,filtered_data_uV(idx,:))
plot(time,filtered_data_uV(idx2,:))
plot(time,filtered_data_uV(idx3,:))
%% L UE C3-C4
clear all
% Load data file
addpath '/Users/samantha/Desktop/Clinical_data'

% Import csv file set - this is the only thing to change
trace_data = csvread('L UEC3-C4_trace_data.csv');
timestamp_labels = csvread('L UEC3-C4_timestamp_labels.csv');
trace_data_scalar = csvread('L UEC3-C4_trace_data_scalar.csv');
timestamp_list = csvread('L UEC3-C4_timestamp_list.csv');
load('messages')

% Filtering
Fs = 640/0.05;
time = (1:length(trace_data(1,:)))./Fs*1000; % in ms

% Filter data and plot
for ch=1:length(trace_data(:,1))
data_seg = trace_data(ch,:);
RawFilt = data_seg;
     for nf = [60,120,180,240,300,360,420,480]
        Wo = nf/(Fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        RawFilt = filtfilt(b,a,RawFilt')';
     end
     filtered_data1(ch,:)= RawFilt;
end

%Run bandpass filter
[b,a] = butter(4,[30 750]./(Fs/2));
filtered_data=filtfilt(b,a,filtered_data1')';

% Apply smoothing filter
upper_limit = Fs*0.0156; % SF5 = 200 Hz
[b,a] = butter(4,upper_limit./(Fs/2),'high');
filtered_data_smooth=filtfilt(b,a,filtered_data')';

%
figure()
hold on
for i=160 %:40:length(filtered_data(:,1))
plot(filtered_data(i,:))
plot(filtered_data_smooth(i,:))
end

% Converting to uV with trace data scalar
% Voltage units: 2*10^4*10^-10 in uV
% Trace Data Scalar: 3.6439453833736478E-10
filtered_data_uV = zeros(length(filtered_data_smooth(:,1)),length(filtered_data_smooth(1,:)));
for row = 1:length(trace_data_scalar)
   filtered_data_uV(row,:) = filtered_data_smooth(row,:).*trace_data_scalar(row)*10^6;
end
%% L UE CV-FZ
clear all
% Load data file
addpath '/Users/samantha/Desktop/Clinical_data'

% Import csv file set - this is the only thing to change
trace_data = csvread('L UECV-FZ_trace_data.csv');
timestamp_labels = csvread('L UECV-FZ_timestamp_labels.csv');
trace_data_scalar = csvread('L UECV-FZ_trace_data_scalar.csv');
timestamp_list = csvread('L UECV-FZ_timestamp_list.csv');
load('messages')

% Filtering
Fs = 640/0.05;
time = (1:length(trace_data(1,:)))./Fs*1000; % in ms

% Filter data and plot
for ch=1:length(trace_data(:,1))
data_seg = trace_data(ch,:);
RawFilt = data_seg;
     for nf = [60,120,180,240,300,360,420,480]
        Wo = nf/(Fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        RawFilt = filtfilt(b,a,RawFilt')';
     end
     filtered_data1(ch,:)= RawFilt;
end

%Run bandpass filter
[b,a] = butter(4,[30 750]./(Fs/2));
filtered_data=filtfilt(b,a,filtered_data1')';

% Apply smoothing filter
upper_limit = Fs*0.0156; % SF5 = 200 Hz
[b,a] = butter(4,upper_limit./(Fs/2),'high');
filtered_data_smooth=filtfilt(b,a,filtered_data')';

%
figure()
hold on
for i=160 %:40:length(filtered_data(:,1))
plot(filtered_data(i,:))
plot(filtered_data_smooth(i,:))
end

% Converting to uV with trace data scalar
% Voltage units: 2*10^4*10^-10 in uV
% Trace Data Scalar: 3.6439453833736478E-10
filtered_data_uV = zeros(length(filtered_data_smooth(:,1)),length(filtered_data_smooth(1,:)));
for row = 1:length(trace_data_scalar)
   filtered_data_uV(row,:) = filtered_data_smooth(row,:).*trace_data_scalar(row)*10^6;
end
%% R UE C3-C4
clear all
% Load data file
addpath '/Users/samantha/Desktop/Clinical_data'

% Import csv file set - this is the only thing to change
trace_data = csvread('R UEC3-C4_trace_data.csv');
timestamp_labels = csvread('R UEC3-C4_timestamp_labels.csv');
trace_data_scalar = csvread('R UEC3-C4_trace_data_scalar.csv');
timestamp_list = csvread('R UEC3-C4_timestamp_list.csv');
load('messages')

% Filtering
Fs = 640/0.05;
time = (1:length(trace_data(1,:)))./Fs*1000; % in ms

% Filter data and plot
for ch=1:length(trace_data(:,1))
data_seg = trace_data(ch,:);
RawFilt = data_seg;
     for nf = [60,120,180,240,300,360,420,480]
        Wo = nf/(Fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        RawFilt = filtfilt(b,a,RawFilt')';
     end
     filtered_data1(ch,:)= RawFilt;
end

%Run bandpass filter
[b,a] = butter(4,[30 750]./(Fs/2));
filtered_data=filtfilt(b,a,filtered_data1')';

% Apply smoothing filter
upper_limit = Fs*0.0156; % SF5 = 200 Hz
[b,a] = butter(4,upper_limit./(Fs/2),'high');
filtered_data_smooth=filtfilt(b,a,filtered_data')';

%
figure()
hold on
for i=160 %:40:length(filtered_data(:,1))
plot(filtered_data(i,:))
plot(filtered_data_smooth(i,:))
end

% Converting to uV with trace data scalar
% Voltage units: 2*10^4*10^-10 in uV
% Trace Data Scalar: 3.6439453833736478E-10
filtered_data_uV = zeros(length(filtered_data_smooth(:,1)),length(filtered_data_smooth(1,:)));
for row = 1:length(trace_data_scalar)
   filtered_data_uV(row,:) = filtered_data_smooth(row,:).*trace_data_scalar(row)*10^6;
end
%% R UE C3-FZ
clear all
% Load data file
addpath '/Users/samantha/Desktop/Clinical_data'

% Import csv file set - this is the only thing to change
trace_data = csvread('R UEC3-FZ_trace_data.csv');
timestamp_labels = csvread('R UEC3-FZ_timestamp_labels.csv');
trace_data_scalar = csvread('R UEC3-FZ_trace_data_scalar.csv');
timestamp_list = csvread('R UEC3-FZ_timestamp_list.csv');
load('messages')

% Filtering
Fs = 640/0.05;
time = (1:length(trace_data(1,:)))./Fs*1000; % in ms

% Filter data and plot
for ch=1:length(trace_data(:,1))
data_seg = trace_data(ch,:);
RawFilt = data_seg;
     for nf = [60,120,180,240,300,360,420,480]
        Wo = nf/(Fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        RawFilt = filtfilt(b,a,RawFilt')';
     end
     filtered_data1(ch,:)= RawFilt;
end

%Run bandpass filter
[b,a] = butter(4,[30 750]./(Fs/2));
filtered_data=filtfilt(b,a,filtered_data1')';

% Apply smoothing filter
upper_limit = Fs*0.0156; % SF5 = 200 Hz
[b,a] = butter(4,upper_limit./(Fs/2),'high');
filtered_data_smooth=filtfilt(b,a,filtered_data')';

%
figure()
hold on
for i=160 %:40:length(filtered_data(:,1))
plot(filtered_data(i,:))
plot(filtered_data_smooth(i,:))
end

% Converting to uV with trace data scalar
% Voltage units: 2*10^4*10^-10 in uV
% Trace Data Scalar: 3.6439453833736478E-10
filtered_data_uV = zeros(length(filtered_data_smooth(:,1)),length(filtered_data_smooth(1,:)));
for row = 1:length(trace_data_scalar)
   filtered_data_uV(row,:) = filtered_data_smooth(row,:).*trace_data_scalar(row)*10^6;
end
%% R UE CV-FZ
clear all
% Load data file
addpath '/Users/samantha/Desktop/Clinical_data'

% Import csv file set - this is the only thing to change
trace_data = csvread('R UECV-FZ_trace_data.csv');
timestamp_labels = csvread('R UECV-FZ_timestamp_labels.csv');
trace_data_scalar = csvread('R UECV-FZ_trace_data_scalar.csv');
timestamp_list = csvread('R UECV-FZ_timestamp_list.csv');
load('messages')

% Filtering
Fs = 640/0.05;
time = (1:length(trace_data(1,:)))./Fs*1000; % in ms

% Filter data and plot
for ch=1:length(trace_data(:,1))
data_seg = trace_data(ch,:);
RawFilt = data_seg;
     for nf = [60,120,180,240,300,360,420,480]
        Wo = nf/(Fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        RawFilt = filtfilt(b,a,RawFilt')';
     end
     filtered_data1(ch,:)= RawFilt;
end

%Run bandpass filter
[b,a] = butter(4,[30 750]./(Fs/2));
filtered_data=filtfilt(b,a,filtered_data1')';

% Apply smoothing filter
upper_limit = Fs*0.0156; % SF5 = 200 Hz
[b,a] = butter(4,upper_limit./(Fs/2),'high');
filtered_data_smooth=filtfilt(b,a,filtered_data')';

%
figure()
hold on
for i=160 %:40:length(filtered_data(:,1))
plot(filtered_data(i,:))
plot(filtered_data_smooth(i,:))
end

% Converting to uV with trace data scalar
% Voltage units: 2*10^4*10^-10 in uV
% Trace Data Scalar: 3.6439453833736478E-10
filtered_data_uV = zeros(length(filtered_data_smooth(:,1)),length(filtered_data_smooth(1,:)));
for row = 1:length(trace_data_scalar)
   filtered_data_uV(row,:) = filtered_data_smooth(row,:).*trace_data_scalar(row)*10^6;
end
%% Important part
% Load all left side data files and plot the sampe timepoints:
clear all
load('LUEC4-FZ_data')

stamp = timestamp_labels(16);
[val,idx]=min(abs(timestamp_list-stamp));

% R
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));

stamp3 = timestamp_labels(38);
[val,idx3]=min(abs(timestamp_list-stamp3));

stamp4 = timestamp_labels(50);
[val,idx4]=min(abs(timestamp_list-stamp4));

% Difficulty with right SSEP
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));

% Important plots
figure()
hold on
plot(time,filtered_data_uV(idx,:),'Linewidth',2)
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)
plot(time,filtered_data_uV(idx3,:),'Linewidth',2)
plot(time,filtered_data_uV(idx4,:),'Linewidth',2)
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)
legend('baseline','after exposure, pre-laminectomy','surgery paused for pathology','poor RUE SSEP','difficulty to read SSEP')
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title('L UE C4-FZ')
set(gca,'Fontsize',14)

clear all
load('LUEC3-C4_data')

stamp = timestamp_labels(16);
[val,idx]=min(abs(timestamp_list-stamp));

% R
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));

stamp3 = timestamp_labels(38);
[val,idx3]=min(abs(timestamp_list-stamp3));

stamp4 = timestamp_labels(50);
[val,idx4]=min(abs(timestamp_list-stamp4));

% Difficulty with right SSEP
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));

% Important plots
figure()
hold on
plot(time,filtered_data_uV(idx,:),'Linewidth',2)
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)
plot(time,filtered_data_uV(idx3,:),'Linewidth',2)
plot(time,filtered_data_uV(idx4,:),'Linewidth',2)
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)
legend('baseline','after exposure, pre-laminectomy','surgery paused for pathology','poor RUE SSEP','difficulty to read SSEP')
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title('L UE C3-C4')

clear all
load('LUECV-FZ_data')

stamp = timestamp_labels(16);
[val,idx]=min(abs(timestamp_list-stamp));

% R
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));

stamp3 = timestamp_labels(38);
[val,idx3]=min(abs(timestamp_list-stamp3));

stamp4 = timestamp_labels(50);
[val,idx4]=min(abs(timestamp_list-stamp4));

% Difficulty with right SSEP
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));

% Important plots
figure()
hold on
plot(time,filtered_data_uV(idx,:),'Linewidth',2)
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)
plot(time,filtered_data_uV(idx3,:),'Linewidth',2)
plot(time,filtered_data_uV(idx4,:),'Linewidth',2)
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)
legend('baseline','after exposure, pre-laminectomy','surgery paused for pathology','poor RUE SSEP','difficulty to read SSEP')
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title('L UE CV-FZ')
set(gca,'Fontsize',14)
set(gca,'Fontsize',14)

% Right side
clear all
load('RUEC3-FZ_data')

stamp = timestamp_labels(13);
[val,idx]=min(abs(timestamp_list-stamp));

% R
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));

stamp3 = timestamp_labels(38);
[val,idx3]=min(abs(timestamp_list-stamp3));

stamp4 = timestamp_labels(50);
[val,idx4]=min(abs(timestamp_list-stamp4));

% Difficulty with right SSEP
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));

% Important plots
figure()
hold on
plot(time,filtered_data_uV(idx,:),'Linewidth',2)
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)
plot(time,filtered_data_uV(idx3,:),'Linewidth',2)
plot(time,filtered_data_uV(idx4,:),'Linewidth',2)
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)
legend('baseline','after exposure, pre-laminectomy','surgery paused for pathology','poor RUE SSEP','difficulty to read SSEP')
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title('R UE C3-FZ')
set(gca,'Fontsize',14)

clear all
load('RUEC3-C4_data')

stamp = timestamp_labels(13);
[val,idx]=min(abs(timestamp_list-stamp));

% R
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));

stamp3 = timestamp_labels(38);
[val,idx3]=min(abs(timestamp_list-stamp3));

stamp4 = timestamp_labels(50);
[val,idx4]=min(abs(timestamp_list-stamp4));

% Difficulty with right SSEP
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));

% Important plots
figure()
hold on
plot(time,filtered_data_uV(idx,:),'Linewidth',2)
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)
plot(time,filtered_data_uV(idx3,:),'Linewidth',2)
plot(time,filtered_data_uV(idx4,:),'Linewidth',2)
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)
legend('baseline','after exposure, pre-laminectomy','surgery paused for pathology','poor RUE SSEP','difficulty to read SSEP')
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title('R UE C3-C4')

clear all
load('RUECV-FZ_data')

stamp = timestamp_labels(13);
[val,idx]=min(abs(timestamp_list-stamp));

%
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));

stamp3 = timestamp_labels(38);
[val,idx3]=min(abs(timestamp_list-stamp3));

stamp4 = timestamp_labels(50);
[val,idx4]=min(abs(timestamp_list-stamp4));

% Difficulty with right SSEP
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));

% Important plots
figure()
hold on
plot(time,filtered_data_uV(idx,:),'Linewidth',2)
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)
plot(time,filtered_data_uV(idx3,:),'Linewidth',2)
plot(time,filtered_data_uV(idx4,:),'Linewidth',2)
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)
legend('baseline','after exposure, pre-laminectomy','surgery paused for pathology','poor RUE SSEP','difficulty to read SSEP')
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title('R UE CV-FZ')
set(gca,'Fontsize',14)
set(gca,'Fontsize',14)

%% Plots: left side
addpath '/Users/samantha/Downloads/new_filtered_30-300Hz'
addpath '/Users/samantha/Desktop/Clinical_data'

figure()
% Subdural 1
load('subdural1_LU_30_data_30-300')
load('subdural1_goodchs')
ch = find(goodchs==360);
t = (1:2201)./20-10;

subplot(2,2,1)
y=plot(t,averages(ch,:),'k','Linewidth',2)
xlim([0 30])
ylim([-5 5])
yticks([-5 -2.5 0 2.5 5])
yticklabels({'-5' '-2.5' '0' '2.5' '5'})
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',20)

% Subdural 3
load('subdural3_LU_30_averages_30-300')
load('subdural2_goodchs')
ch = find(goodchs==360);
t = (1:2201)./20-10;

subplot(2,2,3)
y= plot(t,averages(ch,:),'k','Linewidth',2)
xlim([0 30])
ylim([-5 5])
yticks([-5 -2.5 0 2.5 5])
yticklabels({'-5' '-2.5' '0' '2.5' '5'})
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',20)

% Clinical pre-resection
subplot(2,2,2)
hold on

load('LUEC4-FZ_data')
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)

load('LUEC3-C4_data')
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)

load('LUECV-FZ_data')
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)

legend('C4-FZ','C3-C4','CV-FZ')
xlim([0 30])
ylim([-0.4 0.4])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',20)

% Clinical post-resection
subplot(2,2,4)
hold on

load('LUEC4-FZ_data')
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)

load('LUEC3-C4_data')
stamp5 = timestamp_labels(62);
[~,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)

load('LUECV-FZ_data')
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)

legend('C4-FZ','C3-C4','CV-FZ')
xlim([0 30])
ylim([-0.4 0.4])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',34)

set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)

%% Plots: right side
addpath '/Users/samantha/Downloads/new_filtered_30-300Hz'
addpath '/Users/samantha/Desktop/Clinical_data'

figure()
% Subdural 1
load('subdural1_RU_30_data_30-300')
load('subdural1_goodchs')
ch = find(goodchs==533);
t = (1:2201)./20-10;

subplot(2,2,1)
y=plot(t,averages(ch,:),'k','Linewidth',2)
xlim([0 30])
ylim([-5 5])
yticks([-5 -2.5 0 2.5 5])
yticklabels({'-5' '-2.5' '0' '2.5' '5'})
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',20)

% Subdural 3
load('subdural3_RU_30_averages_30-300')
load('subdural2_goodchs')
ch = find(goodchs==533);
t = (1:2201)./20-10;

subplot(2,2,3)
y= plot(t,averages(ch,:),'k','Linewidth',2)
xlim([0 30])
ylim([-5 5])
yticks([-5 -2.5 0 2.5 5])
yticklabels({'-5' '-2.5' '0' '2.5' '5'})
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',20)

% Clinical pre-resection
subplot(2,2,2)
hold on

load('RUEC3-FZ_data')
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)

load('RUEC3-C4_data')
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)

load('RUECV-FZ_data')
stamp2 = timestamp_labels(28);
[val,idx2]=min(abs(timestamp_list-stamp2));
plot(time,filtered_data_uV(idx2,:),'Linewidth',2)

legend('C3-FZ','C3-C4','CV-FZ')
xlim([0 30])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',20)

% Clinical post-resection
subplot(2,2,4)
hold on

load('RUEC3-FZ_data')
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)

load('RUEC3-C4_data')
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)

load('RUECV-FZ_data')
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'Linewidth',2)

legend('C3-FZ','C3-C4','CV-FZ')
xlim([0 30])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'Fontsize',14)

set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)

%% Post-resection
close all
l1= 14;
l2= 22;
figure()
load('RUEC3-FZ_data')
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
set(gca,'visible','off')
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'RUEC3-FZ_data.png')

load('RUEC3-C4_data')
figure()
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'RUEC3-C4_data.png')

load('RUECV-FZ_data')
figure()
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
set(gca,'visible','off')
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'RUECV-FZ_data.png')

figure()
load('LUEC4-FZ_data')
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
set(gca,'visible','off')
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'LUEC4-FZ_data.png')

load('LUEC3-C4_data')
figure()
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
set(gca,'visible','off')
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'LUEC3-C4_data.png')

load('LUECV-FZ_data')
figure()
stamp5 = timestamp_labels(62);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
set(gca,'visible','off')
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'LUECV-FZ_data.png')

%% Pre-resection
close all
l1= 14;
l2= 22;
figure()
load('RUEC3-FZ_data')
stamp5 = timestamp_labels(28);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)
set(gca,'visible','off')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'RUEC3-FZ_data_pre.png')

load('RUEC3-C4_data')
figure()
stamp5 = timestamp_labels(28);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.52 1.48])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)
set(gca,'visible','off')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'RUEC3-C4_data_pre.png')

load('RUECV-FZ_data')
figure()
stamp5 = timestamp_labels(28);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)
set(gca,'visible','off')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'RUECV-FZ_data_pre.png')

figure()
load('LUEC4-FZ_data')
stamp5 = timestamp_labels(28);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)
set(gca,'visible','off')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'LUEC4-FZ_data_pre.png')

load('LUEC3-C4_data')
figure()
stamp5 = timestamp_labels(28);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)
set(gca,'visible','off')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'LUEC3-C4_data_pre.png')

load('LUECV-FZ_data')
figure()
stamp5 = timestamp_labels(28);
[val,idx5]=min(abs(timestamp_list-stamp5));
plot(time,filtered_data_uV(idx5,:),'k','Linewidth',25)
xlim([l1 l2])
ylim([-1.5 1.5])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gcf,'Position',[465   344   775   650])
set(gca,'Fontsize',20)
set(gca,'visible','off')
%set(gcf,'Position',[465   344   1200   650])
set(gca,'Fontsize',26)
size_ = get(gcf,'Position')
size = 10 * size_(end-1:end)
res = 600;
set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
saveas(gcf,'LUECV-FZ_data_pre.png')
