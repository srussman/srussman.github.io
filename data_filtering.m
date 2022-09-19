
%% Select raw data files
clear all
load('subdural1_goodchs') % load goodchs only, generated from impedancemappingcode
fileList = dir('subdural1*.mat');

%% Select only good channels
%cd 'epidural_goodchs'
for i=1:length(fileList) % run through entire fileList of raw data files
filename=fileList(i).name
file = load(filename);
amplifier_data = file.amplifier_data(goodchs,:);
board_adc_data = file.board_adc_data;
t_amplifier = file.t_amplifier;
dotLocations = find(filename == '.');
name = [filename(1:dotLocations(1)-1) '_goodchannels'];
save(name, 'amplifier_data', 'board_adc_data', 't_amplifier', 'goodchs', '-v7.3');
clear amplifier_data board_adc t_amplifier
end
%% Filter
fileList = dir('subdural*goodchannels.mat');

for i=1:length(fileList)
filename=fileList(i).name
file = load(filename);
amplifier_data=file.amplifier_data;
board_adc_data = file.board_adc_data;
t_amplifier = file.t_amplifier;
nfft = [1:1:250];
filtered_data = zeros(length(goodchs),length(amplifier_data(1,:)));

for ch=1:length(goodchs)
raw_data = amplifier_data(ch,:);
IntanFS = 20000;
RawFilt = raw_data;
     for nf = [60,120,180,240,300,360,420,480] % remove 60 Hz noise
        Wo = nf/(IntanFS/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        RawFilt = filtfilt(b,a,RawFilt')';
     end
     filtered_data(ch,:)= RawFilt;
end

IntanFS = 20000;

[b,a] = butter(4,[30 300]./(IntanFS/2)); % filter 30-300 Hz

filtered_data2=filtfilt(b,a,filtered_data')';
filtered_data = filtered_data2;

dotLocations = find(filename == '.');
name = [ filename(1:dotLocations(1)-13) 'filtered_30-300']
save(name, 'filtered_data', 'board_adc_data', 't_amplifier', 'goodchs', '-v7.3');
end

