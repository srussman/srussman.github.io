% This code plots the phase gradients and streamlines of the traveling
% SSEPs


% Author: srussman@ucsd.edu

%% Baseline - Spontaneous activities
clear all
close all
addpath('/wave/analysis')
addpath('/wave/plotting')

% Electrode mapping information (X, Y)
addpath '/codes_for_wave_front_analysis_UCSD/electrodeimpedancemappingcode_SC';
load('padCoords');
load('diameters_neuromon_Sam.mat');
elecCoords = diameters_neuromon_Sam(:,1:2);
diams = diameters_neuromon_Sam(:,3);
elecX = elecCoords(:,1)';
elecY = elecCoords(:,2)';
% k vector contains the channel numbers with impedance below 100kOhm @ 1kHz
load('epidural_goodchs');

% Load data
load('epidural_LU_30_data_30-300') % filtered 30-300 Hz, sampled at 20 kHz

goodchs2=goodchs;
indexMap = [];
for i=1:length(goodchs) % only consider electrodes in main array
    if elecY(goodchs(i))>1500 & elecY(goodchs(i)) <100000
        indexMap = [indexMap;i];
    end
end
goodchs2= goodchs(indexMap); % channels in main array
averages2 =averages(indexMap,:); % averages of main array (10 ms pre stim to 100 ms post)
k=goodchs2';
%% Interpolation
% Data length
d_length = int32(length(averages2(1,:)));

num_interp_x = 12;
num_interp_y = 31;

% Interpolation
interp_x_segments = num_interp_x; % Number of segments in x-axis
interp_y_segments = num_interp_y; % Number of segments in y-axis
% set up new X, Y, and Z coordinates for interpolated output
[xmin,xmax] = bounds(elecX);
[ymin,ymax] = bounds(elecY);
[xGrid,yGrid] = meshgrid(linspace(xmin,xmax,interp_x_segments),linspace(ymin,ymax,interp_y_segments));

zGrid = zeros(num_interp_x,num_interp_y, d_length);
zGrid_amp = zeros(num_interp_x,num_interp_y, d_length);
zGrid_amp_smooth = zeros(num_interp_x,num_interp_y, d_length);

for t_loop = 1:d_length
    Finter = scatteredInterpolant(elecX(k)',elecY(k)',averages2(:, t_loop));
    zGrid(:,:,t_loop) = Finter(xGrid,yGrid)'; % interpolated array
    display(sprintf('Interpolation in progress.. %.1f%% completed.', 100*t_loop/d_length));
end
%% Phase Gradient Vector field plot
addpath('wave')
addpath('wave\analysis')
addpath('wave\plotting')

% parameters
fs = 20*1000; % Sampling rate in Hz
time_ini = 1; % in ms
image_size = 16; % px
pixel_spacing = 1; % a.u.
direction = +1; % +1/-1

% generate data
xf = zGrid;
% z-score data
xf = zscore_independent( xf );
% form analytic signal
xph = analytic_signal( xf );
% calculate instantaneous frequency
[wt,signIF] = instantaneous_frequency( xph, fs );
% calculate phase gradient
[pm,pd,dx,dy] = phase_gradient_complex_multiplication( xph, pixel_spacing, signIF );
display(sprintf('Phase gradient calculation completed.'));

%% RMS video GIF plotting
fName = 'Wave_front'; % File name
v_max = 8; % Max voltage plot range
v_min = 0; % min voltage plot range
pixel_size = 150; % Scatter plot pixel size
frame_rate = 10; % GIF frame rate

t_interval = 1; % Time interval between the frame in 'ms'
t_window = [1:t_interval:800]; % Time window in 'ms', decrease because response only in first 20 ms
figure2 = figure('Position', [0, 0, 0+800, 0+800]);
set(gcf,'color','w');

filename = strcat(fName, '.gif');
for  i = 1:length(t_window)
    time_stamp = t_window(i); % time in 'ms'
    % Scatter plot
    scatter((-elecX(k)+1925+350)./350,(-elecY(k)+1.35971349e+04+400)./400, pixel_size, averages2(:,20*10+ time_stamp*20),'filled', 'o'); % directions flipped to account for electrode orientation, also skip first 10 ms of pre-stim data, also adjust for 1:13,1:31 axes?
    %colormap(flipud(hot));
    
    % Scatter plot: Settings
    caxis([v_min v_max]);
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    txtTitle = strcat('Epidural pre-resection phase gradient',{' '}, 'at', {' '},num2str(time_stamp),{' '}, 'ms' );
    ht = title(txtTitle,'interpreter', 'none');
    ylabel(h, 'Voltage (\muV)');
    pbaspect([1 1 1]);
    %set(gca,'YTick',0:4:31);
    %set(gca,'YTickLabel',((0:4:31)-16));
    %set(gca,'XTick',0:4:31);
    %set(gca,'XTickLabel',((0:4:31)-16));
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',1);
    %axis([-2500 2500 -12500 -1600])
    axis([0 12 0 31])
    box on;
    % Colorbar: setting
    set(h,'linewidth',1);
    set(h,'FontSize',10);
    set(h,'FontWeight','bold');
    set(ht,'FontWeight','bold');
    
    hold on;
    
    % Plot Vector Field
    ph = exp( 1i .* (pd(:,:,time_stamp)) )'; % put transverse here, I think dimensions make sense now?
    assert( ~isreal(ph), 'complex-valued input required, ph' );
    
    % Vector direction
    [XX,YY] = meshgrid( 1:size(ph,2), 1:size(ph,1) );
    M = real( exp( 1i * angle(ph) ) );
    N = imag( exp( 1i * angle(ph) ) );
    
    % Streamlines Red
%     [sx,sy] = meshgrid(1:31, 1:12); % Setup Steamline Origins, according to Youngbin's code: 1:31, 1:12, but axes don't fit
%     hlines_sens = streamline(stream2(XX,YY,M,N,sx,sy)) ;
%     set(hlines_sens,'LineWidth',2,'Color','r') % Red
    hold on;
    
    % quiver plot / Vector field
    quiver( XX, YY, M, N, 0.75, 'k',  'linewidth', 1); % axes don't match the amplitude plots
    
    drawnow
    
    % Capture each frame as an image
    frame = getframe(figure2);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    dtime = 10./frame_rate;
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',dtime);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',dtime);
    end
    hold off;
    
end
close(figure2)