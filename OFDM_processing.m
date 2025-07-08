clear 
close all

filename = 'test_upDown2.wav';

%% Constants
Fs = 48000; % sampling rate (Hz)
symbol_len = 128; % OFDM symbol length (samples)

cyclic_suffix_len = 20; % (samples)
separation_len = 400; % each pulse separated by 400 samples

pulse_len = symbol_len + cyclic_suffix_len + separation_len;

%% Transmitted Signal
[tx, ~] = audioread('ofdm_pulse_548.wav');
ref_pulse = tx(1:pulse_len);

% t = (0:length(ref_pulse)-1)/Fs*1000;  % in milliseconds
% figure(1);
% plot(t, ref_pulse);
% xlabel('Time (ms)');
% ylabel('Amplitude');
% title('One OFDM Pulse (with cyclic suffix)');
% grid on;

%% read and filter received signal
[rx, ~] = audioread(filename);

rx = rx(:,6); % only  using one channel for now

% Bandpass filter for isolation
bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, ...
    'CutoffFrequency1', 17500, ...
    'CutoffFrequency2', 20500, ...
    'SampleRate', Fs);
rx_filtered = filter(bpFilt, rx);
rx_filtered = rx_filtered ./ vecnorm(rx_filtered);  % normalize each channel

%% detect pulse starts (with one channel)
corr = conv(rx_filtered, flipud(ref_pulse));  % equivalent to xcorr
threshold = mean(corr) + 2*std(corr);
% each peak in corr is start of a chirp
[~, locs] = findpeaks(corr, 'MinPeakHeight', threshold, ...
                        'MinPeakDistance', round(548*0.8));  % ~1 pulse apart
locs = locs - length(ref_pulse) + 1;
num_pulses = length(locs);

%% Generate ECHO PROFILES for all reflections at mic (STEP ONE)
corr_all = [];
for k = 1:num_pulses
    start_index = locs(k);
    end_index = start_index + pulse_len - 1;
    if start_index <= 0 || end_index > length(rx_filtered)
        continue
    end
    one_pulse = rx_filtered(start_index:end_index);
    [corr_output, lags] = xcorr(one_pulse,ref_pulse);
    corr_output = abs(corr_output)';
    corr_all = [corr_all; corr_output];
end

% Convert lag to distance
c = 343;  % speed of sound in m/s
sec_to_sample = 1 / Fs;
distance = lags * sec_to_sample * c / 2 * 100;  % round trip to one-way, convert to cm

% valid_idx = lags >= 0;
% distance = distance(valid_idx);

% mean_corr = mean(corr_all, 1);
% mean_corr = mean_corr(valid_idx);
% 
% plot(distance, mean_corr);
% xlabel('Distance (cm)');
% ylabel('Mean Correlation Magnitude');
% title('Average Echo Profile');
% grid on;

%% Identify the Echo corresponding to moving object (STEP TWO)

% computing threshold for motion change
avg_corr = mean(corr_all(1:3,:), 1); % using average of first three pulses
direct_amp = max(avg_corr); % max correlation should correspond to signal itself (not reflection)
motion_threshold = 0.15 * direct_amp;

[numEchoProfiles, length_corr] = size(corr_all);

mov_locs = [];
for k = 1:numEchoProfiles-1
    diff = corr_all(k+1,:)-corr_all(k,:);
    for j = 1:length_corr
        if diff(j) > motion_threshold
            mov_locs = [mov_locs, j];
            break;
        end
    end
end

figure;
moving_distances = distance(mov_locs);
plot(1:length(moving_distances), moving_distances, 'o-');
xlabel('Pulse Index');
ylabel('Moving Object Distance (cm)');
title('Motion Over Time');

