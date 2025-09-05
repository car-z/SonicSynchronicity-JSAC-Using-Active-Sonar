clear
close all
%% -------- CONSTANTS ------------------------------
c = 343;                % Speed of sound (m/s)
B = 4000;               % Bandwidth of chirp (Hz)
T = 0.0213;             % Duration of chirp (s)
N = 1024;               % Samples per chirp
Fc = 18000;             % Chirp start frequency (Hz)
ref_mic = 6;            % mic being processed on UMA-8 Mic Aray

%% -------- TRANSMITTED SIGNAL ----------------------
[tx, ~] = audioread('chirp21.wav');     % transmitted chirp
tx_chirp = tx(1:N);                     % extracting one singular chirp
tx_chirp = tx_chirp / norm(tx_chirp);

%% --------- RECEIVED SIGNAL ------------------------
[rx, Fs] = audioread('hand gestures/by speaker/chirp21-stationaryHigh-uma8.wav');   
if Fs ~= 48000                          % sample must be recorded at 48 kHz
    error('Sample rate mismatch: recorded at %.0f Hz.', Fs);
end

rx = rx(:, 1:7);                        % removing empty 8th channel on UM8 Mic Array

% Bandpass filter to frequency bands of interest
bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, ...
    'CutoffFrequency1', Fc - 500, ...
    'CutoffFrequency2', Fc + B + 500, ...
    'SampleRate', Fs);
rx_filtered = filter(bpFilt, rx);
rx_filtered = rx_filtered ./ vecnorm(rx_filtered);  % normalize each mic channel

[~, num_channels] = size(rx_filtered);              % num_channels is number of microphone channels recorded

%% ----- DETECT CHIRP START LOCATIONS USING ONE CHANNEL ONLY -------
ref_channel = rx_filtered(:,ref_mic);      % mic channel to use to detect start location of each chirp

% doing convolution (instead of xcorr) because xcorr gives negative lags,
% which we do not care about (we only care about delayed signal)
corr = conv(ref_channel, flipud(tx_chirp));

% each peak in corr is the start of a chirp
% to count as a peak: value is at least of threshold amplitude and ~1 chirp apart from previous peak
threshold = mean(corr) + 2*std(corr);
[~, locs] = findpeaks(corr, 'MinPeakHeight', threshold, ...
                        'MinPeakDistance', round(N*0.8)); 
num_chirps = length(locs);

time_axis = (0:num_chirps - 1) * T; % convert x-axis from chirp number to time

% % --- Plot correlation with peaks ---
% figure;
% plot(corr);
% hold on;
% plot(locs, pks, 'rx', 'LineWidth', 2);
% title('Cross-correlation with Transmitted Chirp');
% xlabel('Sample Index');
% ylabel('Correlation');
% legend('Correlation', 'Detected Chirps');
% 
% % --- Plot each extracted chirp ---
% figure;
% for i = 1:num_chirps
%     idx = locs(i) - N + 1;  % adjust to get start of chirp
%     if idx < 1 || (idx+N-1) > length(rx_filtered)
%         continue;  % skip if index out of bounds
%     end
%     chirp_segment = rx_filtered(idx:idx+N-1);
% 
%     subplot(ceil(num_chirps/4), 4, i);  % adjust layout as needed
%     plot(chirp_segment);
%     title(sprintf('Chirp #%d at %d', i, idx));
%     xlabel('Sample');
%     ylabel('Amplitude');
% end

%% ----------------- SET UP RANGE BINS -----------------
max_distance = 1;                          % meters
f_max = (2 * B * max_distance) / (c * T);  % max beat frequency for 1m target
max_bin = floor(f_max * N / Fs);           % only up to max_bin number of FFT bins are relevant to 1m target

min_distance = 0.05;                        % meters
f_min = (2 * B * min_distance) / (c * T);   % min beat frequency for 0.05m target
min_bin = floor(f_min * N / Fs);            % only after min_bin number of FFT bins are relevant

dist_axis = (min_bin:max_bin-1) * Fs / N * c * T / (2 * B);  % maps frequency bins to distance in m

%% ------------ LOOP OVER CHANNELS & CHIRPS --------------

% Preallocate storage matrix for FFT results
% rows: frequency bins, columns: each chirp, 3rd dimension: each channel
range_map_move = zeros(max_bin-min_bin+1, num_chirps, num_channels); 

% for each mic channel
for j = 1:num_channels
    % for each detected chirp
    for k = 1:num_chirps
        % extract the chirp start and end locations
        start_idx = locs(k);
        end_idx = start_idx + N - 1;
        if end_idx > length(rx_filtered) % if chirp is cut off by end of signal
            continue;
        end
        % extract singular chirp
        rx_chirp = rx_filtered(start_idx:end_idx, j);
        % compute beat
        beat = rx_chirp .* tx_chirp;
        % find frequency of beat
        beat_fft = abs(fft(beat, N));
        % store the result, only for relevant bins
        range_map_move(:, k,j) = beat_fft(min_bin:max_bin);
    end
end

range_map_dB = mag2db(range_map_move); % getting decibel values of created range map


%% --------- DO SAME PROCESSING FOR STATIC SIGNAL -----------------------------
% read in .wav file for received signal in static environment
[sx, ~] = audioread('hand gestures/by speaker/chirp21-empty-uma8.wav');

sx = sx(:, 1:7);
bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, ...
    'CutoffFrequency1', Fc - 500, ...
    'CutoffFrequency2', Fc + B + 500, ...
    'SampleRate', Fs);
sx_filtered = filter(bpFilt, sx);
sx_filtered = sx_filtered ./ vecnorm(sx_filtered);
[~, sx_num_channels] = size(sx_filtered);              
sx_ref_channel = sx_filtered(:,ref_mic);
sx_corr = conv(sx_ref_channel, flipud(tx_chirp));
[~, sx_locs] = findpeaks(sx_corr, 'MinPeakHeight', threshold, ...
                        'MinPeakDistance', round(N*0.8)); 
sx_num_chirps = length(sx_locs);
range_map_static = zeros(max_bin-min_bin+1, sx_num_chirps, num_channels); 
for j = 1:num_channels
    for k = 1:sx_num_chirps
        start_idx = sx_locs(k);
        end_idx = start_idx + N - 1;
        if end_idx > length(sx_filtered)
            continue;
        end
        sx_chirp = sx_filtered(start_idx:end_idx, j);
        beat = sx_chirp .* tx_chirp;
        beat_fft = abs(fft(beat, N));
        range_map_static(:, k,j) = beat_fft(min_bin:max_bin);
    end
end

%% ---------- SUBTRACTING STATIC NOISE ------------------
background_noise = median(range_map_static, 2);  % average across all chirps in static signal
range_map_sub = range_map_move - background_noise;  % subtracting static noise
range_map_sub(range_map_sub < 0) = 0;  % clip negatives
range_map_sub_dB = mag2db(range_map_sub + eps); % convert to dB

%% range time plot (no static noise subtracted)
figure;
set(gcf, 'Position', [200, 200, 1200, 600]);
tiledlayout(2,4);

for i = 1:7
    nexttile;
    imagesc(time_axis, dist_axis, range_map_dB(:,:,i));
    axis xy;
    xlabel('Time (s)');
    ylabel('Distance (m)');
    title(sprintf('Range-Time Plot (Mic %d)', i-1));
    colorbar;
    colormap jet;
    ylim([0.05, 1]);
end

% graphing trajectory by following peaks in range map
nexttile;
traj = zeros(1,num_chirps);
for i = 1:num_chirps
    [~, peak_bin] = max(range_map_dB(:, i, ref_mic)); % find freq bin of max motion
    range = dist_axis(peak_bin); % find corresponding distance
    traj(i) = range;
end
smoothed_traj = smoothdata(traj, 'movmean', 10);
plot(time_axis, smoothed_traj);
xlabel('Time (s)');
ylabel('Distance (m)');
title('Trajectory');

%% ----- range-time plots (background noise subtracted) -----------
figure;
set(gcf, 'Position', [300, 300, 1200, 600]);
tiledlayout(2,4);

for i = 1:7
    nexttile;
    imagesc(time_axis, dist_axis, range_map_sub_dB(:,:,i));
    axis xy;
    xlabel('Time (s)');
    ylabel('Distance (m)');
    title(sprintf('Motion only (Mic %d)', i-1));
    colorbar;
    colormap jet;
    clim([-100 -20]);
    ylim([0.05, 1]);
end

% graphing trajectory (motion only)
nexttile;
motion = zeros(1,num_chirps);
for i = 1:num_chirps
    [~, peak_bin] = max(range_map_sub_dB(:, i, ref_mic));
    range = dist_axis(peak_bin);
    motion(i) = range;
end
smoothed_motion = smoothdata(motion, 'movmean', 10);
plot(time_axis, smoothed_motion);
xlabel('Time (s)');
ylabel('Distance (m)');
title('Trajectory');


%% ----------- AoA Estimation (3 mic) ----------------
% Mics must be in a straight line (left, center, right)

d = 0.045;  % spacing between adjacent microphones, assumed, in m

% Select three microphones: left (mic 2), center (mic 0), right (mic 5)
mic_indices = [3, 1, 6];

% Use detected chirp locations to extract beat signals at each chirp
% allocating storage for each chirp
aoa_estimates = zeros(1, num_chirps); 

% for each chirp
for i = 1:num_chirps
    idx = locs(i); % start index of chirp
    end_idx = idx + N - 1;
    if end_idx > size(rx_filtered, 1)
        continue;
    end
    beat_matrix = zeros(N, 3); % rows: samples in chirp, columns: channels
    % for each channel of the mics
    for m = 1:3
        rx_chirp = rx_filtered(idx:end_idx, mic_indices(m)); % extract singular chirp 
        beat = rx_chirp .* tx_chirp; % calculate beat chirp
        beat_matrix(:, m) = beat; % store
    end

    % FFT and find max bin (most energy)
    fft_matrix = fft(beat_matrix, N);
    [~, peak_bin] = max(mean(abs(fft_matrix).^2, 2)); % find freq bin of highest energy among all 3 channels
    snapshot = fft_matrix(peak_bin, :).'; % extract bin of highest energy

    % AoA estimation via phase difference
    phase_diff_1 = angle(snapshot(1)) - angle(snapshot(2)); % between left and center
    phase_diff_2 = angle(snapshot(3)) - angle(snapshot(2)); % between right and center

    tau_1 = phase_diff_1 / (2 * pi * (Fc + B/2)); % calculating time delay based on phase delay
    tau_2 = phase_diff_2 / (2 * pi * (Fc + B/2)); % center frequency approximates beat f

    % Estimate angle using arcsin (only works for linear mic config)
    sin_theta_1 = tau_1 * c / d;
    sin_theta_2 = tau_2 * c / d;

    % Clamp (asin can only accept inputs between [-1, 1]
    sin_theta_1 = min(max(sin_theta_1, -1), 1);
    sin_theta_2 = min(max(sin_theta_2, -1), 1);

    % average the two angles, then find angle in degrees
    aoa_estimates(i) = asind((sin_theta_1 + sin_theta_2)/2);
end

% Plot AoA over time
figure;
set(gcf, 'Position', [400, 400, 560, 425]);  % [left, bottom, width, height]
plot(time_axis, aoa_estimates, '-o');
xlabel('Time (s)');
ylabel('Estimated AoA (degrees)');
title('AoA Estimation using 3-Microphone Linear Array');
grid on;

%% ------------------- 3D Motion Plot ---------------------

% Preallocate vectors
x = zeros(1, num_chirps);
z = zeros(1, num_chirps);
range_vals = zeros(1, num_chirps);

% for every chirp
for i = 1:num_chirps
    if i > size(range_map_sub_dB, 2)
        continue;
    end

    % Find the range bin with the strongest motion energy at that chirp
    % THIS IS WHERE THERE IS MOTION!
    [~, peak_bin] = max(range_map_sub_dB(:, i, ref_mic));

    % Find corresponding distance using dist_axis
    range = dist_axis(peak_bin);
    theta = aoa_estimates(i);  % angle estimate for this chirp in degrees

    % Convert polar to Cartesian
    x(i) = range * sind(theta);   % horizontal (side-to-side)
    z(i) = range * cosd(theta);   % vertical (height)
    range_vals(i) = range;        % for coloring
end

% Optional smoothing
x_smooth = smoothdata(x, 'movmean', 2);
z_smooth = smoothdata(z, 'movmean', 2);

% Plot 3D trajectory (time vs. position)
figure;
set(gcf, 'Position', [600, 400, 600, 450]);
cmap = jet(length(x));
hold on;

% Draw connecting lines
for i = 1:length(x)-1
    plot3([time_axis(i), time_axis(i+1)], ...
          [x_smooth(i), x_smooth(i+1)], ...
          [z_smooth(i), z_smooth(i+1)], ...
          'Color', cmap(i,:), 'LineWidth', 2);
end

% Overlay scatter
scatter3(time_axis, x_smooth, z_smooth, 36, range_vals, 'filled');

xlabel('Time (s)');
ylabel('X (Horizontal) [m]');
zlabel('Z (Height) [m]');
title('3D Motion Trajectory (Time vs. AoA vs. Range)');
ylim([-0.5 0.5])
colormap jet;
colorbar;
grid on;
view(45, 25);
