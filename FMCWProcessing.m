function [] = FMCWProcessing(inputArg)

if ischar(inputArg) || isstring(inputArg)
    [rx, Fs] = audioread(inputArg);   
    if Fs ~= 48000
        error('Sample rate mismatch: recorded at %.0f Hz.', Fs);
    end
elseif isnumeric(inputArg)
    rx = inputArg;
    Fs = 48000;
else 
    error("invalid input. Must be a filename (string) or a audio matrix");
end

%% -------- CONSTANTS ------------------------------
c = 343;                % Speed of sound (m/s)
B = 4000;               % Bandwidth of chirp (Hz)
T = 0.0213;             % Duration of chirp (s)
N = 1024;               % Samples per chirp
Fc = 18000;             % Chirp start frequency (Hz)
ref_mic = 6;

%% -------- TRANSMITTED SIGNAL ----------------------
[tx, ~] = audioread('chirp21.wav');
tx_chirp = tx(1:N);
tx_chirp = tx_chirp / norm(tx_chirp);

%% --------- RECEIVED SIGNAL ------------------------

rx = rx(:, 1:7); % removing empty 8th channel

% Bandpass filter for isolation
bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, ...
    'CutoffFrequency1', Fc - 500, ...
    'CutoffFrequency2', Fc + B + 500, ...
    'SampleRate', Fs);
rx_filtered = filter(bpFilt, rx);
rx_filtered = rx_filtered ./ vecnorm(rx_filtered);  % normalize each channel

[range_map_movement, time_axis, dist_axis, locs] = range_map(rx_filtered, B, N, T, tx_chirp);

num_chirps = length(locs);

range_map_dB = mag2db(range_map_movement);

%% --------- STATIC SIGNAL -----------------------------
[sx, ~] = audioread('hand gestures/by speaker/chirp21-empty-uma8.wav');
sx = sx(:, 1:7); % removing empty 8th channel
% Bandpass filter for isolation
bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, ...
    'CutoffFrequency1', Fc - 500, ...
    'CutoffFrequency2', Fc + B + 500, ...
    'SampleRate', Fs);
sx_filtered = filter(bpFilt, sx);
sx_filtered = sx_filtered ./ vecnorm(sx_filtered);  % normalize each channel
[range_map_static, ~, ~, ~] = range_map(sx_filtered, B, N, T, tx_chirp);

%% ---------- SUBTRACTING STATIC NOISE ------------------
background_profile = median(range_map_static, 2);  % size: [num_bins, 1, num_channels], averaging across all chirps
range_map_sub = range_map_movement - background_profile;  % subtracting
range_map_sub(range_map_sub < 0) = 0;  % clip negatives
range_map_sub_dB = mag2db(range_map_sub + eps); % convert to dB
% motion_mask = range_map_sub_dB > -60;  % motion detection thresholding

%% ----------------- PLOTTING --------------------------
figure;
set(gcf, 'Position', [200, 200, 1200, 600]);  % [left, bottom, width, height]
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

% graphing trajectory
nexttile;
traj = zeros(1,num_chirps);
for i = 1:num_chirps
    [~, peak_bin] = max(range_map_dB(:, i, ref_mic)); % find freq bin of max motion
    range = dist_axis(peak_bin); % find corresponding distance
    traj(i) = range;
end
smoothed_traj = smoothdata(traj, 'movmean', 10);  % or 'gaussian', window=5
plot(time_axis, smoothed_traj);
xlabel('Time (s)');
ylabel('Distance (m)');
title('Trajectory');

figure;
set(gcf, 'Position', [300, 300, 1200, 600]);  % [left, bottom, width, height]
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
    [~, peak_bin] = max(range_map_sub_dB(:, i, ref_mic)); % find freq bin of max motion
    range = dist_axis(peak_bin); % find corresponding distance
    motion(i) = range;
end
smoothed_motion = smoothdata(motion, 'movmean', 10);  % or 'gaussian', window=5
plot(time_axis, smoothed_motion);
xlabel('Time (s)');
ylabel('Distance (m)');
title('Trajectory');


%% ----------- AoA Estimation (3 mic) ----------------
% mics are in a straight line (left, center, right)

d = 0.045;  % 45 mm, spacing between adjacent microphones, assumed

% Select three microphones: left (mic 2), center (mic 0), right (mic 5)
mic_indices = [3, 1, 6];

% Use detected chirp locations to extract beat signals
aoa_estimates = zeros(1, num_chirps); % allocating storage for each chirp

% for each chirp
for i = 1:num_chirps
    idx = locs(i); % start index of chirp
    if idx + N - 1 > size(rx_filtered, 1)
        continue;
    end
    beat_matrix = zeros(N, 3); % rows: samples in chirp, columns: channels
    % for each channel of the mics
    for m = 1:3
        rx_chirp = rx_filtered(idx:idx+N-1, mic_indices(m)); % extract singular chirp 
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

end

