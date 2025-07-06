function [range_map,time_axis, dist_axis, locs] = range_map(rx_filtered, B, N, T, tx_chirp)
%RANGE_MAP Summary of this function goes here
%   Detailed explanation goes here

fs = 48000;
c = 343;

[~, num_channels] = size(rx_filtered);

%% ----- DETECT CHIRP LOCATIONS IN ONE CHANNEL ONLY -------
ref_channel = rx_filtered(:,6);

corr = conv(ref_channel, flipud(tx_chirp));  % equivalent to xcorr
threshold = mean(corr) + 2*std(corr);
% each peak in corr is start of a chirp
[~, locs] = findpeaks(corr, 'MinPeakHeight', threshold, ...
                        'MinPeakDistance', round(N*0.8));  % ~1 chirp apart
num_chirps = length(locs);

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
max_distance = 1;  % meters
f_max = (2 * B * max_distance) / (c * T);  % max beat frequency for 1m target
max_bin = floor(f_max * N / fs);  % only up to max_bin number of FFT bins are relevant to 1m target

min_distance = 0.05;
f_min = (2 * B * min_distance) / (c * T);  % min beat frequency for 0.25m target
min_bin = floor(f_min * N / fs);  % only beginning at min_bin number of FFT bins are relevant to 1m target

dist_axis = (min_bin:max_bin-1) * fs / N * c * T / (2 * B);  % maps frequency bins to distance in m

%% ------------ LOOP OVER CHANNELS & CHIRPS --------------

% Preallocate storage matrix for FFT results
range_map = zeros(max_bin-min_bin+1, num_chirps, num_channels); % rows - frequency (y-axis), columns - chirp (x-axis), 3D - each channel

% for each channel
for j = 1:num_channels
    % for each detected chirp
    for k = 1:num_chirps
        % extract the singular chirp
        start_idx = locs(k);
        if start_idx + N - 1 > length(rx_filtered)
            continue;
        end
        rx_chirp = rx_filtered(start_idx:start_idx + N - 1, j);
        % compute beat
        beat = rx_chirp .* tx_chirp;
        % find frequency of beat
        beat_fft = abs(fft(beat, N));
        % store the result in the appropriate chirp, only up to max_bin (1 m)
        range_map(:, k,j) = beat_fft(min_bin:max_bin);
    end
end

%% --- Convert to Range-Time Axes ---
time_axis = (0:num_chirps - 1) * T; % convert x-axis from chirp number to time

end

