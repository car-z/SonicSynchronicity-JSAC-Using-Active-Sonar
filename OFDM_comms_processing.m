clear 
close all

filename = 'received OFDM/prefix_comms6_updown3_548.wav';

%% Constants
Fs = 48000; % sampling rate (Hz)
symbol_len = 128; % OFDM symbol length (samples)

cyclic_len = 20; % (samples)
separation_len = 400; % each pulse separated by 400 samples

pulse_len = symbol_len + cyclic_len + separation_len;

% Convert lag to distance
    c = 343;  % speed of sound in m/s
    sec_to_sample = 1 / Fs;

subcarrier_width = Fs / (symbol_len); % = 375 Hz
f = (0:symbol_len-1) * subcarrier_width; % vector spanning 0-48 kHz in 128 subcarriers

% identify indices in f corresponding to 18-20 kHz
lower_k = ceil(18000 / subcarrier_width);
upper_k = floor(20000 / subcarrier_width);
active_bins = lower_k+1 : upper_k + 1; % MATLAB uses 1 based indexing

%% Transmitted Signal
[tx, ~] = audioread('OFDM signals/ofdm_pulse_101101_548_prefix.wav');
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

rx = rx(:,1:7); % removing empty channel

% Bandpass filter for isolation
bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, ...
    'CutoffFrequency1', 17500, ...
    'CutoffFrequency2', 20500, ...
    'SampleRate', Fs);
rx_filtered = filter(bpFilt, rx);
rx_filtered = rx_filtered ./ vecnorm(rx_filtered);  % normalize each channel

%% 
figure(1);
set(gcf, 'Position', [100, 100, 1200, 600]);
tiledlayout(2,4);

figure(2);
set(gcf, 'Position', [100, 100, 1200, 600]);
tiledlayout(2,4);

for i = 1:7
    %% detect pulse starts (with one channel)
    rx_one = rx_filtered(:,i);

    corr = conv(rx_one, flipud(ref_pulse));  % equivalent to xcorr
    
   
    [~, firstPeak] = findpeaks(corr, 'NPeaks', 1);
    num_pulses = length(rx_one) / pulse_len;

    locs = firstPeak + (0:num_pulses - 1) * pulse_len;
    % num_pulses = length(locs);
    
    % corr = conv(rx_one, flipud(ref_pulse));  % equivalent to xcorr
    % threshold = mean(corr) + 2*std(corr);
    % % each peak in corr is start of a chirp
    % [~, locs] = findpeaks(corr, 'MinPeakHeight', threshold, ...
    %                         'MinPeakDistance', round(pulse_len*0.8));  % ~1 pulse apart
    % locs = locs - length(ref_pulse) + 1;
    % num_pulses = length(locs);
    % 
    % figure(3);
    % plot(rx_one);
    % hold on;
    % plot(locs,rx_one(locs),'ro', 'MarkerSize',6,'LineWidth',1.2);
    % xlabel('sample');
    % xlim([min(locs)-500, max(locs)+500]);  % show ~500 samples around the pulse region
    % ylabel('amplitude');
    % title(sprintf('Channel %d - raw detected pulse starts',i));
    % legend('filtered signal','pulse starts');
    % grid on;
    % hold off;

    %% Generate ECHO PROFILES for all reflections at mic (STEP ONE)
    corr_all = []; % holds echo profiles for all pulses, where each row is a pulse
    for k = 1:num_pulses
        start_index = locs(k);
        end_index = start_index + pulse_len - 1;
        if start_index <= 0 || end_index > length(rx_one)
            continue
        end
        one_pulse = rx_one(start_index:end_index);
        [corr_output, lags] = xcorr(one_pulse,ref_pulse);
        corr_output = abs(corr_output)';

        valid_lag_idx = lags >= 0;
        lags = lags(valid_lag_idx); % only keeping corr_output and lags lags where lags >= 0
        corr_output = corr_output(valid_lag_idx);

        corr_all = [corr_all; corr_output];
    end
    
    distance = lags * sec_to_sample * c / 2 * 100;

    %% Identify the Echo corresponding to moving object (STEP TWO)

    % computing threshold for motion change
    avg_corr = mean(corr_all(1:3,:), 1); % using average of first three pulses
    direct_amp = max(avg_corr); % max correlation should correspond to signal itself (not reflection)
    
    motion_threshold = 0.2 * direct_amp;
    
    [numEchoProfiles, length_corr] = size(corr_all);
    
    mov_locs = zeros(1,numEchoProfiles); % index (within echo profile) of detected movement for all echo profiles
    for k = 1:numEchoProfiles-1
        % [~, mov_locs(k)] = max(abs(corr_all(k+1,:) - corr_all(k,:)));
        diff = corr_all(k+1,:)-corr_all(k,:);
        for j = 1:length_corr
            if abs(diff(j)) > motion_threshold
                mov_locs(k) = j;
                break;
            end
        end       
    end

    % Find indices where motion was detected
    valid_idx = mov_locs > 0;
    
    % Get pulse indices and corresponding distances
    valid_pulse_indices = find(valid_idx);  % x-axis
    moving_distances = distance(mov_locs(valid_idx));  % y-axis

    % Plot
    figure(1);
    nexttile;
    plot(valid_pulse_indices, moving_distances, 'o-');
    xlabel('Pulse Index');
    ylabel('Moving Object Distance (cm)');
    title(sprintf('Motion vs Time - Channel %d', i));
    grid on;

    %% STEP THREE or something
    ref_symbol = ref_pulse(cyclic_len+1:cyclic_len+symbol_len); % transmitted symbol alone (no prefix)    
    ref_spectrum = fft(ref_symbol,symbol_len)'; % frequency spectrum of transmitted signal

    fine_delays = zeros(1,length(mov_locs(valid_idx)));
    index = 1;

    for k = 1:numEchoProfiles

        comms_start = locs(k) + cyclic_len + 1;
        comms_end = comms_start + symbol_len - 1;
        if comms_end > length(rx_filtered)
            continue;
        end

        comms_spectrum = fft(rx_filtered(comms_start:comms_end),symbol_len);
        rx_bits = real(comms_spectrum(active_bins)) > 0;
        disp('Received bits:');
        disp(rx_bits);

        if mov_locs(k) == 0 % no motion occurred
            continue;
        end
        echo_start = locs(k) + mov_locs(k) - 1+cyclic_len; % index in rx_filtered = pulse start time + move start time - 1
        echo_end = echo_start + symbol_len - 1;
        if echo_end > length(rx_one)
            continue;
        end

        echo_segment = rx_one(echo_start:echo_end);
        echo_spectrum = fft(echo_segment,symbol_len)'; % frequency spectrum of the echo from moving object

        phase_diff = angle(echo_spectrum) - angle(ref_spectrum); % subtract known from received phase, giving you extra phase shift per subcarrier because of delay
        phase_diff = unwrap(phase_diff); % avoids 2pi jumps

        % phase shift as a function of frequency = -2pi*f*delay in time
        % therefore: delay in time 
        % = phase shift as function of frequency/(frequency*-2pi)
        % delay in time = (slope of phase shift vs f)/ -2pi
        
        f_bins = (0:symbol_len-1) * Fs/symbol_len; % calculating frequency
        p = polyfit(f_bins, phase_diff, 1); % finding linear approximation of phase shift as function of frequency
        delta_t = -p(1) / (2*pi);  % p(1) holds slope

        fine_delays(index) = delta_t;
        index = index+1;
    end

    distances_fine = fine_delays * c * 100 / 2;  % cm, round-trip halved

    pulse_times = locs(1:length(fine_delays)) / Fs;  % seconds

    figure(2);
    nexttile;
    plot(1:length(distances_fine), distances_fine, 'o-');
    xlabel('Pulse Index');
    ylabel('Fine-tuned Distance (cm)');
    title(sprintf('Precision Tracking - Channel %d',i));
    grid on;

end