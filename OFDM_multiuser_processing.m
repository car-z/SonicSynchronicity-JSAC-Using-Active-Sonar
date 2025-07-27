clear
close all

filename = ['test_new_setup_multihand.wav'];

motion_threshold_ratio = .15;

multihand = true;

%% Constants -----------------
Fc = 18000;
B = 2000;
Fs = 48000;
preamble_len = 64;
cyclic_len = 20; 
symbol_len = 128;
separation_len = 200;

pulse_len = preamble_len + cyclic_len + symbol_len + separation_len;

c = 343;
sec_to_sample = 1 / Fs;

subcarrier_width = Fs / (symbol_len);
lower_k = ceil(Fc / subcarrier_width);
upper_k = floor((Fc+B) / subcarrier_width);
active_bins = lower_k+1 : upper_k + 1;

%% Transmitted Signal ---------------
[tx, ~] = audioread('OFDM signals/ofdm_pulse_348_prefix_preamble.wav');
ref_pulse = tx(1:pulse_len);

%% Received Signal --------------
[rx, ~] = audioread(filename);

rx = rx(:,1:7); % removing empty channel

% Bandpass filter for isolation
bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, ...
    'CutoffFrequency1', Fc - 500, ...
    'CutoffFrequency2', Fc + B + 500, ...
    'SampleRate', Fs);
rx_filtered = filter(bpFilt, rx);
rx_filtered = rx_filtered ./ vecnorm(rx_filtered);

%% Set Up Figures ------------
figure(1);
set(gcf, 'Position', [100, 100, 1200, 600]);
tiledlayout(2,4);

figure(2);
set(gcf, 'Position', [100, 100, 1200, 600]);
tiledlayout(2,4);

for i = 1:7
    rx_one = rx_filtered(:,i);

    %% Detect Pulse Starts
    preamble = ref_pulse(1:preamble_len);

    % corr = conv(rx_one, flipud(preamble)); 
    % threshold = mean(corr) + 3 * std(corr);
    % [~, locs] = findpeaks(corr, ...
    %     'MinPeakHeight', threshold, ...
    %     'MinPeakDistance', round(pulse_len * 0.7));
    % locs = locs - length(preamble) + 1;
    % if locs(1) <= 0 
	%     locs(1) = [];
    % end
    % num_pulses = length(locs);

	% ALTERNATE METHOD
    corr = conv(rx_one, flipud(ref_pulse));
    [~, firstPeak] = findpeaks(corr, 'NPeaks', 1);
    num_pulses = floor(length(rx_one) / pulse_len);
    locs = firstPeak + (0:num_pulses - 1) * pulse_len;

    %% STEP ONE - Generate ECHO PROFILES for all reflections at mic ----------------
    echoProfiles = [];
    for k = 1:num_pulses
        start_index = locs(k);
        end_index = start_index + pulse_len - 1;
        if end_index > length(rx_one)
            continue
        end
        one_pulse = rx_one(start_index:end_index);
        [corr_output, lags] = xcorr(one_pulse,ref_pulse);
        corr_output = abs(corr_output)';

        valid_lag_idx = lags >= 0;
        lags = lags(valid_lag_idx);
        corr_output = corr_output(valid_lag_idx);

        echoProfiles = [echoProfiles; corr_output];
    end
    
    distance = lags * sec_to_sample * c / 2 * 100;

    %% STEP TWO - Identify moving object ----

    % computing threshold for motion change
    avg_corr = mean(echoProfiles(1:3,:), 1);
    [direct_amp, static_loc] = max(avg_corr);
    motion_threshold = motion_threshold_ratio * direct_amp;
        
    [numEchoProfiles, lagsPerEcho] = size(echoProfiles);
    
    mov_locs_1 = zeros(1,numEchoProfiles); 
    mov_locs_2 = zeros(1, numEchoProfiles);
    min_separation_samples = 28;
    
    for k = 1:numEchoProfiles-1
        diff = echoProfiles(k+1,:) - echoProfiles(k,:);
        found_first = false;
        for j = 1:lagsPerEcho
            if abs(diff(j)) > motion_threshold
                if ~found_first
                    mov_locs_1(k) = j;
                    found_first = true;
                else
                    if abs(j-mov_locs_1(k)) > min_separation_samples
                        mov_locs_2(k) = j;
                        break;
                    end
                end
            end
        end 
    end

    static_distance = distance(static_loc);
	motion_1 = NaN(1, numEchoProfiles);
	motion_2 = NaN(1, numEchoProfiles);
	for k = 1:numEchoProfiles
		if mov_locs_1(k) > 0
			motion_1(k) = distance(mov_locs_1(k));
		end
		if mov_locs_2(k) > 0
			motion_2(k) = distance(mov_locs_2(k));
		end
    end

    if isnan(motion_1(1))
        motion_1(1) = static_distance;
    end
    if isnan(motion_2(1))
        motion_2(1) = static_distance;
    end
		
	motion_1_filled = fillmissing(motion_1, 'previous');
	motion_1_filtered = medfilt1(motion_1_filled, 11); 
    motion_2_filled = fillmissing(motion_2, 'previous');
	motion_2_filtered = medfilt1(motion_2_filled, 7); 
		
    % PLOT
    figure(1);
    nexttile;
    hold on
    plot(1:numEchoProfiles, motion_1_filtered, 'o-');
    if multihand
        plot(1:numEchoProfiles, motion_2_filtered, 'o-');
    end
    xlabel('Pulse Index');
    ylabel('Distance (cm)');
    title(sprintf('Motion vs Time - Channel %d', i));
    grid on;
    hold off

    %% STEP THREE --- use OFDM properties for PRECISE trajectory ----- 
    ref_symbol = ref_pulse(preamble_len+cyclic_len+1:preamble_len+cyclic_len+symbol_len);    
    ref_symbol = hilbert(ref_symbol,symbol_len);
    ref_spectrum = fft(ref_symbol,symbol_len)'; 
    %ref_spectrum(symbol_len/2+2:end) = 0;

    precise_time_1 = NaN(1, numEchoProfiles);
    precise_time_2 = NaN(1, numEchoProfiles);
    f_bins = (active_bins - 1) * Fs / symbol_len;

    for k = 1:numEchoProfiles 

        % comms_start = locs(k) + preamble_len + cyclic_len + 1;
        % comms_end = comms_start + symbol_len - 1;
        % if comms_end > length(rx_one)
        %     continue;
        % end
        % 
        % comms_spectrum = fft(rx_one(comms_start:comms_end),symbol_len);
        % rx_bits = real(comms_spectrum(active_bins)) > 0;
        % disp('Received bits:');
        % disp(rx_bits);

        if mov_locs_1(k) > 0
            echo_start = locs(k) + mov_locs_1(k) - 1; 
            echo_end = echo_start + symbol_len - 1;
            if echo_end <= length(rx_one)
                echo_segment = rx_one(echo_start:echo_end);
                echo_segment = hilbert(echo_segment,symbol_len);
                echo_spectrum = fft(echo_segment,symbol_len)'; 
                phase_diff = angle(echo_spectrum) - angle(ref_spectrum);
                phase_diff = unwrap(phase_diff);
                p = polyfit(f_bins, phase_diff(active_bins), 1); 
                precise_time_1(k) = abs(p(1) / (2*pi));
            end
        end

        if mov_locs_2(k) > 0
            echo_start = locs(k) + mov_locs_2(k) - 1; 
            echo_end = echo_start + symbol_len - 1;
            if echo_end <= length(rx_one)
                echo_segment = rx_one(echo_start:echo_end);
                echo_segment = hilbert(echo_segment,symbol_len);
                echo_spectrum = fft(echo_segment,symbol_len)'; 
                phase_diff = angle(echo_spectrum) - angle(ref_spectrum);
                phase_diff = unwrap(phase_diff);
                p = polyfit(f_bins, phase_diff(active_bins), 1); 
                precise_time_2(k) = -p(1) / (2*pi);
            end
        end  
        
    end
		
	motion_precise_1 = NaN(1, numEchoProfiles);
	motion_precise_2 = NaN(1, numEchoProfiles);
    static_precise_distance = lags(static_loc) * sec_to_sample * c * 100 / 2;
	for k = 1:numEchoProfiles
		if ~isnan(precise_time_1(k))
			motion_precise_1(k) = precise_time_1(k) * c * 100 / 2;
		end
		if ~isnan(precise_time_2(k))
			motion_precise_2(k) = precise_time_2(k) * c * 100 / 2;
		end
	end
    
    if isnan(motion_precise_1(1))
        motion_precise_1(1) = static_precise_distance;
    end
    if isnan(motion_precise_2(1))
        motion_precise_2(1) = static_precise_distance;
    end
    motion_precise_1_filled = fillmissing(motion_precise_1, 'previous');
    motion_precise_1_filtered = medfilt1(motion_precise_1_filled, 15);
    motion_precise_2_filled = fillmissing(motion_precise_2, 'previous');
    motion_precise_2_filtered = medfilt1(motion_precise_2_filled, 7);
    
    figure(2);
    nexttile;
    plot(1:numEchoProfiles, motion_precise_1_filtered, 'o-');
    hold on;
    if multihand
        plot(1:numEchoProfiles,motion_precise_2_filtered,'o-');
        legend("closer", "farther");
    end
    xlabel('Pulse Index');
    ylabel('Precise Distance (cm)');
    title(sprintf('Precision Tracking - Channel %d',i));
    grid on;
    hold off
end