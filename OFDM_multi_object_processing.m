clear
close all

% .wav file of received signal
filename = 'test.wav';

% ratio to use as threshold for motion
motion_threshold_ratio = .01; % NEED TO ADJUST

multihand = false; % TRUE if trying to process multiple moving objects
% processing multiple moving objects doesn't quite work yet

%% Constants -----------------
Fc = 18000;         % signal starting frequency: 18Hz
B = 2000;           % signal bandwidth (2000 Hz)
Fs = 48000;         % sampling rate (48 kHz)
preamble_len = 64;  % # samples in preamble
cyclic_len = 20;    % # samples in cyclic prefix
symbol_len = 128;   % # samples in actual OFDM symbol
separation_len = 200; % # samples of separation between each OFDM pulse

pulse_len = preamble_len + cyclic_len + symbol_len + separation_len;

c = 343; % speed of sound in m/s
sec_to_sample = 1 / Fs;

subcarrier_width = Fs / (symbol_len); % width of each OFDM subcarrier in Hz
% calculating which subcarriers are of interest (desired frequency range)
lower_k = ceil(Fc / subcarrier_width);
upper_k = floor((Fc+B) / subcarrier_width);
active_bins = lower_k+1 : upper_k + 1;

%% Transmitted Signal ---------------
% read in one pulse of the transmitted signal
[tx, ~] = audioread('OFDM signals/ofdm_pulse_348_prefix_preamble.wav');
ref_pulse = tx(1:pulse_len);

%% Received Signal --------------
[rx, ~] = audioread(filename);

rx = rx(:,1:7); % removing empty channel on uma-8 mic

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

% for every channel on the UMA-8 MIC
for i = 1:7
    rx_one = rx_filtered(:,i);

    %% Detect Indices of Pulse Starts
    % use only the preamble to find pulse starts
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

	% ALTERNATE METHOD - FOUND THROUGH TESTING TO BE MORE ACCURATE
    corr = conv(rx_one, flipud(ref_pulse));
    [~, firstPeak] = findpeaks(corr, 'NPeaks', 1);
    num_pulses = floor(length(rx_one) / pulse_len);
    % calculate index positions of the start of all pulses
    locs = firstPeak + (0:num_pulses - 1) * pulse_len;

    %% ----- Generate ECHO PROFILES for all reflections at mic ----------------
    echoProfiles = []; % each row is an echo profile for a pulse in the received signal
    for k = 1:num_pulses
        start_index = locs(k);
        end_index = start_index + pulse_len - 1;
        if end_index > length(rx_one)
            continue
        end
        one_pulse = rx_one(start_index:end_index); % extract singular pulse
        % calculate cross-correlation between transmitted and received
        % pulse
        [corr_output, lags] = xcorr(one_pulse,ref_pulse);
        corr_output = abs(corr_output)';

        % only care about DELAYS in time (lags >= 0)
        valid_lag_idx = lags >= 0;
        lags = lags(valid_lag_idx);
        corr_output = corr_output(valid_lag_idx);

        echoProfiles = [echoProfiles; corr_output];
    end
    
    % map delays to distance
    distance = lags * sec_to_sample * c / 2 * 100;

    %% --- Identify moving object(s) based on echo profile ----

    % computing threshold for motion change
    avg_corr = mean(echoProfiles(1:3,:), 1);
    [direct_amp, static_loc] = max(avg_corr);
    motion_threshold = motion_threshold_ratio * direct_amp;
        
    [numEchoProfiles, lagsPerEcho] = size(echoProfiles);
    
    % arrays to hold the index within each pulse where the moving object is
    mov_locs_1 = zeros(1,numEchoProfiles); 
    mov_locs_2 = zeros(1, numEchoProfiles);
    min_separation_samples = 28; % FINE TUNE THIS BASED ON HOW FAR APART THE TWO MOVING OBJECTS ARE
    
    for k = 1:numEchoProfiles-1
        % calculate difference between consecutive echo profiles
        diff = echoProfiles(k+1,:) - echoProfiles(k,:);
        found_first = false;
        for j = 1:lagsPerEcho
            % if difference is greater than the motion threshold -> moving
            % object is found
            if abs(diff(j)) > motion_threshold
                if ~found_first
                    mov_locs_1(k) = j;
                    found_first = true;
                else
                    % if the distance between this location and the
                    % first moving object is > min_separation_samples than
                    % another distinct moving object has been found
                    if abs(j-mov_locs_1(k)) > min_separation_samples
                        mov_locs_2(k) = j;
                        break;
                    end
                end
            end
        end 
    end

    % finding distance of largest reflector (baseline distance)
    static_distance = distance(static_loc);
    % arrays to hold actual distances of moving objects
	motion_1 = NaN(1, numEchoProfiles);
	motion_2 = NaN(1, numEchoProfiles);
	for k = 1:numEchoProfiles
		if mov_locs_1(k) > 0
            % map lag index within each pulse to distance value
			motion_1(k) = distance(mov_locs_1(k));
		end
		if mov_locs_2(k) > 0
			motion_2(k) = distance(mov_locs_2(k));
		end
    end

    % set initial position, if not already
    if isnan(motion_1(1))
        motion_1(1) = static_distance;
    end
    if isnan(motion_2(1))
        motion_2(1) = static_distance;
    end
	
    % smooth
	motion_1_filled = fillmissing(motion_1, 'previous');
	motion_1_filtered = medfilt1(motion_1_filled, 11); 
    motion_2_filled = fillmissing(motion_2, 'previous');
	motion_2_filtered = medfilt1(motion_2_filled, 7); 

    t_motion = locs(1:numEchoProfiles) / Fs; % time array
    rel_motion_1 = motion_1_filtered - motion_1_filtered(1); % adjust by initial position
		
    % PLOT
    figure(1);
    nexttile;
    hold on
    plot(t_motion, rel_motion_1, 'o-');
    if multihand
        plot(t_motion, motion_2_filtered, 'o-');
    end
    xlabel('Time (s)','FontSize',18);
    ylabel('Motion (cm)','FontSize',18);
    title(sprintf('Motion vs Time - Channel %d', i));
    grid on;
    hold off

    %% --- ADDITIONAL: use OFDM properties for PRECISE trajectory ----- 
    % extract the SYMBOL only (no preamble, no cyclic prefix)
    ref_symbol = ref_pulse(preamble_len+cyclic_len+1:preamble_len+cyclic_len+symbol_len);    
    ref_symbol = hilbert(ref_symbol,symbol_len); % turns real signal into complex analytic
    ref_spectrum = fft(ref_symbol,symbol_len)'; % perform fourier transform to get frequencies

    % arrays to hold the precise timings
    precise_time_1 = NaN(1, numEchoProfiles);
    precise_time_2 = NaN(1, numEchoProfiles);
    f_bins = (active_bins - 1) * Fs / symbol_len; % actual used frequencies

    for k = 1:numEchoProfiles 

        % using received OFDM symbol to derive the transmitted message
        comms_start = locs(k) + preamble_len + cyclic_len + 1;
        comms_end = comms_start + symbol_len - 1;
        if comms_end > length(rx_one)
            continue;
        end
        comms_spectrum = fft(rx_one(comms_start:comms_end),symbol_len);
        rx_bits = real(comms_spectrum(active_bins)) > 0;
        disp('Received bits:');
        disp(rx_bits);

        % if movement occurred during this echo profile
        if mov_locs_1(k) > 0
            % starting index of echo of the movement
            echo_start = locs(k) + mov_locs_1(k) - 1; 
            echo_end = echo_start + symbol_len - 1;
            if echo_end <= length(rx_one)
                % extract movement echo 
                echo_segment = rx_one(echo_start:echo_end);
                % convert the real signal into a complex anlytic signal
                echo_segment = hilbert(echo_segment,symbol_len);
                % take FFT of the signal
                echo_spectrum = fft(echo_segment,symbol_len)'; 
                % calculate the phase difference 
                % between the received and transmitted signals
                phase_diff = angle(echo_spectrum) - angle(ref_spectrum);
                phase_diff = unwrap(phase_diff);

                % phase shift as a function of frequency = -2pi*f*delay in time
                 % therefore: delay in time 
                % = phase shift as function of frequency/(frequency*-2pi)
                % delay in time = (slope of phase shift vs f)/ -2pi
                p = polyfit(f_bins, phase_diff(active_bins), 1); % finding phase shift vs frequency function
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

    % calculating the *precise* distance of the closest static reflector
    static_precise_distance = lags(static_loc) * sec_to_sample * c * 100 / 2;
    % mapping precise time to precise distance
	for k = 1:numEchoProfiles
		if ~isnan(precise_time_1(k))
			motion_precise_1(k) = precise_time_1(k) * c * 100 / 2;
		end
		if ~isnan(precise_time_2(k))
			motion_precise_2(k) = precise_time_2(k) * c * 100 / 2;
		end
	end
    
    % setting initial position
    if isnan(motion_precise_1(1))
        motion_precise_1(1) = static_precise_distance;
    end
    if isnan(motion_precise_2(1))
        motion_precise_2(1) = static_precise_distance;
    end

    % smooth
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
    title(sprintf('Precise Motion vs Time - Channel %d',i));
    grid on;
    hold off
end