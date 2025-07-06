clc;
clear;
close all;   

%% OFDM Paramters
Fs = 48000; % sampling rate (Hz)
T = 0.00592; % duration of one pulse (seconds)
total_time = 120; % duration of signal (seconds)
symbol_len = 128; % OFDM symbol length (samples)
cyclic_prefix_len = 40; % (samples)
separation_len = 400; % each pulse separated by 400 samples
len = symbol_len + cyclic_prefix_len + separation_len; % total length of pulse (samples)
num_repititions = floor(total_time/T); % number of OFDM pulses per signal

% ------ Fading Configuration ------
fading_ratio = 0.2;

%% Frequency Parameters
subcarrier_width = Fs / symbol_len; % = 375 Hz
f = (0:symbol_len-1) * subcarrier_width; % vector spanning 0-24 kHz in 64 subcarriers

% identify indices in f corresponding to 18-20 kHz
lower_k = ceil(18000 / subcarrier_width);
upper_k = floor(20000 / subcarrier_width);
active_bins = lower_k+1 : upper_k + 1; % MATLAB uses 1 based indexing

%% Hermitian-symmetric frequency-domain vector
X = zeros(1, symbol_len);
X(active_bins) = (-1).^(0:length(active_bins)-1);

% Now apply Hermitian symmetry
for k = active_bins
    mirror_idx = symbol_len - k + 2;
    X(mirror_idx) = conj(X(k));
end

% % Zero out other unwanted bins
% X(1:lower_k) = 0;            % Remove low-freq bins
% X(symbol_len/2 + 2:end) = 0; % Remove post-Nyquist bins
% X(1) = 0;                    % Remove DC

%% IFFT to get time-domain OFDM symbol
% win = hann(length(active_bins))';
% X(active_bins) = win .* (-1).^(0:length(active_bins)-1);
% X(1:lower_k) = 0;        % zero out low bins
% X(symbol_len/2 + 2:end) = 0; % (optional extra check)
% X(1) = 0;                % kill DC

x_complex = ifft(X, symbol_len);
x_real = real(x_complex);

%% append cyclic prefix and separation
cyclic_prefix = x_real(1:cyclic_prefix_len);
ofdm_pulse = [cyclic_prefix, x_real]';  % 84 samples total

ofdm_smooth = ApplyFadingInStartAndEndOfSignal(ofdm_pulse, fading_ratio);

separation = zeros(separation_len,1);
full_pulse = [ofdm_smooth; separation];

%% generate repeated OFDM pulses
signal = repmat(full_pulse, 1, num_repititions);
signal = reshape(signal, [], 1);

%% audiowrite
signal = signal / max(abs(signal));  % Normalize to full scale
audiowrite('ofdm_pulse_rand.wav', signal, Fs);

%% plot
t = (0:length(full_pulse)-1)/Fs*1000;  % in milliseconds
figure(1);
plot(t, full_pulse);
xlabel('Time (ms)');
ylabel('Amplitude');
title('One OFDM Pulse (with cyclic prefix)');
grid on;

figure(2);
pspectrum(signal(1:10*len),Fs,'spectrogram', ...
    'FrequencyLimits',[0 24000]);
title('Generated Chirp');
set(gca,'linewidth',2,'fontsize',26,'fontname','Arial');