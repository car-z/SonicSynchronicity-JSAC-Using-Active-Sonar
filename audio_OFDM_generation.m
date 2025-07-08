clc;
clear;
close all;   

%% OFDM Paramters
Fs = 48000; % sampling rate (Hz)
total_time = 120; % duration of signal (seconds)
symbol_len = 128; % OFDM symbol length (samples)
cyclic_suffix_len = 20; % (samples)
separation_len = 400; % each pulse separated by 400 samples
len = symbol_len + cyclic_suffix_len + separation_len; % total length of pulse (samples)
T = len / Fs; % duration of one pulse (seconds)
num_repititions = floor(total_time/T); % number of OFDM pulses per signal

% ------ Fading Configuration ------
% fading_ratio = 0.2;

%% Frequency Parameters
subcarrier_width = Fs / (symbol_len); % = 375 Hz
f = (0:symbol_len-1) * subcarrier_width; % vector spanning 0-48 kHz in 128 subcarriers

% identify indices in f corresponding to 18-20 kHz
lower_k = ceil(18000 / subcarrier_width);
upper_k = floor(20000 / subcarrier_width);
active_bins = lower_k+1 : upper_k + 1; % MATLAB uses 1 based indexing

%% IFFT to get time-domain OFDM symbol
X = zeros(1, symbol_len);
win = hann(length(active_bins))';
X(active_bins) = win .* (-1).^(0:length(active_bins)-1);

% X(active_bins) = (-1).^(0:length(active_bins)-1);
x_complex = ifft(X);
x_real = real(x_complex);

%% append cyclic suffix and separation
cyclic_suffix = x_real(1:cyclic_suffix_len);
ofdm_pulse = [x_real, cyclic_suffix]';  % 168 samples total

% % ofdm_smooth = ApplyFadingInStartAndEndOfSignal(ofdm_pulse, fading_ratio);

separation = zeros(separation_len,1);
full_pulse = [ofdm_pulse; separation];

%% generate repeated OFDM pulses
signal = repmat(full_pulse, 1, num_repititions);
signal = reshape(signal, [], 1);

%% audiowrite
signal = signal / max(abs(signal));  % Normalize to full scale
% audiowrite('ofdm_pulse_548.wav', signal, Fs);

%% plot
t = (0:len-1)/Fs*1000;  % in milliseconds
figure(1);
plot(t, full_pulse);
xlabel('Time (ms)');
ylabel('Amplitude');
title('One OFDM Pulse (with cyclic suffix)');
grid on;

figure(2);
pspectrum(signal(1:10*len),Fs,'spectrogram', ...
    'FrequencyLimits',[0 24000]);
title('Generated Chirp');
set(gca,'linewidth',2,'fontsize',26,'fontname','Arial');