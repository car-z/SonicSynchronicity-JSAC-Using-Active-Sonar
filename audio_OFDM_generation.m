clc;
clear;
close all;   

%% OFDM Paramters
Fs = 48000; % sampling rate (Hz)
total_time = 120; % duration of signal (seconds)
symbol_len = 128; % OFDM symbol length (samples)
cyclic_len = 20; % (samples)
separation_len = 200; % each pulse separated by 400 samples
preamble_len = 64; % samples, adjust as needed
len = symbol_len + cyclic_len + separation_len; % total length of pulse (samples)
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
% X(active_bins) = win .* (-1).^(0:length(active_bins)-1);

X(active_bins) = (-1).^(0:length(active_bins)-1);
x_complex = ifft(X);
x_real = real(x_complex);

%% append cyclic prefix and separation
% cyclic_suffix = x_real(1:cyclic_suffix_len);
% ofdm_pulse = [x_real, cyclic_suffix]';  % 168 samples total

cyclic_prefix = x_real(end-cyclic_len+1:end);

% % Generate Zadoff-Chu sequence
% % zc_root = 25;
% % zc_len = length(active_bins);
% % zc_freq = zadoffChuSeq(zc_root, zc_len).'; % complex - zc seq in freq domain
% % ZC_vector = zeros(1,symbol_len);
% % ZC_vector(active_bins) = zc_freq;
% % zc_preamble = ifft(ZC_vector, symbol_len); % ZC in time domain
% 
% % Zadoff-Chu preamble
% zc_len = 63; % must be odd
% zc_root = 25;
% zc_base = zadoffChuSeq(zc_root, zc_len);
% 
% % Bandshift ZC to ~19 kHz
% fc = 19000; % center freq
% n = (0:zc_len-1)';
% zc_shifted = zc_base .* exp(1j * 2 * pi * fc * n / Fs);
% 
% % Take real part
% zc_real = real(zc_shifted)';

ofdm_pulse = [cyclic_prefix, x_real]';

% % ofdm_smooth = ApplyFadingInStartAndEndOfSignal(ofdm_pulse, fading_ratio);

separation = zeros(separation_len,1);
full_pulse = [ofdm_pulse; separation];

%% Generate strong known preamble: linear chirp 18–20 kHz
f_start = 18000; % Hz
f_end = 20000;   % Hz
t_preamble = (0:preamble_len-1) ./ Fs;
preamble = chirp(t_preamble, f_start, t_preamble(end), f_end, 'linear')';
preamble = ApplyFadingInStartAndEndOfSignal(preamble, 0.8);
    
pulse_and_preamble = [preamble; full_pulse];

%% generate repeated OFDM pulses
signal = repmat(pulse_and_preamble, 1, num_repititions);
signal = reshape(signal, [], 1);

%% audiowrite
signal = signal / max(abs(signal));  % Normalize to full scale
%audiowrite('OFDM signals/ofdm_pulse_348_prefix_preamble_no_fading.wav', signal, Fs);

%% plot
t = (0:len + preamble_len -1)/Fs*1000;  % in milliseconds
figure(1);
plot(t, pulse_and_preamble);
xlabel('Time (ms)');
ylabel('Amplitude');
title('One OFDM Pulse (with cyclic prefix)');
grid on;

figure(2);
pspectrum(signal(1:10*(preamble_len +len)),Fs,'spectrogram', ...
    'FrequencyLimits',[0 24000]);
title('Generated Chirp');
set(gca,'linewidth',2,'fontsize',26,'fontname','Arial');