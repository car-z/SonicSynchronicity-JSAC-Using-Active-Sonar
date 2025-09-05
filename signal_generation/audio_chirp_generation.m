clc;
clear;
close all;   

%% Global Paramters
% ------ FMCW configuration ------
Fc = 18000; % start frequency of chirp (Hz)
B = 4000; % bandwidth of chirp 
Fs = 48000; % sampling rate (Hz)
T = 0.01067; % duration of one chirp (seconds)
total_time = 120; % duration of signal (seconds)
single_chirp_len = 512; % samples per chirp
num_of_repititions = floor(total_time/T); % number of chirps per signal

% ------ Fading Configuration ------
fading_ratio = 0.2;

%% Transmitter
time = (0:single_chirp_len-1)./Fs;
trans_sw_fmcw = chirp(time, Fc, time(end), Fc+B).';

trans_sw_fmcw = ApplyFadingInStartAndEndOfSignal(trans_sw_fmcw, fading_ratio);

trans_sw_fmcw = repmat(trans_sw_fmcw,1,num_of_repititions);
trans_sw = reshape(trans_sw_fmcw,[],1);


% audiowrite('chirp10.wav',trans_sw,Fs);

figure;
pspectrum(trans_sw(1:10*single_chirp_len),Fs,'spectrogram', ...
    'FrequencyLimits',[0 24000]);
title('Generated Chirp');
set(gca,'linewidth',2,'fontsize',26,'fontname','Arial');

player = audioplayer(trans_sw, Fs);
play(player);