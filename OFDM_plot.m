clear
close all

[x ~] = audioread('ofdm_pulse_548.wav');
len = 548;

Fs = 48000;

% x = x(:,6);

t = (0:len-1)/Fs*1000;  % in milliseconds
figure(1);
plot(t, x(1:len));
xlabel('Time (ms)');
ylabel('Amplitude');
title('One OFDM Pulse (with cyclic suffix)');
grid on;

figure(2);
pspectrum(x(1:10*len),Fs,'spectrogram', ...
    'FrequencyLimits',[0 24000]);
title('Generated Chirp');
set(gca,'linewidth',2,'fontsize',26,'fontname','Arial');