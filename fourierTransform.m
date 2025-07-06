clear
close all

[ref, fs] = audioread('hand gestures/18khz-uma8.wav');
[x, fs] = audioread('hand gestures/fingerSingleCircle-uma8.wav');

n = length(ref);
f = (0:floor(n/2)-1) * fs / n;     % Frequency vector in Hz

% 
% figure(1);
% set(gcf, 'Position', [100, 100, 1200, 600]);
% tiledlayout(2,4);
% for i = 1:8
%     nexttile;
%     one_channel = x(:,i);
%     spectrogram(one_channel, hamming(2048), 1536, 2048, fs, 'yaxis');
%     title(sprintf('Spectrogram - Channel %d', i));
%     xlabel('Time (s)');
%     ylabel('Frequency (kHz)');
%     ylim([15 22]);
% end
% 
% figure(2);
% set(gcf, 'Position', [200, 100, 1200, 600]);
% tiledlayout(2,4);
% for i = 1:8
%     nexttile;
%     one_channel = x(:,i);
%     ref_one_channel = ref(:,i);
% 
%     Y = mag2db(abs(fft(one_channel)));
%     Y_ref = mag2db(abs(fft(ref_one_channel)));
%     Y_half = Y(1:floor(n/2));
%     ref_half = Y_ref(1:floor(n/2));
% 
%     hold on
%     plot(f, ref_half);
%     plot(f, Y_half);
%     legend('18khz pure', 'moving');
%     grid on
%     xlabel("Frequency (Hz)");
%     xlim([15000 22000]);
%     ylabel("Decibel");
%     title(sprintf('Fourier Transform - Channel %d', i));
% end