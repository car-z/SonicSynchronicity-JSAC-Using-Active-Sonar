clear;
close all;

[x,fs] = audioread('hand gestures/chirp21-upDownBySpeaker-uma8.wav');

figure(1);
set(gcf, 'Position', [100, 100, 1200, 600]);
tiledlayout(2,4);
for i = 1:8
    nexttile;
    one_channel = x(:,i);
    spectrogram(one_channel, hamming(2048), 1536, 2048, fs, 'yaxis');
    title(sprintf('Spectrogram - Channel %d', i));
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');
    % xlim([.05 0.35]);
    % ylim([15 24]);
end