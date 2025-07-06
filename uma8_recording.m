clear
close all

% Setup parameters
fs = 48000;         % Sampling rate (UMA-8 default)
frameLength = 2048; % Block size for reading
duration = 4;      % Recording duration in seconds
numChannels = 8;    % UMA-8 has 7 microphones + 1 ref channel (if present)

% Create device reader
deviceReader = audioDeviceReader('SampleRate', fs, ...
                                 'NumChannels', 8, ...
                                 'Device', 'micArray RAW SPK', ...
                                 'SamplesPerFrame', frameLength);

% Preallocate buffer
numFrames = floor((fs * duration) / frameLength);
audioData = zeros(numFrames * frameLength, numChannels);

% Record loop
disp('Recording...');
for i = 1:numFrames
    frame = deviceReader();  % Read samples from all channels
    audioData((i-1)*frameLength+1:i*frameLength, :) = frame;
end
disp('Recording finished.');

% %Compute time vector and average signal
% t = (0:size(audioData,1)-1) / fs;
% x_combined = mean(audioData, 2);
% 
% % Side-by-side layout: left = waveforms, right = spectrogram
% figure;
% set(gcf, 'Position', [200, 100, 1400, 700]);  % Wider figure for side-by-side layout
% 
% % Left column: channel waveforms
% for ch = 1:numChannels
%     subplot(numChannels, 2, (ch-1)*2 + 1);  % Left side of each row
%     plot(t, audioData(:, ch));
%     xlabel('Time (s)');
%     ylabel(sprintf('Ch %d', ch));
%     title(sprintf('Waveform - Channel %d', ch));
% end
% % 
% % Right column: large spectrogram (spans rows)
% subplot(1, 2, 2);  % Right half of figure
% spectrogram(x_combined, hamming(1024), 768, 1024, fs, 'yaxis');
% title('Spectrogram - Averaged Across 8 Channels');
% xlabel('Time (s)');
% ylabel('Frequency (kHz)');
% 
% sgtitle('Waveforms and Spectrogram Side-by-Side');
% 
% figure;
% spectrogram(x_combined, hamming(2048), 1536, 2048, fs, 'yaxis');
% title('Spectrogram - Averaged Across 8 Channels');
% xlabel('Time (s)');
% ylabel('Frequency (kHz)');
% ylim([15 22])

% [ref, fs] = audioread('hand gestures/18khz-today-uma8.wav');
% n = length(ref);
% f = (0:floor(n/2)-1) * fs / n;     % Frequency vector in Hz

% figure(1);
% set(gcf, 'Position', [100, 100, 1200, 600])f;
% tiledlayout(2,4);
% for i = 1:8
%     nexttile;
%     one_channel = audioData(:,i);
%     spectrogram(one_channel, hamming(2048), 1536, 2048, fs, 'yaxis');
%     title(sprintf('Spectrogram - Channel %d', i));
%     xlabel('Time (s)');
%     ylabel('Frequency (kHz)');
%     % xlim([.05 0.35]);
%     % ylim([15 24]);
% end

% figure(2);
% set(gcf, 'Position', [200, 100, 1200, 600]);
% tiledlayout(2,4);
% for i = 1:8
%     nexttile;
%     one_channel = audioData(:,i);
%     % ref_one_channel = ref(:,i);
% 
%     Y = mag2db(abs(fft(one_channel)));
%     % Y_ref = mag2db(abs(fft(ref_one_channel)));
%     Y_half = Y(1:floor(n/2));
%     % ref_half = Y_ref(1:floor(n/2));
% 
%     hold on
%     % plot(f, ref_half);
%     f = (0:floor(length(one_channel)/2)-1) * fs / n;     % Frequency vector in Hz
%     plot(f, Y_half);
%     % legend('18khz pure', 'moving');
%     grid on
%     xlabel("Frequency (Hz)");
%     xlim([15000 24000]);
%     ylabel("Decibel");
%     title(sprintf('Fourier Transform - Channel %d', i));
% end

%audiowrite('hand gestures/chirp21-stationaryHighBySpeaker-uma8.wav', audioData, fs);

FMCWProcessing(audioData);