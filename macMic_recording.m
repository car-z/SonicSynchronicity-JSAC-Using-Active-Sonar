clear
close all

% Setup parameters
fs = 48000;         % Sampling rate (UMA-8 default)
frameLength = 2048; % Block size for reading
duration = 2;      % Recording duration in seconds
numChannels = 1;    % UMA-8 has 7 microphones + 1 ref channel (if present)

% Create device reader
deviceReader = audioDeviceReader('SampleRate', fs, ...
                                 'NumChannels', 1, ...
                                 'Device', 'Macbook Air Microphone', ...
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

figure(1);
spectrogram(audioData(:,1), hamming(2048), 1536, 2048, fs, 'yaxis');
title(sprintf('Spectrogram - Channel %d', i));
xlabel('Time (s)');
ylabel('Frequency (kHz)');
% ylim([15 22]);

audiowrite('hand gestures/chirp-pushIn-macMic.wav', audioData, fs);