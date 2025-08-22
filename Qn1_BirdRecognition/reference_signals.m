% Plotting spectrum of given audio files

% clc, clearvars, close all;

% List of reference bird audio files
audio_files = {'./GivenSignals/Project_BirdRecognition/Reference/bird1.wav',...
    './GivenSignals/Project_BirdRecognition/Reference/bird2.wav',... 
    './GivenSignals/Project_BirdRecognition/Reference/bird3.wav'};
% audio_files = {'./Signals/Project_BirdRecognition/Reference/bird1.wav'};

% Loop through each file
for i = 1:length(audio_files)
    [audio_signal, Fs] = audioread(audio_files{i});
    
    % Convert to mono if stereo
    if size(audio_signal, 2) > 1
        audio_signal = mean(audio_signal, 2); % Average the two channels
    end
    
    % Compute FFT
    N = length(audio_signal);
    fft_signal = fft(audio_signal, N);
    f = (-N/2:N/2-1) * (Fs / N); % Frequency range: [-Fs/2, Fs/2)
    fft_shifted = fftshift(fft_signal);
    
    magnitude = abs(fft_shifted);
    
    Y = magnitude/max(magnitude);
    [peaks, locs] = findpeaks(Y, 'MinPeakHeight', 0.3, 'MinPeakDistance', 750);
    peakFreqs = f(locs);
   
    mainPeakFreq = peakFreqs(peaks == max(peaks));
    numSignificantPeaks = sum(peaks > 0.5);
    freqSpread = std(peakFreqs);
    
    % Plot magnitude spectrum
    figure;
    subplot(1, 1, 1);
    plot(f, magnitude);
    title(['Magnitude Spectrum of ', 'Bird - ' , num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
end