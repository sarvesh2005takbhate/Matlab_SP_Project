% Plotting spectrum of given audio files

% clc, clearvars, close all;

% List of task bird audio files
audio_files = {'./GivenSignals/Project_BirdRecognition/Task/F1.wav',...
    './GivenSignals/Project_BirdRecognition/Task/F2.wav',... 
    './GivenSignals/Project_BirdRecognition/Task/F3.wav',...
    './GivenSignals/Project_BirdRecognition/Task/F4.wav',...
    './GivenSignals/Project_BirdRecognition/Task/F5.wav',...
    './GivenSignals/Project_BirdRecognition/Task/F6.wav',...
    './GivenSignals/Project_BirdRecognition/Task/F7.wav',...
    './GivenSignals/Project_BirdRecognition/Task/F8.wav'};
Fs = 0;

function plotting_task(start_idx, end_idx, audio_files)
    for i_loop = start_idx:end_idx
    idx = mod(i_loop-1, 4) +1;
    % if idx > 4
    %     figure;
    %     flag = 0;
    %     idx = idx - 4;
    %     disp(idx);
    % end
    
    [audio_signal, Fs] = audioread(audio_files{i_loop});

    if size(audio_signal, 2) > 1
        audio_signal = mean(audio_signal, 2); 
    end
    
    % Compute FFT
    N = length(audio_signal);
    fft_signal = fft(audio_signal, N);
    f = (-N/2:N/2-1) * (Fs / N); % Frequency range: [-Fs/2, Fs/2)
    fft_shifted = fftshift(fft_signal);
    magnitude = abs(fft_shifted);

    % Plot magnitude spectrum
    subplot(2, 2, idx);
    plot(f, magnitude);
    title(['Magnitude Spectrum of ', 'F - ' , num2str(i_loop)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    
    % Frequency-Domain Features
    N = length(audio_signal);
    fft_signal = fft(audio_signal, N);
    fft_magnitude = abs(fft_signal(1:N/2)); 
    f = (0:N/2-1) * (Fs / N); 
    
    % Spectral centroid
    spectral_centroid = sum(f .* fft_magnitude') / sum(fft_magnitude); 
    
    % Loudness (RMS)
    rms_power = sqrt(mean(audio_signal.^2));
    
    % Display Results
    fprintf('File: %s\n', num2str(i_loop));
    fprintf('Spectral Centroid: %.2f Hz\n', spectral_centroid);
    fprintf('RMS Power: %.2f\n\n', rms_power);
    end
end

figure;
plotting_task(1, 4, audio_files);
figure;
plotting_task(5, 8, audio_files);
