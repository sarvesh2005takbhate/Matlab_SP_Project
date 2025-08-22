function filteredAudio = filtering_bird2(location)
    BandpassFilter1 = Bird2_Bandpass1;  
    BandpassFilter2 = Bird2_Bandpass2;

    [audioIn, fs] = audioread(location);    
    audioIn = audioIn(:, 1);    

    filteredAudio1 = filter(BandpassFilter1, audioIn);  
    filteredAudio2 = filter(BandpassFilter2, audioIn);
    filteredAudio = filteredAudio1 + filteredAudio2;
        
    N = length(filteredAudio);          
    Y = fft(filteredAudio, N);          
    Y = fftshift(Y);                    
    magnitudeSpectrum = abs(Y);         
    freq = fs * (-N/2:N/2-1) / N;    
    
    % figure;
    % subplot(1,1,1);
    % plot(freq, magnitudeSpectrum);
    % xlabel('Frequency (Hz)');
    % ylabel('Magnitude');
    % title('Magnitude Spectrum of the Filtered Signal of Bird-2');
    % grid on;
end

