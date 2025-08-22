function filteredAudio = filtering_bird3(location)
    BandpassFilter1 = Bird3_Bandpass1;  
    BandpassFilter2 = Bird3_Bandpass2;
    BandpassFilter3 = Bird3_Bandpass3;
    
    [audioIn, fs] = audioread(location);    
    audioIn = audioIn(:, 1);    
    
    filteredAudio1 = filter(BandpassFilter1, audioIn);
    filteredAudio2 = 2*filter(BandpassFilter2, audioIn);
    filteredAudio3 = filter(BandpassFilter3, audioIn);
    filteredAudio = filteredAudio1 + filteredAudio2 + filteredAudio3;
    
    N = length(filteredAudio);        
    Y = fft(filteredAudio, N);        
    Y = fftshift(Y);                  
    magnitudeSpectrum = abs(Y);       
    freq = fs * (-N/2:N/2-1) / N;     
    
    % figure;
    % subplot(2,1,1);
    % plot(freq, magnitudeSpectrum);
    % xlabel('Frequency (Hz)');
    % ylabel('Magnitude');
    % title('Magnitude Spectrum of the Filtered Signal of Bird-3');
    % grid on;
end

