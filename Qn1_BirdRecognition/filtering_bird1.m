function filteredAudio = filtering_bird1(location)
    BandstopFilter = Bird1_Bandstop();
    
    [audioIn, fs] = audioread(location); 
    audioIn = audioIn(:, 1); 
    
    filteredAudio = filter(BandstopFilter, audioIn); 
    
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
    % title('Magnitude Spectrum of the Filtered Signal of Bird-1');
    % grid on;
end



