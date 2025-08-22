% clc, clearvars, close all;
% clc, clearvars;

function species = identifyBird(y, fs, test)    
    N = length(y);
    Y = fft(y, N);
    Y_shifted = fftshift(Y);
    Y = abs(Y_shifted);
    f = (-N/2:N/2-1) * (fs / N);
    
    % Normalize the spectrum
    Y = Y / max(Y);

    % figure;
    % plot(f, Y);
    % title(['Magnitude Spectrum after applying Filter of ', 'Bird - ' , num2str(test)]);
    % xlabel('Frequency (Hz)');
    % ylabel('Magnitude');

    % Find peaks in the spectrum
    [peaks, locs] = findpeaks(Y, 'MinPeakHeight', 0.3, 'MinPeakDistance', 750);
    peakFreqs = f(locs);
   
    mainPeakFreq = peakFreqs(peaks == max(peaks));
    numSignificantPeaks = sum(peaks > 0.5);
    freqSpread = std(peakFreqs);
    
    
    % Decision logic based on observed characteristics
    species = 0; % Bird Not familiar

    switch test
    case 1
        if any(abs(peakFreqs - 3068) < 100) && any(abs(peakFreqs - 6008) < 100)
            species = 1;  % Bird 1
        end
    case 2
        if freqSpread < 4000 && all(peakFreqs < 4500)
            species = 2;  % Bird 2
        end
    case 3
        if any(abs(peakFreqs - 7000) < 100) || any(abs(peakFreqs - 4128) < 100)
            species = 3;  % Bird 3
        end
    end
end

function analyzeReference()
    % Analyze reference files
    for idx = 1:3
        filename = sprintf('./GivenSignals/Project_BirdRecognition/Reference/bird%d.wav', idx);
        [~, fs] = audioread(filename);

        switch idx
        case 1
            filteredAudio = filtering_bird1(filename);
        case 2
            filteredAudio = filtering_bird2(filename);
        case 3
            filteredAudio = filtering_bird3(filename);
        end

        species = identifyBird(filteredAudio, fs, idx);
        fprintf('Reference Bird %d identified as Species %d\n', idx, species);
    end
end

function analyzeTask()
    files = dir('./GivenSignals/Project_BirdRecognition/Task/*.wav');
    
    % Analyze each task file
    for idx = 1:length(files)
        filename = fullfile('./GivenSignals/Project_BirdRecognition/Task', files(idx).name);
        [~, fs] = audioread(filename);
        ans = 0;

        for test = 1:3
            switch test
            case 1
                filteredAudio = filtering_bird1(filename);
            case 2
                filteredAudio = filtering_bird2(filename);
            case 3
                filteredAudio = filtering_bird3(filename);
            end

            ans = identifyBird(filteredAudio, fs, test);
            disp(['ans = ', num2str(ans)])
            if ans ~= 0
                break;
            end
        end
       
        fprintf('Task file %s identified as Bird- %d\n', files(idx).name, ans);
    end
end

analyzeTask();
% analyzeReference();