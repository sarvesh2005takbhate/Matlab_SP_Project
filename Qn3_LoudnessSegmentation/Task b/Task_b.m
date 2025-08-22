clearvars, clc, close all;

% List of audio and text file names
audioFiles = {'./audios/1.wav', './audios/2.wav', './audios/3.wav', './audios/4.wav',  './audios/5.wav',  './audios/6.wav',  './audios/7.mp3',  './audios/8.mp3', './audios/9.mp3'};

% audioFiles = {'./audios/2.wav'};

for k = 1:length(audioFiles)

    [audioData, fs] = audioread(audioFiles{k});
    % If stereo, convert to mono
    if size(audioData, 2) > 1
        audioData = mean(audioData, 2);
    end

    audioData = audioData';

    N = length(audioData);
    mag_spectra = fftshift(abs(fft(audioData)));

    % Defining segment/window
    segment = ones(1, ceil(N/1e3));
    seglen = length(segment);

    seg_energy = zeros(1, N);

    % Append zeros of length seglen-1 to allow overflow of window (slide all the way to the end)
    sq_signal = [audioData.^2 zeros(1, seglen-1)];

    for ix = 1:N
        window = sq_signal(ix:ix+seglen-1);
        seg_energy(ix) = sum(window);
        avg_seg_energy(ix) = mean(window);
        seg_peak(ix) = max(window);
    end

    threshold = 0.05*mean(seg_energy);  % determine whether a word or not
    mingap = 0.2;   % minimum gap between two words. Equivalently, maximum gap between two syllables
    mindur = 0.1;   % minimum duration of the word, i.e, shortest word duration
    starttimes = [0];
    endtimes = [0];

    % Logic to segregate noise and speech (segmentation)
    for ix = 2:N
        if (seg_energy(ix) > threshold && seg_energy(ix-1) <= threshold && abs(ix/fs - starttimes(length(starttimes)))>mingap)
            if (length(starttimes)==length(endtimes))
                starttimes = [starttimes; ix/fs];
            end 
        elseif (seg_energy(ix) < threshold && seg_energy(ix-1) >= threshold && abs(ix/fs - endtimes(length(endtimes)))>mingap)
            if (ix/fs - starttimes(length(starttimes)) > mindur)
                if (length(endtimes)+1==length(starttimes))
                    endtimes = [endtimes; ix/fs];
                end
            end
        end
    end
    
    if (length(endtimes) < length(starttimes))
        endtimes = [endtimes; N/fs];
    end

    starttimes = starttimes(2:length(starttimes));
    endtimes = endtimes(2:length(endtimes));

    % Initialize loudness array
        loudness = zeros(length(starttimes), 1);

        for iy = 1:length(starttimes)
            % Calculate loudness values for each segment
                % Frequency Weighting
                    aWeighting = weightingFilter('A-weighting', fs);
                % Filter the signal
                weighted_signal = aWeighting((audioData(fs*starttimes(iy):fs*endtimes(iy))));
                loudness = rms(weighted_signal);
        end

        for i = 1:length(starttimes)
            loudness(i) = rms(audioData(starttimes(i)*fs:endtimes(i)*fs));
        end

        loudness = loudness';
        matrix = [starttimes endtimes loudness];
        sortrows(matrix, 3, 'descend');

        cutoff = mean(loudness) + 1.2*std(loudness);
        isloud = zeros(length(loudness), 1);

        counter = 0;
        for iz = 1:length(loudness)
            if (loudness(iz) > cutoff)
                isloud(iz) = 1;
            else
                isloud(iz) = 0;
            end
        end

        disp("File "+ num2str(k));
        matrix = [matrix isloud];
        matrix = sortrows(matrix, 3, 'descend');
        matrix

    figure;
        subplot(2, 1, 1)
            plot(linspace(0, N/fs, N), audioData);
            title("Time domain Speech signal of file " + num2str(k));
        subplot(2, 1, 2)
            plot(linspace(0, N/fs, N), seg_energy);
            hold on;
            plot(linspace(0, N/fs, N), threshold*ones(1, N));
            hold off;
            title("Energy of segment of file " + num2str(k));
end