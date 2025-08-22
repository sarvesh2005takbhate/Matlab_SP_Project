clearvars, clc, close all;

% List of audio and text file names
audioFiles = {'./audios/1.wav', './audios/2.wav', './audios/3.wav', './audios/4.wav',  './audios/5.wav',  './audios/6.wav',  './audios/7.mp3',  './audios/8.mp3', './audios/9.mp3'};
textFiles = {'./text/1.txt', './text/2.txt', './text/3.txt',  './text/4.txt',  './text/5.txt',  './text/6.txt',  './text/7.txt',  './text/8.txt', './text/9.txt'};

for i = 1:length(audioFiles)

    [audioData, fs] = audioread(audioFiles{i});
    N = length(audioData);
    mag_spectra = fftshift(abs(fft(audioData)));

    % Convert to mono if stereo
    if size(audioData, 2) > 1
        audioData = mean(audioData, 2);
    end

    fileID = fopen(textFiles{i}, 'r');
    textData = textscan(fileID, '%s %f %f');
    fclose(fileID);
    
    words = textData{1}(1:2:length(textData{1}));
    startTimes = textData{2}(1:2:length(textData{2}));
    endTimes = textData{3}(1:2:length(textData{3}));
    
    allIntervals = [];
    for j = 1:length(endTimes)
        % Extract segment based on start and end times
        startSample = max(1, floor(startTimes(j) * fs));
        endSample = min(length(audioData), ceil(endTimes(j) * fs));
        segment = audioData(startSample:endSample);

        % Calculate loudness
            % Frequency Weighting
                aWeighting = weightingFilter('A-weighting', fs);
            % Filter the signal
                weighted_signal = aWeighting(segment);
                loudness = rms(weighted_signal);
        
        % Store results as [start time, end time, loudness]
        allIntervals = [allIntervals; [startTimes(j), endTimes(j), loudness]];
    end

    % Sort intervals by loudness in descending order
    sortedIntervals = sortrows(allIntervals, -3);
    
    % Calculate the average loudness
    averageLoudness = mean(sortedIntervals(:, 3));
    stdLoudness = std(sortedIntervals(:, 3));

    cutoff = averageLoudness + 0.5544*stdLoudness;

    % Add a new column: 1 if loudness > average, else 0
    aboveAverage = sortedIntervals(:, 3) > cutoff;
    sortedIntervals = [sortedIntervals, aboveAverage];

    % Display or save results
    disp('Start Time | End Time | Loudness | Is it Loud');
    disp(sortedIntervals);

    disp("dsod");

    % Write loudness to a file
        outputFile = ['Audiofile_' num2str(i) '.txt'];
        fileID = fopen(outputFile, 'w');
        fprintf(fileID, 'Start Time  |    End Time   |    Loudness |  Is it Loud\n');
        fprintf(fileID, '_______________________________________________________\n');
        fprintf(fileID, '%f \t %f \t %f \t %d\n', sortedIntervals');
        fclose(fileID);
end

% disp(['Mean = ', num2str(mean(sortedIntervals(:, 3)))]);
% disp(['Standard Deviation = ', num2str(std(sortedIntervals(:, 3)))]);