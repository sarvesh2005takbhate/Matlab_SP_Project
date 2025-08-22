close all;
clearvars;

fs = 128;
% Sampling rate is given to be 128Hz

Ts = 1/fs;
ecg1 = load('E1.mat').E1;
ecg2 = load('E2.mat').E2;
ecg3 = load('E3.mat').E3;

% Vectors used for plotting the first 10 seconds of time domain signals
axis = 0:1279;
first_ten = axis+1;

% X axis for 100001 point FFT analysis
N = 100001;
w_axis = -pi:2*pi/N:pi-2*pi/N;
f_axis = w_axis*64/pi;

%Plotting the first ten seconds of the raw signals
figure;
subplot(3,1,1);
plot(axis/fs,ecg1(first_ten));
title('Clean E1 signal');
xlabel('Time(seconds)');
ylabel('Amplitude')
subplot(3,1,2);
plot(axis/fs,ecg2(first_ten));
title('Noisy E2 signal');
xlabel('Time(seconds)');
ylabel('Amplitude')
subplot(3,1,3);
plot(axis/fs,ecg3(first_ten));
title('Noisy E3 signal');
xlabel('Time(seconds)');
ylabel('Amplitude')

% Frequency spectra of raw 
figure;
subplot(3,1,1);
plot(f_axis,abs(fftshift(fft(ecg1,N))));
title('Unfiltered spectrum of E1 signal');
xlabel('Frequency(Hz)')
ylabel('Magnitude')
xlim([-64 64])
subplot(3,1,2);
plot(f_axis,abs(fftshift(fft(ecg2,N))));
title('Unfiltered spectrum of E2 signal');
xlabel('Frequency(Hz)')
ylabel('Magnitude')
xlim([-64 64])
subplot(3,1,3);
plot(f_axis,abs(fftshift(fft(ecg3,N))));
title('Unfiltered spectrum of E3 signal');
xlabel('Frequency(Hz)')
ylabel('Magnitude')
xlim([-64 64])

%High pass filtering to remove base band "wander" noise - cutoff frequency 2Hz
% Using a 100 coefficient FIR filter

f_cutoff = 2;
order = 100;
b = fir1(order, f_cutoff / (fs / 2), 'high'); 

ecg1 = filtfilt(b, 1, ecg1);
ecg2 = filtfilt(b, 1, ecg2);
ecg3 = filtfilt(b, 1, ecg3);

%Low pass filtering to remove high frequency noise
% Cutoff frequency 40Hz, implemented with Butterworth filter
[c, d] = butter(30, 40/(fs/2), 'low');
ecg2 = filtfilt(c,d,ecg2);
ecg3 = filtfilt(c,d,ecg3);

%IIR Notch filter to remove select 22Hz noise in E3
fn = 22;
Q = 450;
[bn, an] = iirnotch(fn/(fs/2), fn/(fs/2)/Q);
ecg3 = filter(bn,an,ecg3);


% Implementing Pan-Tompkins QRS complex detection

%Derivative - Amplifies high frequency spikes
z1 = zeros(1,length(ecg1));
z2 = zeros(1, length(ecg2));
z3 = zeros(1, length(ecg3));
for ii = 5:length(ecg2)
    z1(ii) = (1 / 8) * (2 *ecg1(ii) + ecg1(ii - 1) - ecg1(ii - 3) - 2 * ecg1(ii - 4));
    z2(ii) = (1 / 8) * (2 *ecg2(ii) + ecg2(ii - 1) - ecg2(ii - 3) - 2 * ecg2(ii - 4));
    z3(ii) = (1 / 8) * (2 *ecg3(ii) + ecg3(ii - 1) - ecg3(ii - 3) - 2 * ecg3(ii - 4));
end

% Squaring the signal increases amplitude differences
z1 = z1.*z1;
z2 = z2.*z2;
z3 = z3.*z3;

%Integration, here with a 3 length moving average window
% Attenuates spread out energy in P and T waves, doesn't affect concentrated peaks 

w1 = zeros(1,length(z1));
w2 = zeros(1, length(z2));
w3 = zeros(1, length(z3));
for ii = 11:length(z2)
    for jj = 0:3
        w1(ii) = w1(ii) + z1(ii - jj) / 2;
        w2(ii) = w2(ii) + z2(ii - jj) / 2;
        w3(ii) = w3(ii) + z3(ii - jj) / 2;
    end
end

%Plotting the time domain signals post filtering
figure;
subplot(3,1,1);
plot(axis/fs, ecg1(first_ten))
title('Filtered signal E1')
xlabel('Time(seconds)');
ylabel('Amplitude')
subplot(3,1,2);
plot(axis/fs, ecg2(first_ten))
title('Filtered signal E2');
xlabel('Time(seconds)');
ylabel('Amplitude')
subplot(3,1,3);
plot(axis/fs, ecg3(first_ten))
title('Filtered signal E3');
xlabel('Time(seconds)');
ylabel('Amplitude')

%Signals after applying Pan-Tompkins algorithm
figure;
subplot(3,1,1);
plot(axis/fs, w1(first_ten))
title('E1 after applying QRS detection')
xlabel('Time(seconds)');
ylabel('Amplitude')
subplot(3,1,2);
plot(axis/fs, w2(first_ten))
title('E2 after applying QRS Detection');
xlabel('Time(seconds)');
ylabel('Amplitude')
subplot(3,1,3);
plot(axis/fs, w3(first_ten))
title('E3 after applying QRS detection');
xlabel('Time(seconds)');
ylabel('Amplitude')

%Plotting frequency spectra post filtering
figure;
s1 = abs(fftshift(fft(ecg1,N)));
subplot(3,1,1)
plot(f_axis,s1);
title('Frequency Spectrum of E1 post filtering');
xlabel('Frequency(Hz)')
ylabel('Magnitude')
xlim([-64 64])
s2 = abs(fftshift(fft(ecg2,N)));
subplot(3,1,2)
plot(f_axis,s2);
title('Frequency Spectrum of E2 post filtering');
xlabel('Frequency(Hz)')
ylabel('Magnitude')
xlim([-64 64])
s3 = abs(fftshift(fft(ecg3,N)));
subplot(3,1,3)
plot(f_axis,s3);
title('Frequency Spectrum of E3 post filtering');
xlabel('Frequency(Hz)')
ylabel('Magnitude')
xlim([-64 64])

% A simple peak detection algorithm based on amplitude thresholding
% Here we set the amplitude threshold for a peak and the minimum distance
% between peaks

threshold = 0.1;
dist = 20;

% Initialize arrays to store the sample number of the detected peaks
peak_samplesE1 = [];
peak_samplesE2 = [];
peak_samplesE3 = [];

% Peak detection
for i = 2:length(w1)-1
    if w1(i) > w1(i-1) && w1(i) > w1(i+1) && w1(i) > threshold
        if ~isempty(peak_samplesE1)
            p1 = peak_samplesE1(end);
        else 
            p1 = 0;
        end
        if i - p1 <= dist
            if w1(i) > w1(p1)
                peak_samplesE1(end) = i;
            end
        else 
            peak_samplesE1 = [peak_samplesE1, i];
        end
    end
    if w2(i) > w2(i-1) && w2(i) > w2(i+1) && w2(i) > threshold
        if ~isempty(peak_samplesE2)
            p2 = peak_samplesE2(end);
        else 
            p2 = 0;
        end
        if i - p2 <= dist
            if w2(i) > w2(p2)
                peak_samplesE2(end) = i;
            end
        else 
            peak_samplesE2 = [peak_samplesE2, i];
        end
    end
    if w3(i) > w3(i-1) && w3(i) > w3(i+1) && w3(i) > threshold
        if ~isempty(peak_samplesE3)
            p3 = peak_samplesE3(end);
        else 
            p3 = 0;
        end
        if i - p3 <= dist
            if w3(i) > w3(p3)
                peak_samplesE3(end) = i;
            end
        else 
            peak_samplesE3 = [peak_samplesE3, i];
        end
    end
end

% Convert peak sample indices to real times (in seconds)
real_times1 = (peak_samplesE1 - 1) / fs;
real_times2 = (peak_samplesE2 - 1) / fs;
real_times3 = (peak_samplesE3 - 1) / fs;

% Initialize arrays to store heart rate data
heart_rate_bpm1 = [];
heart_rate_bpm2 = [];
heart_rate_bpm3 = [];

% Count number of peaks per minute for each signal - total length is 13 min
for ii = 1:13
    count1 = sum(real_times1 >= (ii-1)*60 & real_times1 < ii*60);
    heart_rate_bpm1 = [heart_rate_bpm1, count1];
    count2 = sum(real_times2 >= (ii-1)*60 & real_times2 < ii*60);
    heart_rate_bpm2 = [heart_rate_bpm2, count2];
    count3 = sum(real_times3 >= (ii-1)*60 & real_times3 < ii*60);
    heart_rate_bpm3 = [heart_rate_bpm3, count3];
end

% Plotting the heart rate data for each ECG signal

%Calculating the average over the entire signal
avg_bpm1 = mean(heart_rate_bpm1);
avg_bpm2 = mean(heart_rate_bpm2);
avg_bpm3 = mean(heart_rate_bpm3);

figure;

% BPM for E1
subplot(3,1,1);
stem(1:13, heart_rate_bpm1);
hold on;
plot(1:13, heart_rate_bpm1, 'r--', 'LineWidth', 1); 
yline(avg_bpm1, 'g-', 'LineWidth', 1); 
hold off;
ylim([90 125]);
title('Heart rate in Beats per Minute - Signal E1');
xlabel('Time(minutes)')
ylabel('HR')
legend(sprintf('Average Value: %.2f bpm', avg_bpm1), 'Location', 'northeast');

% BPM for E2
subplot(3,1,2);
stem(1:13, heart_rate_bpm2);
hold on;
plot(1:13, heart_rate_bpm2, 'r--', 'LineWidth', 1); 
yline(avg_bpm2, 'g-', 'LineWidth', 1); 
hold off;
ylim([90 125]);
title('Heart rate in Beats per Minute - Signal E2');
xlabel('Time(minutes)')
ylabel('HR')
legend( sprintf('Average Value: %.2f bpm', avg_bpm2), 'Location', 'northeast');

% BPM for E3
subplot(3,1,3);
stem(1:13, heart_rate_bpm3);
hold on;
plot(1:13, heart_rate_bpm3, 'r--', 'LineWidth', 1); 
yline(avg_bpm3, 'g-', 'LineWidth', 1);
hold off;
ylim([90 125]);
title('Heart rate in Beats per Minute - Signal E3');
xlabel('Time(minutes)')
ylabel('HR')
legend(sprintf('Average Value: %.2f bpm', avg_bpm3), 'Location', 'northeast');

