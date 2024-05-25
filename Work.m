clc
clear

%% part 1 heart rate interval %%
Fs = 250; % sampling frequency: ECG signal is sampled at a rate of 250 samples per second.
Fn = 50*pi; % Nyquist frequency: half of the sampling frequency. It is calculated as 50 times pi = 157 Hz. Fn = maximum frequency that can be represented in the sampled ECG signal without aliasing
T_new_s = 1;% store the duration of the ECG signal in seconds.  
T = 1/Fs; % sampling time: the time interval between successive samples.
x1 = 0:452323 - 1; %indices ranging from 0 to 452,322, with each index corresponding to a sample in the ECG signal of old individuals (The array is indexed from 0 to 452,322 to match the length of the ECG signal).
x2 = 0:452696 - 1; % For young
t1 = x1*T; % time vector for the ecg signal % multiplying each index in x1 by the sampling time T. It represents the time instances corresponding to each sample in the ECG signal of old individuals.
t2 = x2*T; % time vector

O = importdata('ECG_O.txt');
Y = importdata('ECG_Y.txt');

%a) remove 8V offset values
O = detrend(O); % remove linear trends from the ECG signal stored in variable O
Y = detrend(Y);
O = O - mean(O);% centers the signals around zero, which can be beneficial for certain analyses.
Y = Y - mean(Y);
N = 7; % order (7 coefficients in its transfer function) of the low-pass filter # higher-order filters provide steeper roll-off rates and better attenuation of out-of-band frequencies.
Fp = 100; % sets the passband frequency of the low-pass filter to 100 Hz.
Ap = 1; % sets the maximum passband ripple of the low-pass filter to 1 dB.
h = fdesign.lowpass('N,Fp,Ap',N,Fp,Ap,Fs); %
d = design(h,'cheby1'); % Chebyshev Type I filter: optimize between passband ripple and stopband attenuation.
xfilter = filter(d,O); % applies the designed filter d to the signal O
xfiltfilt = filtfilt(d.sosMatrix,d.ScaleValues,O); % filtfilt is a function used for zero-phase digital filtering. It applies a filter to a signal twice, once forwards and once backwards, to remove phase distortion. This ensures that the filtered signal has zero phase delay.
yfilter = filter(d,Y);
yfiltfilt = filtfilt(d.sosMatrix,d.ScaleValues,Y); % filters the signal Y in both forward and reverse directions to remove phase distortion.
O = xfiltfilt;
Y = yfiltfilt; % contain the filtered signal after zero-phase filtering.

% design a High Pass Filter
b=[-1/32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/32]; % numerator coefficients  
a=[1 -1]; % denominator coefficients
O=filter(b,a,O); % apply the designed high-pass filter defined by the numerator b and denominator a coefficients 
Y=filter(b,a,Y);

figure(1)
subplot(2,1,1)
plot(t1,O); % t1: time vector
title('ECG sinal for old');
xlabel('time(s)');
ylabel('Voltage (V)');
subplot(2,1,2)
plot(t2,Y);
title('ECG sinal for young');
xlabel('time(s)');
ylabel('Voltage (V)');
sgtitle('Preprocessed ECG singal');



% b)
figure(2)
subplot(2,2,1)
stem(t1,O); % plots the ECG signal against time (t1) using stems, which are vertical lines, for each data point in O.
xlim([0 5]);
title('ECG sinal for old');
xlabel('time(s)');
ylabel('Voltage (V)');
grid on
subplot(2,2,2)
stem(t2,Y); % stem plot is a type of plot that displays data as points along a line and extends a vertical line (a stem) from each data point down to a horizontal baseline. 
xlim([0 5]);
grid on
title('ECG sinal for young');
xlabel('time(s)');
ylabel('Voltage (V)');

subplot(2,2,3)
stem(t1,O);
xlim([6 10]);
grid on
title('ECG sinal for old');
xlabel('time(s)');
ylabel('Voltage (V)');
subplot(2,2,4)
stem(t2,Y);
xlim([6 10]);
title('ECG sinal for young');
xlabel('time(s)');
ylabel('Voltage (V)');
grid on

% b) amplitude spectrum : the absolute value of voltage = amplitude (=magnitude or strength of the signal: distance from the baseline to the peak)
figure(3)
subplot(2,1,1)
dt1 = diff(t1); % Computes the differences between consecutive elements of the time vector t1, which gives the time intervals between samples.
nfft1 = length(O); % Determines the length of the signal O, which is used for the number of points in the FFT (Fast Fourier Transform).
A1 = fft(O,nfft1); % Computes the FFT of the signal O with nfft1 points.
A1 = A1(1:nfft1/2); % Extracts the first half of the FFT result, as the second half is the mirror image due to the symmetry of the FFT output for real-valued input.
m1 = abs(A1); % Computes the magnitude spectrum by taking the absolute values of the FFT result.
f1 = (0:nfft1/2-1)*Fs/nfft1; % Computes the frequency axis for plotting. The frequencies are calculated from 0 to the Nyquist frequency (Fs/2) with a resolution of Fs/nfft1.
plot(f1,m1); % This plot provides insight into the frequency content of the ECG signal, showing how much of each frequency component is present in the signal.
title('Amplitude spectrum of ECG sinal for old');
xlabel('frequency (Hz)');
ylabel('|V(f)|');
grid on
%to reconstruct a continuous-time signal from its samples, fs must be at least twice the highest frequency component (Nyquist frequency) present in the signal. So, in FFT, it's customary to display the frequency axis from 0 to Fs/2 to cover the entire frequency range of interest while avoiding aliasing.
subplot(2,1,2)
dt2 = diff(t2);
nfft2 = length(Y);
A2 = fft(Y,nfft2);
A2 = A2(1:nfft2/2);
m2 = abs(A2);
f2 = (0:nfft2/2-1)*Fs/nfft2;
plot(f2,m2);
title('Amplitude spectrum of ECG sinal for young');
xlabel('freqeuncy (Hz)');
ylabel('|V(f)|');
grid on

% figure(3)
% subplot(2,1,1)
% odft = fft(O);
% odft = odft(1:length(O)/2+1);
% freq = 0:Fs/length(O): Fs/2;
% plot(freq,abs(odft));
% subplot(2,1,2)
% ydft = fft(Y);
% ydft = ydft(1:length(Y)/2+1);
% freq1 = 0:Fs/length(Y): Fs/2;
% plot(freq1,abs(ydft));

% finding the time difference between consecutive R-peaks in ecg signal
% c) NN interval by hand calculation & d) plot RR
int_O = 1.012; % s for old % the average NN (normal to normal) interval duration: time duration between consecutive R-peaks in ecg, which is a measure of the heart rate variability
int_Y = 0.86; % s for young %HR=60/NN-interval=60/0.86=69.77bpm=average heart rate 
%300 represents maximum expected distance between consecutive peaks in samples: limit search space for peak detection and improve efficiency of the algorithm. Peaks that are farther apart than this threshold will not be considered as part of the same RR interval.
%t_QRS contains time points corresponding to the R-peaks detected in the signal. Each element of t_QRS represents the time point (in seconds) where an R-peak occurs.
%QRS contains the amplitude values of the R-peaks detected in the signal. Each element of QRS corresponds to the voltage (or amplitude) of the R-peak at the corresponding time point in t_QRS.
%RR contains the durations of the RR intervals, which are the time intervals between consecutive R-peaks. Each element of RR represents duration (in seconds) of an RR interval.
%find_RR detects R-peaks in the input signal O using the provided parameters and returns the time points of the R-peaks (t_QRS_O), their corresponding amplitudes (QRS_O), and the durations of the RR intervals (RR_O).

[t_QRS_O,QRS_O,RR_O] = find_RR(O,t1,0,300);
[t_QRS_Y,QRS_Y,RR_Y] = find_RR(Y,t2,0,300);

figure(4)
subplot(2,1,1)
plot(t_QRS_O,QRS_O); %plotting the amplitude of the R-peaks against time. This visualization can provide insights into the amplitude variations of the R-peaks over time.
title('R peaks indentification represented by line for old');
xlabel('time(s)');
ylabel('Voltage (V)');
subplot(2,1,2)
plot(t_QRS_Y,QRS_Y);
title('R peaks indentification represented by line for young');
xlabel('time(s)');
ylabel('Voltage (V)');

figure(5) % variability in the NN intervals over time % more variable tachogram suggests greater heart rate variability, which is often associated with better cardiovascular health and autonomic function.
subplot(2,1,1)
plot(RR_O); % contains the durations of the NN intervals, representing the time intervals between consecutive R-peaks (heartbeats) in the ECG signal O. Each element of RR_O corresponds to the duration of an NN interval.
title('tachogram of ECG signal for old') % tachogram: graphical representation of heart rate variability (HRV) over time % displays the NN (normal-to-normal) intervals, which are the time intervals between successive heartbeats, plotted against the beat number or time.
xlabel('beat number');
ylabel('NN interval for Old');
subplot(2,1,2)
plot(RR_Y);
title('tachogram of ECG signal for young')
xlabel('beat number');
ylabel('NN interval for Young');
sgtitle('tachogram');

% e) mean plot & range plot % identifies R-peaks and calculates the NN interval durations for segments of the ECG signal, and computes the maximum, minimum, mean, and range of these NN intervals using the Mean_range function. These statistical measures are typically used to assess HRV and can provide valuable insights into cardiac function.
% old % Each segment covers a time interval of 300 seconds (or 5 minutes).Each segment is analyzed to examine how the statistical properties of the NN intervals (like mean) vary across different time periods. This approach allows for the identification of any temporal patterns or changes in HRV within the signal. By analyzing multiple segments, researchers can gain a more comprehensive understanding of the dynamics of the cardiovascular system over time.
%t_QRS_1: Time points of the R-peaks. QRS_1: Amplitudes of the R-peaks. RR_1: Durations of the NN intervals (time intervals between consecutive R-peaks).
%max_1: Maximum NN interval duration in segment 1. min_1: Minimum NN interval duration in segment 1. M_1: Mean (average) of the NN interval durations in segment 1. range_1: Range of the NN interval durations in segment 1 (the difference between the maximum and minimum NN intervals).
[t_QRS_1,QRS_1,RR_1] = find_RR(O,t1,0,300);% segment 1 % detect R-peaks in the ECG signal O within the time interval [0, 300] seconds.
[max_1,min_1,M_1,range_1] = Mean_range(RR_1);

[t_QRS_2,QRS_2,RR_2] = find_RR(O,t1,300,600);% segment 2 
[max_2,min_2,M_2,range_2] = Mean_range(RR_2);

[t_QRS_3,QRS_3,RR_3] = find_RR(O,t1,600,900);% segment 3
[max_3,min_3,M_3,range_3] = Mean_range(RR_3);

[t_QRS_4,QRS_4,RR_4] = find_RR(O,t1,900,1200);% segment 4
[max_4,min_4,M_4,range_4] = Mean_range(RR_4);

[t_QRS_5,QRS_5,RR_5] = find_RR(O,t1,1200,1500);% segment 5
[max_5,min_5,M_5,range_5] = Mean_range(RR_5);

[t_QRS_6,QRS_6,RR_6] = find_RR(O,t1,1500,t1(end));% segment 6
[max_6,min_6,M_6,range_6] = Mean_range(RR_6);

mean_RRo = 1:6; % segment #s
mean_RRo = ([M_1 M_2 M_3 M_4 M_5 M_6])'; % creates a column vector by concatenating the mean values (M_1, M_2, ..., M_6) computed for each segment.
max_RRo = 1:6;
max_RRo = ([max_1 max_2 max_3 max_4 max_5 max_6])'; % values are in [], and the apostrophe ' transposes the vector to make it a column vector.
min_RRo = 1:6;
min_RRo = ([min_1 min_2 min_3 min_4 min_5 min_6])';
RR_rangeo = 1:6;
RR_rangeo =([range_1 range_2 range_3 range_4 range_5 range_6])';
segmento = 1:6;
segmento =([1 2 3 4 5 6])';
intervalo = (["0-300","300-600","600-900","900-1200","1200-1500","1500-end"])'; % containing strings representing the time intervals for each segment.
tablevalue_old = table(segmento,intervalo,RR_rangeo,min_RRo,max_RRo,mean_RRo);
tablevalue_old(1:6,:); %displays the first 6 rows of the table tablevalue_old. It selects rows from 1 to 6 and all columns (:) for display.

% young %
[t_QRS_1y,QRS_1y,RR_1y] = find_RR(Y,t2,0,300);% segment 1
[max_1y,min_1y,M_1y,range_1y] = Mean_range(RR_1y);

[t_QRS_2y,QRS_2y,RR_2y] = find_RR(Y,t2,300,600);% segment 2
[max_2y,min_2y,M_2y,range_2y] = Mean_range(RR_2y);

[t_QRS_3y,QRS_3y,RR_3y] = find_RR(Y,t2,600,900);% segment 3
[max_3y,min_3y,M_3y,range_3y] = Mean_range(RR_3y);

[t_QRS_4y,QRS_4y,RR_4y] = find_RR(Y,t2,900,1200);% segment 4
[max_4y,min_4y,M_4y,range_4y] = Mean_range(RR_4y);

[t_QRS_5y,QRS_5y,RR_5y] = find_RR(Y,t2,1200,1500);% segment 5
[max_5y,min_5y,M_5y,range_5y] = Mean_range(RR_5y);

[t_QRS_6y,QRS_6y,RR_6y] = find_RR(Y,t2,1500,t2(end));% segment 6
[max_6y,min_6y,M_6y,range_6y] = Mean_range(RR_6y);

mean_RR = 1:6;
mean_RRy = ([M_1y M_2y M_3y M_4y M_5y M_6y])';
max_RRy = 1:6;
max_RRy = ([max_1y max_2y max_3y max_4y max_5y max_6y])';
min_RRy = 1:6;
min_RRy = ([min_1y min_2y min_3y min_4y min_5y min_6y])';
RR_rangey = 1:6;
RR_rangey =([range_1y range_2y range_3y range_4y range_5y range_6y])';
segmenty = 1:6;
segmenty =([1 2 3 4 5 6])';
intervaly = (["0-300","300-600","600-900","900-1200","1200-1500","1500-end"])';
tablevalue_young = table(segmenty,intervaly,RR_rangey,min_RRy,max_RRy,mean_RRy);
tablevalue_young(1:6,:);

figure(7)
subplot(2,2,1)
stem(segmento,mean_RRo);
title('mean value  for old');
xlabel('segment number');
ylabel('mean');
subplot(2,2,2)
stem(segmento,RR_rangeo);
title('range value for old');
xlabel('segment number');
ylabel('range');
subplot(2,2,3)
stem(segmenty,mean_RRy);
title('mean value for young');
xlabel('segment number');
ylabel('mean');
subplot(2,2,4)
stem(segmenty,RR_rangey);
title('range value for young');
xlabel('segment number');
ylabel('range');
%% part 2 time domain anaylsis %%
% old %
so1 = std(RR_1);
so2 = std(RR_2);
so3 = std(RR_3);
so4 = std(RR_4);
so5 = std(RR_5);
so6 = std(RR_6);
so = ([so1 so2 so3 so4 so5 so6])';

sy1 = std(RR_1y);
sy2 = std(RR_2y);
sy3 = std(RR_3y);
sy4 = std(RR_1y);
sy5 = std(RR_5y);
sy6 = std(RR_6y);
sy = ([sy1 sy2 sy3 sy4 sy5 sy6])';

figure(8)
subplot(2,1,1)
stem(segmento,so);
title('standard deviation for old');
xlabel('segment number');
ylabel('SDRR');

subplot(2,1,2)
stem(segmenty,sy);
title('standard deviation for young');
xlabel('segment number');
ylabel('SDRR');
%% part 3 frequency domain analysis %%
% old %
[pow_1,f_1,f_new_1,TOT_LF_1,TOT_HF_1,ratio_1] = frequencyanalysis(RR_1,T_new_s)
figure(9)
subplot(2,1,1)
plot(f_1,pow_1);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_1,pow_1);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 1 time recorded 100-300 s');

[pow_2,f_2,f_new_2,TOT_LF_2,TOT_HF_2,ratio_2] = frequencyanalysis(RR_2,T_new_s)
figure(10)
subplot(2,1,1)
plot(f_2,pow_2);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_2,pow_2);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 2 time recorded 300-600 s');

[pow_3,f_3,f_new_3,TOT_LF_3,TOT_HF_3,ratio_3] = frequencyanalysis(RR_3,T_new_s)
figure(11)
subplot(2,1,1)
plot(f_3,pow_3);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_3,pow_3);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 3 time recorded 600-900 s');

[pow_4,f_4,f_new_4,TOT_LF_4,TOT_HF_4,ratio_4] = frequencyanalysis(RR_4,T_new_s)
figure(12)
subplot(2,1,1)
plot(f_4,pow_4);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_4,pow_4);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 4 time recorded 900-1200 s');

[pow_5,f_5,f_new_5,TOT_LF_5,TOT_HF_5,ratio_5] = frequencyanalysis(RR_5,T_new_s)
figure(13)
subplot(2,1,1)
plot(f_5,pow_5);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_5,pow_5);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 5 time recorded 1200-1500 s');

[pow_6,f_6,f_new_6,TOT_LF_6,TOT_HF_6,ratio_6] = frequencyanalysis(RR_6,T_new_s)
figure(14)
subplot(2,1,1)
plot(f_6,pow_6);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_6,pow_6);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 6 time recorded 1500-1800 s');

% young %

[pow_1y,f_1y,f_new_1y,TOT_LF_1y,TOT_HF_1y,ratio_1y] = frequencyanalysis(RR_1y,T_new_s)
figure(15)
subplot(2,1,1)
plot(f_1y,pow_1y);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_1y,pow_1y);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 1 time recorded 100-300 s');

[pow_2y,f_2y,f_new_2y,TOT_LF_2y,TOT_HF_2y,ratio_2y] = frequencyanalysis(RR_2y,T_new_s)
figure(16)
subplot(2,1,1)
plot(f_2y,pow_2y);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_2y,pow_2y);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 2 time recorded 300-600 s');

[pow_3y,f_3y,f_new_3y,TOT_LF_3y,TOT_HF_3y,ratio_3y] = frequencyanalysis(RR_3y,T_new_s)
figure(17)
subplot(2,1,1)
plot(f_3y,pow_3y);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_3y,pow_3y);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 3 time recorded 600-900 s');

[pow_4y,f_4y,f_new_4y,TOT_LF_4y,TOT_HF_4y,ratio_4y] = frequencyanalysis(RR_4y,T_new_s)
figure(18)
subplot(2,1,1)
plot(f_4y,pow_4y);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_4y,pow_4y);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 4 time recorded 900-1200 s');

[pow_5y,f_5y,f_new_5y,TOT_LF_5y,TOT_HF_5y,ratio_5y] = frequencyanalysis(RR_5y,T_new_s)
figure(19)
subplot(2,1,1)
plot(f_5y,pow_5y);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_5y,pow_5y);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 5 time recorded 1200-1500 s');

[pow_6y,f_6y,f_new_6y,TOT_LF_6y,TOT_HF_6y,ratio_6y] = frequencyanalysis(RR_6y,T_new_s)
figure(20)
subplot(2,1,1)
plot(f_6y,pow_6y);
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude in dB')
subplot(2,1,2)
plot(f_new_6y,pow_6y);
title('Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude in dB');
sgtitle('segment 6 time recorded 1500-1800 s');


%% part 4 additional analysis %%
%% part 4 additional analysis %%
% old %
[xm_1, xp_1, SD1_1, SD2_1, AFib_1, RMSSD_1] = Pointcareanalysis(RR_1);
figure(21)
plot(xp_1, xm_1, '.'); % Poincaré plot
title('Poincare Plot for segmentation 1 of old');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_1), 'Units', 'normalized'); % Display AFib status
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_1), 'Units', 'normalized'); % Display RMSSD

[xm_2, xp_2, SD1_2, SD2_2, AFib_2, RMSSD_2] = Pointcareanalysis(RR_2);
figure(22)
plot(xp_2, xm_2, '.');
title('Poincare Plot for segmentation 2 of old');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_2), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_2), 'Units', 'normalized');

[xm_3, xp_3, SD1_3, SD2_3, AFib_3, RMSSD_3] = Pointcareanalysis(RR_3);
figure(23)
plot(xp_3, xm_3, '.');
title('Poincare Plot for segmentation 3 of old');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_3), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_3), 'Units', 'normalized');

[xm_4, xp_4, SD1_4, SD2_4, AFib_4, RMSSD_4] = Pointcareanalysis(RR_4);
figure(24)
plot(xp_4, xm_4, '.');
title('Poincare Plot for segmentation 4 of old');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_4), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_4), 'Units', 'normalized');

[xm_5, xp_5, SD1_5, SD2_5, AFib_5, RMSSD_5] = Pointcareanalysis(RR_5);
figure(25)
plot(xp_5, xm_5, '.');
title('Poincare Plot for segmentation 5 of old');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_5), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_5), 'Units', 'normalized');

[xm_6, xp_6, SD1_6, SD2_6, AFib_6, RMSSD_6] = Pointcareanalysis(RR_6);
figure(26)
plot(xp_6, xm_6, '.');
title('Poincare Plot for segmentation 6 of old');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_6), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_6), 'Units', 'normalized');

% young %

[xm_1y, xp_1y, SD1_1y, SD2_1y, AFib_1y, RMSSD_1y] = Pointcareanalysis(RR_1y);
figure(27)
plot(xp_1y, xm_1y, '.');
title('Poincare Plot for segmentation 1 of young');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_1y), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_1y), 'Units', 'normalized');

[xm_2y, xp_2y, SD1_2y, SD2_2y, AFib_2y, RMSSD_2y] = Pointcareanalysis(RR_2y);
figure(28)
plot(xp_2y, xm_2y, '.');
title('Poincare Plot for segmentation 2 of young');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_2y), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_2y), 'Units', 'normalized');

[xm_3y, xp_3y, SD1_3y, SD2_3y, AFib_3y, RMSSD_3y] = Pointcareanalysis(RR_3y);
figure(29)
plot(xp_3y, xm_3y, '.');
title('Poincare Plot for segmentation 3 of young');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_3y), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_3y), 'Units', 'normalized');

[xm_4y, xp_4y, SD1_4y, SD2_4y, AFib_4y, RMSSD_4y] = Pointcareanalysis(RR_4y);
figure(30)
plot(xp_4y, xm_4y, '.');
title('Poincare Plot for segmentation 4 of young');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_4y), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_4y), 'Units', 'normalized');

[xm_5y, xp_5y, SD1_5y, SD2_5y, AFib_5y, RMSSD_5y] = Pointcareanalysis(RR_5y);
figure(31)
plot(xp_5y, xm_5y, '.');
title('Poincare Plot for segmentation 5 of young');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_5y), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_5y), 'Units', 'normalized');

[xm_6y, xp_6y, SD1_6y, SD2_6y, AFib_6y, RMSSD_6y] = Pointcareanalysis(RR_6y);
figure(32)
plot(xp_6y, xm_6y, '.');
title('Poincare Plot for segmentation 6 of young');
xlabel('RR interval(i)');
ylabel('RR interval(i+1)');
text(0.1, 0.9, sprintf('AFib: %d', AFib_6y), 'Units', 'normalized');
text(0.1, 0.8, sprintf('RMSSD: %.2f', RMSSD_6y), 'Units', 'normalized');


%% part 5 additional analysis %%
%% Continuous Wavelet Transform (CWT):

% Call wavelet_analysis function for segment 1
[cwt_coeffs_old, frequencies_old, scales_old] = wavelet_analysis(RR_1, T_new_s); %computes the CWT coefficients and the corresponding frequencies and scales.
[cwt_coeffs_young, frequencies_young, scales_young] = wavelet_analysis(RR_1y, T_new_s);

% Plots show how the frequency content of the signal varies with time. Darker regions in the plot represent higher coefficients or greater energy at corresponding time-frequency locations, indicating significant frequency components in the signal.
% Plot continous wavelet scalogram for segment 1 - old  %color represents the magnitude of the wavelet coefficients.
figure(33)
subplot(2,1,1)
imagesc(1:length(RR_1), frequencies_old, abs(cwt_coeffs_old)); %computes the absolute values of the CWT coefficients to visualize the magnitude.
title('Wavelet Scalogram - Segment 1 (Old)');
xlabel('Time');
ylabel('Frequency (Hz)');
colorbar;

% Plot wavelet scalogram for segment 1 - young
subplot(2,1,2)
imagesc(1:length(RR_1y), frequencies_young, abs(cwt_coeffs_young));
title('Wavelet Scalogram - Segment 1 (Young)');
xlabel('Time');
ylabel('Frequency (Hz)');
colorbar;


%% Discrete Wavelet Transform (DWT):
% Call dwt_analysis function for segment 1 - old
[cA_old, cD_old] = dwt_analysis(RR_1,1); % level 1 (if no level: (RR_1) in both codes)

% Call dwt_analysis function for segment 1 - young
[cA_young, cD_young] = dwt_analysis(RR_1y,1);

% Plot discrete wavelet coefficients for segment 1 - old
figure(34)
subplot(2,1,1)
plot(cA_old);
hold on;
plot(cD_old);
hold off;
title('Discrete Wavelet Coefficients - Segment 1 (Old)');
legend('Approximation Coefficients (cA)', 'Detail Coefficients (cD)');
xlabel('Sample Index');
ylabel('Coefficient Value');

% Plot discrete wavelet coefficients for segment 1 - young
subplot(2,1,2)
plot(cA_young);
hold on;
plot(cD_young);
hold off;
title('Discrete Wavelet Coefficients - Segment 1 (Young)');
legend('Approximation Coefficients (cA)', 'Detail Coefficients (cD)');
xlabel('Sample Index');
ylabel('Coefficient Value');


% level 10 of decomposition:
%% Discrete Wavelet Transform (DWT):
% Call dwt_analysis function for segment 1 - old
[cA_old, cD_old] = dwt_analysis(RR_1,10);

% Call dwt_analysis function for segment 1 - young
[cA_young, cD_young] = dwt_analysis(RR_1y,10);

% Plot discrete wavelet coefficients for segment 1 - old
figure(35)
subplot(2,1,1)
plot(cA_old);
hold on;
plot(cD_old);
hold off;
title('Discrete Wavelet Coefficients - Segment 1 (Old)');
legend('Approximation Coefficients (cA)', 'Detail Coefficients (cD)');
xlabel('Sample Index');
ylabel('Coefficient Value');

% Plot discrete wavelet coefficients for segment 1 - young
subplot(2,1,2)
plot(cA_young);
hold on;
plot(cD_young);
hold off;
title('Discrete Wavelet Coefficients - Segment 1 (Young)');
legend('Approximation Coefficients (cA)', 'Detail Coefficients (cD)');
xlabel('Sample Index');
ylabel('Coefficient Value');
%{
%% Part 6: Additional Analysis % Atrial Fibrillation % Added in Poincare function
% Preprocess ECG signals if necessary (e.g., denoising, filtering)
O_denoised = denoiseSignal(O, Fs);
Y_denoised = denoiseSignal(Y, Fs);

% Plot the denoised ECG signals for visual inspection
figure;
subplot(2,1,1);
plot((0:length(O)-1)/Fs, O_denoised);
title('Old Individual ECG Signal (Denoised)');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot((0:length(Y)-1)/Fs, Y_denoised);
title('Young Individual ECG Signal (Denoised)');
xlabel('Time (s)');
ylabel('Amplitude');

% Call the AFib detection function for old individuals' ECG signal
AFib_O = detectAFib(O_denoised, Fs);
disp(['AFib detected in old individual: ', num2str(AFib_O)]);
% Call the AFib detection function for young individuals' ECG signal
AFib_Y = detectAFib(Y_denoised, Fs);
disp(['AFib detected in young individual: ', num2str(AFib_Y)]);
%}

%% Part 6: Additional Analysis % SVM Trainer not work
%accuracy = SVM_Trainer(O, Y); %SVM classifier to distinguish between ECG signals of old and young individuals and evaluate its performance using a testing dataset.
