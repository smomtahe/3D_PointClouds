function [t_QRS,QRS,RR] = find_RR(sig,t1,L1,L2)
%Fixing threshold to find QRS complex that is peak point
Fs = 250;
T = 1/Fs;
range = find((t1>=L1) & (t1 <= L2));
t1_new = t1(range);
O_new = sig(range);

th = abs(max(O_new));
th = 0.5*th;
c = 0; % counter: keep track of the number of peaks (QRS complexes) detected in the signal.
% peaks are defined as pulses
for i = 1:length(O_new)
    if O_new(i) > th %if the amplitude of the signal at index i (current position in the loop) is greater than the threshold th: If the amplitude exceeds the threshold, it's considered a peak (QRS complex).
        QRS(i) = 1; %presence of a peak (or QRS complex) at index i.
        c = c+1;
        i = i+15;%skip a certain number of samples after detecting a peak: By skipping samples after detecting a peak, you ensure that subsequent peaks are not mistakenly detected as part of the same QRS complex.
    else
        QRS(i) = 0;%no peak at index i in the signal.
    end
end


t_QRS = (1:length(QRS)) * T %This line creates a time vector t_QRS corresponding to the detected peaks.

% figure(4)
% plot(t_QRS,QRS);
% xlim([0 5]);

%These lines remove zero points from the time vector t_QRS where no peaks were detected.
t_QRS_Real = t_QRS;
t_QRS_Real(QRS == 0) = []; % no zero points
%t_QRS_Real(2:2:end) = [];

b_num = 0;
%Calculates RR intervals based on the time differences between consecutive peaks.
for j = 1:length(t_QRS_Real) - 1
        RR(j) = t_QRS_Real(j+1) - t_QRS_Real(j);
        b_num = j;
end

 %t_QRS(RR < 0.03) = [];
 RR(RR < 0.03) = []; %This line removes RR intervals shorter than a certain threshold (0.03 seconds).
 b_num = 292; %This line assigns a fixed beat number (292 in this case), but it seems unused and might be a remnant of earlier code.
 
end

