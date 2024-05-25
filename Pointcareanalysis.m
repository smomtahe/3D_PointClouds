function [xm, xp, SD1, SD2, AFib, RMSSD] = Pointcareanalysis(ECG)
    % Perform a Poincaré plot analysis on a given time series data, such as RR intervals from an ECG signal.

    % Shift RR intervals
    xp = ECG; % RR intervals shifted one sample to the right.
    xp(end) = []; 
    xm = ECG; % RR intervals shifted one sample to the left.
    xm(1) = [];
    
    % Calculate SD1
    SD1 = std(xp - xm) / sqrt(2); % Short-term variability
    
    % Calculate SD2
    SD2 = sqrt((2 * std(ECG) ^ 2) - SD1 ^ 2); % Long-term variability

    % AFib detection
    % Implement your AFib detection algorithm here
    % Example: Check if SD1 exceeds a certain threshold
    threshold_AFib_SD1 = 0.1; % Adjust as needed
    if SD1 > threshold_AFib_SD1
        AFib = 1; % AFib detected
    else
        AFib = 0; % No AFib detected
    end
    
    % Calculate RMSSD
    RR_diff = diff(ECG);
    RMSSD = sqrt(mean(RR_diff .^ 2));
end


%{
Pointcareanalysis calculates SD1 and SD2, and returns them along with the shifted RR intervals for further analysis. 
The Poincaré plot is often used in HRV analysis to visualize and quantify the dynamic behavior of RR interval time series.

Output:
xm: RR intervals shifted one sample to the left.
xp: RR intervals shifted one sample to the right.
SD1: Standard deviation representing short-term variability in the Poincaré plot.
SD2: Standard deviation representing long-term variability in the Poincaré plot.
%}