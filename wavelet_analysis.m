function [cwt_coeffs, frequencies, scales] = wavelet_analysis(signal, sampling_frequency)  
    % Perform continuous wavelet transform (CWT) for wavelet analysis

    % Define wavelet parameters
    wavelet = 'morl'; % Morlet wavelet
    scales = 1:1:128; % Wavelet scales
    frequencies = scal2frq(scales, wavelet, 1/sampling_frequency); % Corresponding frequencies

    % Compute Continuous Wavelet Transform (CWT)
    cwt_coeffs = cwt(signal, scales, wavelet);

end
%{
Continuous Wavelet Transform (CWT) operates by convolving the signal with a scaled and translated version of a continuous wavelet function.
It provides a continuous-time representation of the signal's frequency content across different scales.
The CWT can capture both frequency and time information simultaneously.
The output of the CWT is a continuous 2D plot called a scalogram, which shows how the frequency content of the signal changes over time.

A scalogram is a visualization of the wavelet transform coefficients across time and frequency. It's essentially a two-dimensional plot that represents the intensity or magnitude of wavelet coefficients as a function of both time and frequency.
In the context of wavelet analysis, the scalogram provides a way to visually inspect how the frequency content of a signal changes over time. The x-axis typically represents time, the y-axis represents frequency, and the color or intensity of each point in the plot represents the magnitude or strength of the wavelet coefficients at that particular time and frequency.
Scalograms are useful for analyzing signals with time-varying frequency content, such as non-stationary signals or signals containing transient events. They are commonly used in various fields including signal processing, time series analysis, and image processing.
%}