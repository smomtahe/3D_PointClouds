function [pow,f,f_new,TOT_LF,TOT_HF,ratio] = frequencyanalysis(sig,T_new_s)
NFFT = 2^nextpow2(length(sig)); %Determines the next power of 2 greater than or equal to the length of the input signal sig and assigns it to NFFT. This is used for zero-padding the signal to improve frequency resolution in the FFT.
Y = fft(sig,NFFT)/length(sig); %FFT of the input signal sig with length NFFT %normalizes the FFT result by dividing it by the length of the signal.
f = 250/2*linspace(0,1,NFFT/2 + 1); %Creates the frequency axis (f) corresponding to the FFT result. It spans from 0 to the Nyquist frequency (half of the sampling frequency, 250 Hz in this case) with a resolution determined by NFFT.

pow = 2*10*log10(abs(Y(1:NFFT/2+1))); % PSD: takes the absolute value of the FFT result, squares it (to obtain power), and then converts it to decibels (dB) using a factor of 2 and a base 10 logarithm. Only the first half of the FFT result is considered, as the second half is the mirror image due to the symmetry of the FFT output for real-valued input.

%Frequency Axis Adjustment
f_new_s = 1/T_new_s; % Calculates the sampling frequency (f_new_s) based on the sampling period T_new_s.
f_new = f_new_s/2*linspace(0,1,NFFT/2+1); % Adjusts the frequency axis (f_new) to match the correct frequency range based on the new sampling frequency.

%Frequency Band Power Calculation:
LF = [0.04 0.15]; % LF region
HF = [0.15 0.4]; % HF region
 %Determine the indices of frequencies within the LF and HF regions.
 iLF = (pow>=LF(1)) & (pow<=LF(2)); % index value in LF region
 iHF = (pow>=HF(1)) & (pow<=HF(2)); % index value in HF region

%Compute the total power within the LF and HF regions using numerical integration (trapezoidal rule).
TOT_LF = abs(trapz(LF(1):LF(2),pow)); % calculates the total power within the HF band of the PSD estimate (pow). PSD represents the distribution of power intensity across different frequencies in the signal.%HF(1) represents the lower bound, and HF(2) represents the upper bound of the HF band.
TOT_HF = abs(trapz(HF(1):HF(2),pow)); %trapz computes the integral of data using the trapezoidal method. It estimates the integral of a function over a specified range by approximating the area under the curve using trapezoids.

ratio = TOT_LF/TOT_HF;  %alculate the ratio of LF power to HF power.
end
%Output:
%pow: Power spectral density estimate.
%f: Frequency axis.
%f_new: Adjusted frequency axis based on the new sampling frequency.
%TOT_LF and TOT_HF: Total power within LF and HF regions, respectively.
%ratio: Ratio of LF power to HF power.

% power spectral density estimation and calculation of power in specific frequency bands (Low Frequency, LF, and High Frequency, HF)
% NFFT = 2^nextpow2(length(sig)) is responsible for determining the length of the FFT (Fast Fourier Transform) window, assigned to the variable NFFT. 
% The function nextpow2 computes the next power of 2 greater than or equal to the length of the input signal sig, ensuring that NFFT is a power of 2. 
% This technique is commonly used in signal processing to optimize the efficiency of the FFT algorithm. By using a length that is a power of 2, the FFT computation becomes more efficient. 
% Additionally, this process involves zero-padding the signal, meaning that additional zeros are appended to the end of the signal to reach the desired length. 
% Zero-padding helps to improve the frequency resolution of the FFT output, allowing for more precise analysis of the signal's frequency content. 







