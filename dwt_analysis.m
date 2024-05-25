function [cA, cD] = dwt_analysis(signal, level)
    % Perform Discrete Wavelet Transform (DWT) for signal analysis

    % Define wavelet parameters
    wavelet = 'db4'; % Daubechies wavelet with 4 coefficients
    % Compute Discrete Wavelet Transform (DWT)
    [cA, cD] = dwt(signal, wavelet, level);
end

%{ 
The plot of DWT coefficients visually represents the decomposition of the signal into different frequency bands. 
When you plot the approximation coefficients (`cA`) and detail coefficients (`cD`), you're essentially visualizing how the signal is decomposed across different frequency scales. Here's what each component typically represents:

1. **Approximation Coefficients (`cA`)**:
   - These coefficients represent the low-frequency components of the signal after decomposition.
   - Plotting `cA` shows the overall trend or the "smooth" part of the signal.
   - This component captures the coarse details of the signal.
2. **Detail Coefficients (`cD`)**:
   - These coefficients represent the high-frequency components or details of the signal after decomposition.
   - Plotting `cD` shows the fluctuations or the "rough" part of the signal.
   - Each level of `cD` represents details at different scales or frequencies, with higher levels capturing finer details.
In summary, the plot of DWT coefficients provides insight into how the signal's frequency content is distributed across different scales or levels of detail. It helps in analyzing both the overall trend and the finer variations within the signal.

cA has a higher magnitude compared to cD. This is often observed due to the nature of the decomposition process and the properties of the signal being analyzed. Here are a few reasons why the magnitude of `cA` is generally higher:
1. **Smoothing Effect**: The `cA` coefficients represent the low-frequency components of the signal, capturing the overall trend or smooth part of the signal. In many cases, the low-frequency components tend to have higher magnitudes compared to the high-frequency components, resulting in larger `cA` coefficients.
2. **Energy Concentration**: The majority of the signal's energy is often concentrated in the lower-frequency components, which are captured by `cA`. Therefore, `cA` tends to have larger magnitudes since it contains a significant portion of the signal's energy.
3. **Downsampling**: In the DWT process, downsampling is performed after each level of decomposition, reducing the number of coefficients in each subsequent level. This downsampling can lead to a perceived increase in magnitude for `cA` relative to `cD`.
4. **Wavelet Basis**: The choice of wavelet basis can also influence the magnitudes of `cA` and `cD`. Some wavelets are designed to preserve energy or have specific properties that affect the decomposition coefficients' magnitudes.
%}