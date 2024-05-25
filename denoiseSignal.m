function denoised_signal = denoiseSignal(signal, Fs)
    % Design a low-pass filter
    fc = 50; % Cut-off frequency in Hz
    N = 100; % Filter order
    b = fir1(N, fc / (Fs / 2), 'low'); % Design a low-pass FIR filter

    % Apply the filter to the input signal
    denoised_signal = filtfilt(b, 1, signal);
end
