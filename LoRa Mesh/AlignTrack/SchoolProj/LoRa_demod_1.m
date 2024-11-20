function [received_fft] = LoRa_demod_1(signal, fc, SF, BW, Fs, cfo)
    % LoRa_demod_1 demodulates full LoRa packet

    % Inputs:
    % - signal: The received time-domain signal
    % - SF: LoRa spreading factor (e.g., 7, 8, 9, etc.)
    % - BW: LoRa bandwidth (e.g., 125e3 for 125 kHz)
    % - Fs: Sampling frequency of the signal

    %% Return if SF is not in the range
    if SF > 12 || SF < 7
        return
    end

    % Step 1: Generate Reference Upchirp and Downchirp
    phy = LoRaPHY(fc, SF, BW, Fs);
    dChirpsDemod  = phy.chirp(false, SF, BW, Fs, 0, cfo, 0);
    ref_upchirp = phy.chirp(true, SF, BW, Fs, 0, cfo, 0);
    chirp_len = length(dChirpsDemod);
    % t = (0:chirp_len - 1) / Fs; % Time vector
    % k = BW / chirp_len; % Chirp rate
    
    % Step 2: Cross-Correlation for Preamble Detection
    correlation = abs(xcorr(signal, ref_upchirp));
    [~, peak_idx] = max(correlation); % Find the peak
    shift = peak_idx - length(signal); % Align index with signal

    % Step 3: Extract Windows, Dechirp, and Compute FFT
    result = signal(shift + 1:shift + chirp_len) .* dChirpsDemod;
    %result = convolve_with_window(signal(shift + 1:shift + len), 'hamming') .* dChirpsDemod;

    figure
    NFFT = 2^nextpow2(length(result));
    Y = abs(fft(result, NFFT)) / length(result);
    f = Fs / 2 * linspace(0, 1, NFFT / 2+1);
    plot(f, 2 * abs(Y(1:NFFT / 2 + 1)))
    xlim([0 200000])

    % Step 3.2: Apply FFT to dechirped signal
    received_fft(1, :) = fft(result);
end

% Helper function to convolve with a specified window
function output_signal = convolve_with_window(input_signal, window_type)
    % Define the window size based on the input signal length
    window_length = length(input_signal);
    
    % Generate the window based on the specified type
    switch window_type
        case 'rectangular'
            window = ones(1, window_length);
        case 'linear_increasing'
            window = linspace(0, 1, window_length);
        case 'linear_decreasing'
            window = linspace(1, 0, window_length);
        case 'triangular'
            window = 1 - abs(linspace(-1, 1, window_length));
        case 'welch'
            window = 1 - ((linspace(-1, 1, window_length)).^2);
        case 'sine'
            window = sin(pi * (linspace(0, 1, window_length) - 0.5));
        case 'hann'
            window = 0.5 * (1 - cos(2 * pi * (0:window_length - 1) / (window_length - 1)));
        case 'hamming'
            window = 0.54 - 0.46 * cos(2 * pi * (0:window_length - 1) / (window_length - 1));
        otherwise
            error('Unknown window type: %s', window_type);
    end
    
    % Convolve the input signal with the window
    output_signal = conv(input_signal, window, 'same'); % 'same' to keep the output size same as input, 'full' to show full output
end