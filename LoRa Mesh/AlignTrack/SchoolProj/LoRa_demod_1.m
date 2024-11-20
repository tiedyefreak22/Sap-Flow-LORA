function [received_fft] = LoRa_demod_1(signal, SF, BW, Fs, shift)
    % LoRa_Demodulate_Full demodulates full LoRa packet
    
    %% Return if SF is not in the range
    if SF > 12 || SF < 7
        return
    end
    %% Demodualte
    fc = 915e6;  % carrier center frequency
    cfo = 0;
    phy = LoRaPHY(fc, SF, BW, Fs);
    dChirpsDemod  = phy.chirp(false, SF, BW, Fs, 0, cfo, 0);
    len = length(dChirpsDemod);
    
    result = signal(shift + 1:shift + len) .* dChirpsDemod;
    %result = convolve_with_window(signal(shift + 1:shift + len), 'hamming') .* dChirpsDemod;

    figure
    NFFT = 2^nextpow2(length(result));
    Y = abs(fft(result, NFFT)) / length(result);
    f = Fs / 2 * linspace(0, 1, NFFT / 2+1);
    plot(f, 2 * abs(Y(1:NFFT / 2 + 1)))
    xlim([0 200000])
    received_fft = Y;
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