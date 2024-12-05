function [f, Y, Y2] = LoRa_demod_frft(signal, fc, SF, BW, Fs, cfo, window)
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

    % Generate Reference Downchirp
    phy = LoRaPHY(fc, SF, BW, Fs);
    dChirpsDemod  = phy.chirp(false, SF, BW, Fs, 0, cfo, 0);
    result = multiply_with_window(signal, window) .* dChirpsDemod;
    NFFT = 2^nextpow2(length(result));
    x2 = frft(result, NFFT, NFFT, pi/4);
    Y = abs(x2) / length(result);
    Y2 = imag(x2) / length(result);

    %NFFT=-T:0.01:T;
    figure,
    plot(NFFT,result);xlabel('Time');ylabel('Amplitude');title('Input Function');

    figure,
    plot(Y);
    xlabel('frequency');ylabel('Amplitude');title('Magnitude Response of FRFT');

    figure,
    plot(Y2);
    xlabel('frequency');ylabel('Amplitude');title('Imaginary part of FRFT');

    f = Fs / 2 * linspace(0, 1, NFFT / 2+1);
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

% Helper function to multiply with a specified window
function output_signal = multiply_with_window(input_signal, window_type)
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
    
    % Multiply the input signal with the window
    output_signal = input_signal .* window.';
end