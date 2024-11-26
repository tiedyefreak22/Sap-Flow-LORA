% Apply different windowing functions to an input signal

% Example input signal (sample data)
%input_signal = sin(2 * pi * 0.1 * (0:99)); % A simple sine wave for demonstration
input_signal = ones(1,100); % Constant 1's vector

% Convolve input signal with various window types
output_rectangular = convolve_with_window(input_signal, 'rectangular');
output_linear_increasing = convolve_with_window(input_signal, 'linear_increasing');
output_linear_decreasing = convolve_with_window(input_signal, 'linear_decreasing');
output_triangular = convolve_with_window(input_signal, 'triangular');
output_welch = convolve_with_window(input_signal, 'welch');
output_sine = convolve_with_window(input_signal, 'sine');
output_hann = convolve_with_window(input_signal, 'hann');
output_hamming = convolve_with_window(input_signal, 'hamming');
output_parzen = convolve_with_window(input_signal, 'parzen');
output_blackman = convolve_with_window(input_signal, 'blackman');
output_nuttall = convolve_with_window(input_signal, 'nuttall');
output_kaiser = convolve_with_window(input_signal, 'kaiser');
output_chebyshev = convolve_with_window(input_signal, 'dolph_chebyshev');

% Plot results
figure;
subplot(5, 3, 1); plot(output_rectangular); title('Rectangular Window');
subplot(5, 3, 2); plot(output_linear_increasing); title('Linear Increasing Window');
subplot(5, 3, 3); plot(output_linear_decreasing); title('Linear Decreasing Window');
subplot(5, 3, 4); plot(output_triangular); title('Triangular Window');
subplot(5, 3, 5); plot(output_welch); title('Welch Window');
subplot(5, 3, 6); plot(output_sine); title('Sine Window');
subplot(5, 3, 7); plot(output_hann); title('Hann Window');
subplot(5, 3, 8); plot(output_hamming); title('Hamming Window');
subplot(5, 3, 9); plot(output_parzen); title('Parzen Window');
subplot(5, 3, 10); plot(output_blackman); title('Blackman Window');
subplot(5, 3, 11); plot(output_nuttall); title('Nuttall Window');
subplot(5, 3, 12); plot(output_kaiser); title('Kaiser Window');
subplot(5, 3, 13); plot(output_chebyshev); title('Dolph/Chebyshev Window');

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
            %window = triang(window_length)';
        case 'welch'
            window = 1 - ((linspace(-1, 1, window_length)).^2);
        case 'sine'
            window = sin(pi * (linspace(0, 1, window_length) - 0.5));
        case 'hann'
            window = 0.5 * (1 - cos(2 * pi * (0:window_length-1) / (window_length-1)));
            %window = hann(window_length)';
        case 'hamming'
            window = 0.54 - 0.46 * cos(2 * pi * (0:window_length-1) / (window_length-1));
            %window = hamming(window_length)';
        case 'parzen'
            window = parzenwin(window_length)';
        case 'blackman'
            window = blackman(window_length)';
        case 'nuttall'
            window = nuttallwin(window_length)';
        case 'kaiser'
            window = kaiser(window_length, 3)'; % Beta=3 is an example, can adjust as needed
        case 'dolph_chebyshev'
            window = chebwin(window_length, 60)'; % Ripple=60 dB is an example, can adjust
        otherwise
            error('Unknown window type: %s', window_type);
    end
    
    % Convolve the input signal with the window
    output_signal = conv(input_signal, window, 'same'); % 'same' to keep the output size same as input, 'full' to show full output
end
