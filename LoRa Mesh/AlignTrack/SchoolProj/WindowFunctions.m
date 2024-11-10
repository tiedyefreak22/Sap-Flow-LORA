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

% Plot results
figure;
subplot(4, 2, 1); plot(output_rectangular); title('Rectangular Window');
subplot(4, 2, 2); plot(output_linear_increasing); title('Linear Increasing Window');
subplot(4, 2, 3); plot(output_linear_decreasing); title('Linear Decreasing Window');
subplot(4, 2, 4); plot(output_triangular); title('Triangular Window');
subplot(4, 2, 5); plot(output_welch); title('Welch Window');
subplot(4, 2, 6); plot(output_sine); title('Sine Window');
subplot(4, 2, 7); plot(output_hann); title('Hann Window');
subplot(4, 2, 8); plot(output_hamming); title('Hamming Window');

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
            window = 0.5 * (1 - cos(2 * pi * (0:window_length-1) / (window_length-1)));
        case 'hamming'
            window = 0.54 - 0.46 * cos(2 * pi * (0:window_length-1) / (window_length-1));
        otherwise
            error('Unknown window type: %s', window_type);
    end
    
    % Convolve the input signal with the window
    output_signal = conv(input_signal, window, 'same'); % 'same' to keep the output size same as input, 'full' to show full output
end
