function [recieved_fft] = LoRa_demod_1(signal,SF,BW,Fs,shift)
% LoRa_Demodulate_Full demodulates full LoRa packet

%% Return if SF is not in the range
if SF > 12 || SF < 7
    return
end
%% Demodualte
dChirpsDemod  = loramod(0,SF,BW,Fs,-1);
len=length(dChirpsDemod);

%result=signal(shift+1:shift+len).*dChirpsDemod;
result=convolve_with_window(signal(shift+1:shift+len), 'rectangular').*dChirpsDemod;

fft_signal = (fft(result));% take  fft window
recieved_fft=abs(fft_signal);
end
function [y] = loramod(x,SF,BW,Fs,varargin)
% loramod LoRa modulates a symbol vector specified by x
%
%   in:  x          1xN symbol vector
%                   with values {0,1,2,...,2^(SF)-1}
%        BW         signal bandwidth of LoRa transmisson
%        SF         spreading factor
%        Fs         sampling frequency
%        varargin{1} set polarity of chirp
%
%  out:  y          LoRa IQ waveform
if (nargin < 4)
    error(message('comm:pskmod:numarg1'));
end
if (nargin > 5)
    error(message('comm:pskmod:numarg2'));
end
% Check that x is a positive integer
if (~isreal(x) || any(any(ceil(x) ~= x)) || ~isnumeric(x))
    error(message('comm:pskmod:xreal1'));
end
M       = 2^SF ;
% Check that M is a positive integer
if (~isreal(M) || ~isscalar(M) || M<=0 || (ceil(M)~=M) || ~isnumeric(M))
    error(message('comm:pskmod:Mreal'));
end
% Check that x is within range
if ((min(min(x)) < 0) || (max(max(x)) > (M-1)))
    error(message('comm:pskmod:xreal2'));
end
% Polarity of Chirp
if nargin == 4
    Inv = 1 ;
elseif nargin == 5
    Inv = varargin{1} ;
end
% Symbol Constants
Ts      = 2^SF/BW ;
Ns      = Fs.*M/BW ;
gamma   = x/Ts ;
beta    = BW/Ts ;
time    = (0:Ns-1)'.*1/Fs ;
freq    = mod(gamma + beta.*time,BW) ;
Theta   = cumtrapz(time,Inv.*freq) ;
y       = reshape(exp(j.*2.*pi.*Theta),numel(Theta),1) ;
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
            window = 0.5 * (1 - cos(2 * pi * (0:window_length-1) / (window_length-1)));
        case 'hamming'
            window = 0.54 - 0.46 * cos(2 * pi * (0:window_length-1) / (window_length-1));
        otherwise
            error('Unknown window type: %s', window_type);
    end
    
    % Convolve the input signal with the window
    output_signal = conv(input_signal, window, 'same'); % 'same' to keep the output size same as input, 'full' to show full output
end
