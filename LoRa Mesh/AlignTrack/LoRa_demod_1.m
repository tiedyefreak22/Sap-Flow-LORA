function [received_fft] = LoRa_demod_1(signal,SF,BW,Fs,shift)
% LoRa_Demodulate_Full demodulates full LoRa packet

%% Return if SF is not in the range
if SF > 12 || SF < 7
    return
end
%% Demodulate
dChirpsDemod  = loramod(0,SF,BW,Fs,-1);
len=length(dChirpsDemod);
result=signal(shift+1:shift+len).*dChirpsDemod;
fft_signal = (fft(result));% take fft window
received_fft=abs(fft_signal);
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