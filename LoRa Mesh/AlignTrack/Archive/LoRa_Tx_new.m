function [signal_mod] = LoRa_Tx_new(message,Bandwidth,SF,Pt,Fs,varargin)
% LoRa_Tx emulates a Lora transmission
%
%   in:  message      payload message
%        Bandwidth    signal bandwidth of LoRa transmisson  
%        SF           spreading factor
%        Pt           transmit power in deicbels
%        Fs           sampling frequency
%        dF           frequency offset
%        varargin{1}  code rate
%        varargin{2}  symbols in preamble
%        varargin{3}  sync key
%
%  out:  signal       LoRa IQ waveform

if nargin == 5
    CR = 1 ;
    n_preamble = 8 ;
    SyncKey = 5 ;
elseif nargin == 6
    CR = varargin{1} ;
    n_preamble = 8;
    SyncKey = 5 ;
elseif nargin == 7
    CR = varargin{1} ;
    n_preamble = varargin{2} ;
    SyncKey = 5 ;
elseif nargin == 8
    CR = varargin{1} ;
    n_preamble = varargin{2} ;
    SyncKey = varargin{3} ;
end
payload = LoRa_Encode_Full(message,SF,CR) ; % encode message
signal = LoRa_Modulate_Full(payload,SF,Bandwidth,n_preamble,SyncKey,Fs) ; % LoRa modulate message
signal_mod = 10.^(Pt./20).*signal ;
% signal_mod = 10.^(Pt./20).*signal.*exp(-j.*2.*pi.*df/Fs.*(0:length(signal)-1))' ; % use when CFO is present. frequency shift and convert to power
end
function [payload] = LoRa_Encode_Full(message,SF,CR)
% LoRa_Encode_Full emulates a Lora transmission
%
%   in:  message      payload message
%        SF           spreading factor
%        CR           coding rate 
%
%  out:  packet       encoded lora packet 
CRC_pld = 1 ;  % cyclic rate code flag
imp = 0 ;
opt = 0 ;
%% String to Decimal
message_chr = convertStringsToChars(message) ;
message_dbl = uint8(message_chr) ;
%% Packet Length Calculations
N_pld = (SF == 7).*1 + (SF == 8).*2 + (SF == 9).*3 + (SF == 10).*4 + (SF == 11).*5 + (SF == 12).*6 ;
n_packet = 8 + max([ceil((8*(length(message_dbl) + 5) - 4.*SF + 28 + 16.*CRC_pld - 20.*imp)/(4.*(SF - 2.*opt))).*(CR + 4) 0]) ;
n_wht = SF .* floor((n_packet-8)/(4 + CR)) + N_pld - 1 ;
n_pld = ceil((n_wht + (SF == 7).*0 + (SF == 8).*1 + (SF == 9).*2 + (SF == 10).*3 + (SF == 11).*4 + (SF == 12).*5)/2) ;
n_pad = n_pld - 5 - length(message_dbl) - CRC_pld.*2 ;
%% Create payload message
CRC_dbl = CRC_pld.*[1 1] ; % CRC is not working atm
pad_dbl = zeros(1,n_pad + N_pld - 1) ; % padding
payload = double([255 255 0 0 message_dbl 0 CRC_dbl pad_dbl]);  % LoRa payload
end
function [signal] = LoRa_Modulate_Full(payload,SF,Bandwidth,n_preamble,SyncKey,Fs)
% LoRa_Modulate_Full constructs a lora packet (preamble + sync header + payload)
%
%   in:  payload         payload 1xN symbol vector wher N=1-Inf 
%                       with values {0,1,2,...,2^(SF)-1}
%        SF             spreading factor   
%        Bandwidth      signal bandwidth of LoRa transmisson  
%        n_preamble     number of symbols in the preamble
%        SyncKey        synchronize key
%        Fs             sampling frequency
%
%  out:  signal          LoRa IQ packet
signal_prmb = loramod((SyncKey - 1).*ones(1,n_preamble),SF,Bandwidth,Fs,1) ; % preamble upchirps
signal_sync_u = loramod([0 0],SF,Bandwidth,Fs,1) ; % sync upchirp
signal_sync_d1 = loramod(0,SF,Bandwidth,Fs,-1) ; % header downchirp
signal_sync_d = [signal_sync_d1; signal_sync_d1; signal_sync_d1(1:length(signal_sync_d1)/4)] ; % concatenate header
signal_mesg = loramod(mod(payload + SyncKey,2^SF),SF,Bandwidth,Fs,1) ; % add sync key to payload messaage
signal = [signal_prmb; signal_sync_u; signal_sync_d; signal_mesg];  % concatenate LoRa packet
end
function [y] = loramod(x,SF,BW,fs,varargin)
% loramod LoRa modulates a symbol vector specified by x
%
%   in:  x          1xN symbol vector  
%                   with values {0,1,2,...,2^(SF)-1}
%        BW         signal bandwidth of LoRa transmisson  
%        SF         spreading factor   
%        Fs         sampling frequency
%        varargin{1} polarity of chirp
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
Ns      = fs.*M/BW ;
gamma   = x/Ts ;
beta    = BW/Ts ;
time    = (0:Ns-1)'.*1/fs ;
freq    = mod(gamma + Inv.*beta.*time,BW);
Theta   = cumtrapz(time,freq) ;
y       = reshape(exp(j.*2.*pi.*Theta),numel(Theta),1) ;
end