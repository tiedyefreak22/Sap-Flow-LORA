clear all;
close all;
clc;
SF = 8; % Spreading factor
M = 2^SF; % no. of samples in one symbol
BW = 125e3; % Bandwidth
Fs = 10e6;  % sampling freq
Ts = 2^SF / BW;   % Symbol period
fc = 915e6;  % carrier center frequency
Power = 14;  % Tx power 14 dB
noise_sigma = 1;
numPreambleSymbols = 8; % Standard LoRa preamble symbol count
% message1 = "Hello World!";
% message2 = "Hello AlignTrack!";
% message3 = "Hello Goodbye!";
message1 = "1";
message2 = "2";
message3 = "3";
Fc = 921.5e6; % spectrum center frequency
df = Fc - fc;   % if CFO exists
shift = 0;
payload_offset1 = 25255;
payload_offset2 = 45255; % offset from second signal, not first

phy = LoRaPHY(fc, SF, BW, Fs);
phy.has_header = 1;                         % explicit header mode
phy.cr = 1;                                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                                % enable payload CRC checksum
phy.preamble_len = numPreambleSymbols;      % preamble: 8 basic upchirps

%% CSS modulated Signal
% Encode payload
fprintf("[encode] message 1:\n");
disp(double(char(message1)).')
symbols1 = phy.encode(double(char(message1)).');
fprintf("[encode] symbols 1:\n");
disp(symbols1);

% Encode payload2
fprintf("[encode] message 2:\n");
disp(double(char(message2)).')
symbols2 = phy.encode(double(char(message2)).');
fprintf("[encode] symbols 2:\n");
disp(symbols2);

% Encode payload3
fprintf("[encode] message 3:\n");
disp(double(char(message3)).')
symbols3 = phy.encode(double(char(message3)).');
fprintf("[encode] symbols 3:\n");
disp(symbols3);

% Baseband Modulation
switch find([numel(symbols1), numel(symbols2), numel(symbols3)] == max(max([numel(symbols1), numel(symbols2), numel(symbols3)])), 1, 'first')
    case 1
        signalIQ1 = zeros(numel(phy.modulate(symbols1)) + payload_offset1 + payload_offset2, 1);
        signalIQ2 = signalIQ1;
        signalIQ3 = signalIQ1;
    case 2
        signalIQ1 = zeros(numel(phy.modulate(symbols2)) + payload_offset1 + payload_offset2, 1);
        signalIQ2 = signalIQ1;
        signalIQ3 = signalIQ1;
    case 3
        signalIQ1 = zeros(numel(phy.modulate(symbols3)) + payload_offset1 + payload_offset2, 1);
        signalIQ2 = signalIQ1;
        signalIQ3 = signalIQ1;
    otherwise
        error("Can't determine max signal length");
end

signalIQ1(1:numel(phy.modulate(symbols1))) = phy.modulate(symbols1);
signalIQ2(payload_offset1 + 1:numel(phy.modulate(symbols2)) + payload_offset1) = phy.modulate(symbols2);
signalIQ3(payload_offset1 + payload_offset2 + 1:numel(phy.modulate(symbols3)) + payload_offset1 + payload_offset2) = phy.modulate(symbols3);
signalIQtotal = signalIQ1 + signalIQ2 + signalIQ3;

% spectrogram for one LoRa packet
figure(1)
spectrogram(signalIQtotal, 1000, 0, 1000, Fs, 'yaxis', 'centered')

% %AWGN noise
noise = noise_sigma * randn(length(signalIQtotal), 1);

% signal recieved at  LoRa gateway
%recieved_signal = signalIQtotal + noise;
recieved_signal = signalIQtotal;

% recieved signal after windowing and FFT
[recieved_fft] = LoRa_demod_1(recieved_signal, SF, BW, Fs, 0);
new_recieved_fft(1, :) = recieved_fft.';
[recieved_fft] = LoRa_demod_1(recieved_signal, SF, BW, Fs, payload_offset1);
new_recieved_fft(2, :) = recieved_fft.';
[recieved_fft] = LoRa_demod_1(recieved_signal, SF, BW, Fs, payload_offset2);
new_recieved_fft(3, :) = recieved_fft.';

% peak extraction algorith with align track decoding for complete packet
[row_fft, col_fft] = size(new_recieved_fft);
k = 6;
for fft_row_idx = 0:1:row_fft - 1
    for fft_col_idx = 1:1:col_fft
        [val(fft_row_idx + 1), idx(fft_row_idx + 1)] = max(new_recieved_fft(fft_row_idx + 1, :));
        r(fft_row_idx + 1) = mean(new_recieved_fft(fft_row_idx + 1, :)) + k * std(new_recieved_fft(fft_row_idx + 1, :));
        if val(fft_row_idx + 1) >= r(fft_row_idx + 1) % test against dynamic peak extraction threshold
            I(fft_row_idx + 1, fft_col_idx) = idx(fft_row_idx + 1);
            before = (idx(fft_row_idx + 1) - 1);
            after = (idx(fft_row_idx + 1) + 1);
            % finding local min before index value
            ii1 = 1;
            while ii1 < col_fft
                if before == 0
                    before = idx(fft_row_idx + 1);
                end
                if before-1 <= 0
                    before = idx(fft_row_idx + 1) + 1;
                end
                if recieved_fft(before + 1) > recieved_fft(before) && recieved_fft(before) < recieved_fft(before - 1)
                    b = before;
                    before = before-1;
                    break;
                else
                    before = before - 1;
                end
                ii1 = ii1 + 1;
            end
            % finding local min after index value
            ii2 = 1;
            while ii2 < col_fft
                if after >= col_fft
                    after = idx(fft_row_idx + 1);
                end
                if after + 1 >= col_fft
                    after = idx(fft_row_idx + 1) - 1;
                end
                if recieved_fft(after)<recieved_fft(after - 1) && recieved_fft(after + 1) > recieved_fft(after)
                    a = after;
                    after = after + 1;
                    break;
                else
                    after = after + 1;
                end
                ii2 = ii2 + 1;
            end
            %  remove local min points before and after index value
            for y = before:after
                new_recieved_fft(fft_row_idx + 1, y) = 0;
            end
        else
            break;
        end
    end
end
[m1 n1] = size(I);
% sort the index value for further processing
I = sort(I, 2);
% side lobe elimination
for ee = 1:1:m1
    for kk = 1:1:length(I(ee, :))
        if I(ee, kk) == 0
            continue;
        else
            d1 = kk;
            break;
        end
    end
    A.(sprintf('RandomVariable_%d', ee)) = I(ee, [d1:end]);
    AA = A.(sprintf('RandomVariable_%d', ee));
    issidelobe = zeros(1, length(AA));
    for i = 1:1:length(AA)
        for j = i + 1:1:length(AA)
            for k = 1:1:length(AA)
                if AA(j) - AA(i) == AA(i) - AA(k)
                    if k ~= j && recieved_fft(AA(k)) == recieved_fft(AA(j))
                        issidelobe(k) = 1;
                        issidelobe(j) = 1;
                    end
                end
            end
        end
        i
    end
    AA = AA(issidelobe ~= 1);
end