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
message1 = "Hello World!";
message2 = "Hello AlignTrack!";
message3 = "Hello Goodbye!";
Fc = 921.5e6; % spectrum center frequency
cfo = 0;
df = Fc - fc;   % if CFO exists

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

offset_time = round(max([(length(phy.modulate(symbols1)) / Fs) * 1000, (length(phy.modulate(symbols2)) / Fs) * 1000, (length(phy.modulate(symbols3)) / Fs) * 1000])); % ms

windows = {'rectangular', 'welch', 'sine', 'hann', 'hamming'};
tEnds = [];
trials = 20;

for trial = 1:1:trials
    trial
    
    while true
        try
            payload_offset1 = round(Fs * (randi([10, offset_time]) / 1000)); % samples
            payload_offset2 = round(Fs * (randi([10, offset_time]) / 1000)); % offset from second signal, not first

            % Baseband Modulation
            switch find([numel(symbols1), numel(symbols2), numel(symbols3)] == max(max([numel(symbols1), numel(symbols2), numel(symbols3)])), 1, 'first')
                case 1
                    diff = (max(max([numel(symbols1), numel(symbols2), numel(symbols3)])) - min(min([numel(symbols1), numel(symbols2), numel(symbols3)]))) * Ts * Fs;
                    signalIQ1 = zeros(numel(phy.modulate(symbols1)) + payload_offset1 + payload_offset2 - diff, 1);
                    signalIQ2 = signalIQ1;
                    signalIQ3 = signalIQ1;
                case 2
                    diff = (max(max([numel(symbols1), numel(symbols2), numel(symbols3)])) - min(min([numel(symbols1), numel(symbols2), numel(symbols3)]))) * Ts * Fs;
                    signalIQ2 = zeros(numel(phy.modulate(symbols2)) + payload_offset1 + payload_offset2 - diff, 1);
                    signalIQ1 = signalIQ2;
                    signalIQ3 = signalIQ2;
                case 3
                    diff = (max(max([numel(symbols1), numel(symbols2), numel(symbols3)])) - min(min([numel(symbols1), numel(symbols2), numel(symbols3)]))) * Ts * Fs;
                    signalIQ3 = zeros(numel(phy.modulate(symbols3)) + payload_offset1 + payload_offset2 - diff, 1);
                    signalIQ1 = signalIQ3;
                    signalIQ2 = signalIQ3;
                otherwise
                    error("Can't determine max signal length");
            end
            
            signalIQ1(1:numel(phy.modulate(symbols1))) = phy.modulate(symbols1);
            signalIQ2(payload_offset1 + 1:numel(phy.modulate(symbols2)) + payload_offset1) = phy.modulate(symbols2);
            signalIQ3(payload_offset1 + payload_offset2 + 1:numel(phy.modulate(symbols3)) + payload_offset1 + payload_offset2) = phy.modulate(symbols3);
            signalIQtotal = signalIQ1 + signalIQ2 + signalIQ3;
            break
        end
    end
    
    % Add zero leader
    initial_offset = Fs * (randi([1, 10]) / 1000);
    signalIQtotal = [zeros(initial_offset, 1); signalIQtotal];
    
    % spectrogram for one LoRa packet
    % figure
    % spectrogram(signalIQtotal, 1000, 0, 1000, Fs, 'yaxis', 'centered')
    % ylim([-BW / 2000000, BW / 2000000])
    
    % AWGN noise
    noise = noise_sigma * randn(length(signalIQtotal), 1);
    
    % signal recieved at  LoRa gateway
    received_signal = signalIQtotal + noise;
    %received_signal = signalIQtotal;

    for window = 1:1:length(windows)
        tStart = tic;
        
        % received signal after windowing and FFT
        [received_fft] = LoRa_demod_1(received_signal, fc, SF, BW, Fs, cfo, string(windows(window)));
        for i = 1:1:size(received_fft, 1)
            new_received_fft(i, :) = received_fft(i, :).';
        end
        
        % peak extraction algorithm with AlignTrack decoding for complete packet
        [row_fft, col_fft] = size(new_received_fft);
        k = 6;
        for fft_row_idx = 0:1:row_fft - 1
            for fft_col_idx = 1:1:col_fft
                [val(fft_row_idx + 1), idx(fft_row_idx + 1)] = max(new_received_fft(fft_row_idx + 1, :));
                r(fft_row_idx + 1) = mean(new_received_fft(fft_row_idx + 1, :)) + k * std(new_received_fft(fft_row_idx + 1, :));
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
                        if received_fft(before + 1) > received_fft(before) && received_fft(before) < received_fft(before - 1)
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
                        if received_fft(after) < received_fft(after - 1) && received_fft(after + 1) > received_fft(after)
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
                        new_received_fft(fft_row_idx + 1, y) = 0;
                    end
                else
                    break;
                end
            end
        end
        [m1, n1] = size(I);
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
            A.(sprintf('RandomVariable_%d', ee)) = I(ee, d1:end);
            AA = A.(sprintf('RandomVariable_%d', ee));
            issidelobe = zeros(1, length(AA));
            for i = 1:1:length(AA)
                for j = i + 1:1:length(AA)
                    for k = 1:1:length(AA)
                        if AA(j) - AA(i) == AA(i) - AA(k)
                            if k ~= j && received_fft(AA(k)) == received_fft(AA(j))
                                issidelobe(k) = 1;
                                issidelobe(j) = 1;
                            end
                        end
                    end
                end
            end
            AA = AA(issidelobe ~= 1);
        end
        tEnds(window, trial) = toc(tStart);
    end
end

windows
mean(tEnds, 2).'

figure
hold on
for i = 1:1:size(tEnds, 1)
    x = 1:1:trials;
    plot(x, tEnds(i, :))
    legend(windows)
end
hold off