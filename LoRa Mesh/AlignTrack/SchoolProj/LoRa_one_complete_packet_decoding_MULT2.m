clear all;
close all;
clc;
warning('off', 'all');

%% =============================================================================================
%% START Parameter Definition

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

%% END Parameter Definition
%% =============================================================================================
%% START Symbol Generation

% Encode payload 1
fprintf("[encode] message 1:\n");
disp(double(char(message1)).')
symbols1 = phy.encode(double(char(message1)).');
fprintf("[encode] symbols 1:\n");
disp(symbols1);

% Encode payload 2
fprintf("[encode] message 2:\n");
disp(double(char(message2)).')
symbols2 = phy.encode(double(char(message2)).');
fprintf("[encode] symbols 2:\n");
disp(symbols2);

% Encode payload 3
fprintf("[encode] message 3:\n");
disp(double(char(message3)).')
symbols3 = phy.encode(double(char(message3)).');
fprintf("[encode] symbols 3:\n");
disp(symbols3);

offset_time = round(min([(length(phy.modulate(symbols1)) / Fs) * 1000, (length(phy.modulate(symbols2)) / Fs) * 1000, (length(phy.modulate(symbols3)) / Fs) * 1000])); % ms

windows = {'rectangular', 'welch', 'sine', 'hann', 'hamming'};
tEnds = [];
SER = [];
cum_SER = [];
new_SER = [];
trials = 20;

for trial = 1:1:trials
    trial
    while true
        try
            payload_offset1 = round(Fs * (randi([round(offset_time / 2), offset_time]) / 1000)); % samples
            payload_offset2 = round(Fs * (randi([round(offset_time / 2), offset_time]) / 1000)); % offset from second signal, not first

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
    if trial == trials
        figure
        spectrogram(signalIQtotal, 1000, 0, 1000, Fs, 'yaxis', 'centered')
        ylim([-BW / 2000000, BW / 2000000])
        saveas(gcf,'LoRa_Spectrogram.png')
    end

    % AWGN noise
    noise = noise_sigma * randn(length(signalIQtotal), 1);

    % signal recieved at LoRa gateway
    received_signal = signalIQtotal + noise;

%% END Symbol Generation
%% =============================================================================================
%% START AlignTrack Trials

    for window = 1:1:length(windows)
        tStart = tic;

        % Cross-Correlation for Preamble Detection
        ref_upchirp = phy.chirp(true, SF, BW, Fs, 0, cfo, 0);
        chirp_len = length(ref_upchirp);
        k = BW / (chirp_len / Fs); % Chirp rate
        correlation = abs(xcorr(received_signal, ref_upchirp));
        [pks, locs] = findpeaks(correlation, 'MinPeakHeight', 0.5e4);
        i = 1;
        j = 1;
        shift(i) = locs(1) - length(received_signal);
        final_shift(j) = shift(i);
        while true
            try
                [decoded_message, ~] = phy.demodulate(received_signal(final_shift(j):end));
                decoded_messages = {decoded_message(:, 1).'};
                break
            catch
                while true
                    try
                        payload_offset1 = round(Fs * (randi([round(offset_time / 2), offset_time]) / 1000)); % samples
                        payload_offset2 = round(Fs * (randi([round(offset_time / 2), offset_time]) / 1000)); % offset from second signal, not first
            
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
                if trial == trials
                    figure
                    spectrogram(signalIQtotal, 1000, 0, 1000, Fs, 'yaxis', 'centered')
                    ylim([-BW / 2000000, BW / 2000000])
                    saveas(gcf,'LoRa_Spectrogram.png')
                end
            
                % AWGN noise
                noise = noise_sigma * randn(length(signalIQtotal), 1);
            
                % signal recieved at LoRa gateway
                received_signal = signalIQtotal + noise;
                tStart = tic;

                % Cross-Correlation for Preamble Detection
                ref_upchirp = phy.chirp(true, SF, BW, Fs, 0, cfo, 0);
                chirp_len = length(ref_upchirp);
                k = BW / (chirp_len / Fs); % Chirp rate
                correlation = abs(xcorr(received_signal, ref_upchirp));
                [pks, locs] = findpeaks(correlation, 'MinPeakHeight', 0.5e4);
                i = 1;
                j = 1;
                shift(i) = locs(1) - length(received_signal);
                final_shift(j) = shift(i);
            end
        end

        % Loop over received signal, run through AlignTrack algorithm, then use
        % output to advance to next peak corresponding to next window alignment
        % point
        while (shift(i) + chirp_len - 1) <= length(received_signal)
            % received signal after windowing and FFT
            [f, received_fft] = LoRa_demod_1(received_signal(shift(i):shift(i) + chirp_len - 1), fc, SF, BW, Fs, cfo, string(windows(window)));
            AA = AlignTrack(received_fft);
            
            if length(AA) > 1
                [B, IX] = sort(AA); %order the amplitudes
                % A1 = B(end); %amplitude of second peak
                % A2 = B(end - 1); %amplitude of second peak
                % f1 = f(IX(end)); %frequency of first peak
                f2 = f(IX(end - 1)); %frequency of second peak
            else
                break
            end

            i = i + 1;
            shift(i) = shift(i - 1) + abs(Fs * ((f2 - BW) / k));
            [decoded_message, ~] = phy.demodulate(received_signal(shift(i):end));

            % Append the new row if it's not a duplicate
            if (~isempty(decoded_message)) && (~any(cellfun(@(row) isequal(row, decoded_message(:, 1).'), decoded_messages)))
                decoded_messages{end + 1} = decoded_message(:, 1).';
                final_shift(end + 1) = shift(i);
            end
        end

        final_shift

        % Loop through results, cross-correlate with the inputs to determine
        % closest match, then compare to determine SER
        symbol_corr = [];
        for i = 1:length(decoded_messages)
            fprintf('Row %d: ', i);
            disp(decoded_messages{i}); % Display the contents of the cell
            max_num = 0;

            % Loop through each output row
            for j = 1:3
                switch j
                    case 1
                        symbols = symbols1;
                    case 2
                        symbols = symbols2;
                    case 3
                        symbols = symbols3;
                end
                correlation = abs(xcorr(symbols, decoded_messages{i}));

                if max(correlation) > max_num
                    max_num = max(correlation);
                    symbol_corr(i) = j;
                end
            end

            switch symbol_corr(i)
                case 1
                    symbols = symbols1;
                case 2
                    symbols = symbols2;
                case 3
                    symbols = symbols3;
            end

            % Ensure arrays are the same length
            if length(symbols.') ~= length(decoded_messages{i})
                error('Input and output arrays must have the same length.');
            end

            % Calculate the number of symbol errors
            numErrors = sum(symbols.' ~= decoded_messages{i});

            % Calculate the symbol error rate (SER)
            SER(window, trial, i) = numErrors / length(symbols.');

            % Display results
            fprintf('Correlated Input: %d\n', symbol_corr(i));
            fprintf('Number of symbol errors: %d\n', numErrors);
            fprintf('Symbol Error Rate (SER): %.4f\n', SER(window, trial, i));

        end
        % Calculate the cumulative SER
        if trial == 1
            cum_SER(window, trial) = sum(SER(window, trial, :));
        else
            cum_SER(window, trial) = cum_SER(window, trial - 1) + sum(SER(window, trial, :));
        end
        tEnds(window, trial) = toc(tStart);
    end
end

%% END AlignTrack Trials
%% =============================================================================================
%% START Plot/Display Results

windows
mean(tEnds, 2).'

% Plot time for each window by trial
figure
hold on
title("Time vs. Window by Trial")
for i = 1:1:size(tEnds, 1)
    x = 1:1:trials;
    plot(x, tEnds(i, :))
    legend(windows)
end
saveas(gcf,'Time_vs_Window_by_Trial.png')
hold off

% Plot SER for each window by trial
figure
hold on
title("SER vs. Window by Trial")
for i = 1:1:size(SER, 1)
    x = 1:1:trials;
    plot(x, sum(SER(i, :, :), 3))
    legend(windows)
end
saveas(gcf,'SER_vs_Window_by_Trial.png')
hold off

% Plot cumulative SER for each window
figure
hold on
title("Cumulative SER vs. Window")
for i = 1:1:size(cum_SER, 1)
    x = 1:1:trials;
    plot(x, cum_SER(i, :))
    legend(windows)
end
saveas(gcf,'Cumulative_SER_vs_Window.png')
hold off

%% END Plot/Display Results
%% =============================================================================================